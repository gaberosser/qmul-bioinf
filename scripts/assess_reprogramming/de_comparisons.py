from rnaseq import loader, differential_expression, filter, general
from plotting import common, venn, bar
import pandas as pd
import numpy as np
import os
import re
import collections
from utils import output, log, setops, excel
from utils.string_manipulation import make_tuple
from matplotlib import pyplot as plt
import seaborn as sns
from copy import copy

logger = log.get_console_logger()


class SetMe(object):
    """
    This is useful to determine which fields need to be set
    """


def de_grouped_dispersion(dat, groups, comparisons, min_cpm=1., **de_params):
    cpm = dat.divide(dat.sum(), axis=1) * 1e6
    keep = (cpm > min_cpm).sum(axis=1) > 0
    dat = dat.loc[keep]

    res = differential_expression.run_multiple_de(
        the_data=dat,
        the_groups=groups,
        comparisons=comparisons,
        **de_params
    )

    for the_comparison, this_res in res.items():
        try:
            the_genes = this_res.index
            the_groups = groups[groups.isin(the_comparison)]
            the_cpm = cpm.loc[the_genes, the_groups.index]
            keep = filter.filter_cpm_by_group(
                the_cpm,
                groups,
                min_cpm=min_cpm
            )
            res[the_comparison] = this_res.loc[keep]
        except Exception as exc:
            print repr(exc)

    return res


def load_refs(ref_dict, **load_kwds):
    ref_objs_arr = []

    for k, v in ref_dict.items():
        the_kwds = copy(load_kwds)
        for k1, v1 in the_kwds.items():
            if v1 is SetMe:
                the_kwds[k1] = v.get(k1)
        the_obj = loader.load_references(k, **the_kwds)
        the_obj.batch_id = v['batch']
        ref_objs_arr.append(the_obj)

    if len(ref_objs_arr) == 1:
        return ref_objs_arr[0]
    else:
        return loader.loader.MultipleBatchLoader(ref_objs_arr)


def core_de(ipsc_vs_esc):
    """
    For each result, get the set of core DE genes that are present in all ref comparisons
    :return:
    """
    keep_cols = ['logFC', 'FDR', 'Direction']
    refs = set()
    core_tbl = {}
    core_ix = {}
    for k, v in ipsc_vs_esc.items():
        if k[0] not in core_ix:
            core_ix[k[0]] = v.index
        else:
            core_ix[k[0]] = v.index.intersection(core_ix[k[0]])
        refs.add(k[1])

    refs = sorted(refs)

    for k, v in ipsc_vs_esc.items():
        if k[0] not in core_tbl:
            core_tbl[k[0]] = pd.DataFrame(index=core_ix[k[0]])
        for t in keep_cols:
            core_tbl[k[0]].insert(0, "%s_%s" % (t, k[1]), v[t])

    # run back through and check that direction is consistent
    for k in core_tbl:
        the_cols = ["logFC_%s" % r for r in refs]
        eq = np.sign(core_tbl[k].loc[:, the_cols]).apply(lambda x: all(x == x[0]), axis=1)
        core_tbl[k] = core_tbl[k].loc[eq]

    return core_tbl


def classify_de_residual_denovo(core_tbl, parental, exclude=None):
    """
    :param core_tbl: From core_de() function (e.g. iPSC vs ESC).
    :param parental: Raw results from parental comparison (e.g. iPSC vs FB).
    :param exclude: If present, this contains features (genes) to exclude from the analysis.
    :return:
    """
    res = {}

    for k1, this_tbl in core_tbl.items():
        if k1 not in parental:
            logger.warning("Sample %s has no matching parental sample. Skipping", k1)
            continue

        if exclude is not None:
            this_tbl = this_tbl.loc[this_tbl.index.difference(exclude)]
        mean_logfc = this_tbl.loc[:, this_tbl.columns.str.contains('logFC_')].mean(axis=1)

        # look for the same DMRs in the parental comparison
        this_par = parental[k1]
        p_ix = this_par.index.intersection(this_tbl.index)
        mlfc_lookup = mean_logfc.loc[p_ix]

        res[k1] = this_tbl.copy()
        clas = pd.Series('partial', index=this_tbl.index)

        # de novo: DM in (A vs B) AND (A vs C) with same direction
        ix = (mlfc_lookup > 0.) & (this_par.loc[p_ix, 'logFC'] > 0.)
        clas.loc[ix.index[ix]] = 'up_de_novo'
        ix = (mlfc_lookup < 0.) & (this_par.loc[p_ix, 'logFC'] < 0.)
        clas.loc[ix.index[ix]] = 'down_de_novo'

        # residual: DM in (A vs B) AND NOT in (A vs C)
        not_p_ix = ~this_tbl.index.isin(this_par.index)
        mlfc_lookup = mean_logfc.loc[not_p_ix]
        ix = (mlfc_lookup > 0)
        clas.loc[ix.index[ix]] = 'up_residual'
        ix = (mlfc_lookup < 0)
        clas.loc[ix.index[ix]] = 'down_residual'

        res[k1].insert(2, 'classification', clas)

    return res


def get_dm_associated_de(
        dmr_ids,
        de_res_full,
        dmr_res_full,
        dmr_id_to_ens,
        ens_to_dmr_id,
        ref_labels
):
    all_genes = setops.reduce_union(*[dmr_id_to_ens[i].values for i in dmr_ids])
    this_de_full = {}
    for r in ref_labels:
        tt = de_res_full[r]
        tt = tt.loc[tt['Gene Symbol'].isin(all_genes)]
        this_de_full[r] = tt

    common_ix = sorted(setops.reduce_intersection(*[t.index for t in this_de_full.values()]))
    common_gs = this_de_full.values()[0].loc[common_ix, 'Gene Symbol']

    dmr_median_delta = {}
    for e in common_ix:
        dmr_median_delta[e] = {}
        for r in ref_labels:
            dmr_median_delta[e][r] = np.mean([dmr_res_full[r][i]['median_change'] for i in ens_to_dmr_id[e]])
    dmr_median_delta = pd.DataFrame(dmr_median_delta).transpose().sort_index()

    this_logfc = pd.concat(
        (this_de_full[r].loc[common_ix, 'logFC'] for r in ref_labels),
        axis=1
    )
    this_logfc.columns = ref_labels
    this_logfc.index = common_gs

    de_logfc = this_logfc.sort_index()

    return {'de': de_logfc, 'dmr': dmr_median_delta}


if __name__ == '__main__':
    pids = ['017', '018', '019', '030', '031', '049', '050', '054', '061', '026', '052']
    min_cpm = 1.
    n_above_min = 3
    eps = 0.01  # offset to use when applying log transform
    n_hipsci = 12

    dmr_indir = os.path.join(output.OUTPUT_DIR, "assess_reprogramming_dmr")

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'QLGLM'
    }

    ref_dict = {
        'encode_roadmap/ENCSR000EYP': {'batch': 'ENCODE Wold', 'strandedness': 'u'},
        # 'encode_roadmap/ENCSR670WQY': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        'GSE62772': {'batch': 'Cacchiarelli et al.', 'strandedness': 'u'},
        'GSE97265': {'batch': 'Kogut et al.', 'strandedness': 'r'},
    }

    # discard irrelevant samples
    to_discard = [
        '_rapamycin_',
        'HFF_PD16', 'HFF_PD46', 'HFF_PD64', 'HFF_PD74',
        'BJ_OLD',
        'IMR90_O',
        'MRC_5_PD52', 'MRC_5_PD62', 'MRC_5_PD72',
        'WI_38_O',
        '-EPS-',
        'F50S', 'I50S', 'FN1',
        'DOX', 'hiF_', 'hiF-T_', 'BJ_P12', '_MTG_', '_VTN_', 'hIPSC_P15_Rep1'
    ]

    to_aggr = [
        ('IN2-', 'IN2'),
        ('I50-', 'I50'),
    ]

    nsc_ref_dict = {
        # 'E-MTAB-3867': {'batch': 'Caren et al.', 'strandedness': '?'},
        'GSE61794': {'batch': 'Duan et al.', 'strandedness': '?'},
        # 'GSE64882': {'batch': 'Shahbazi et al.', 'strandedness': '?'},
        # 'encode_roadmap/ENCSR244ISQ': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR291IZK': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR572EET': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR977XUX': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
    }

    # to_aggr_nsc = [
    #     (r'H9_NSC_[12]', 'H9 NSC'),
    #     # (r'Pollard NSC [12]', 'Fetal NSC'),
    # ]

    outdir = output.unique_output_dir("assess_reprogramming_de")

    load_kwds = {
        'source': 'star',
        'alignment_subdir': SetMe,
        'strandedness': SetMe,
    }

    # our data (everything)
    obj = loader.load_by_patient(pids, source='star')

    # HipSci data
    hip_obj = loader.hipsci_ipsc(aggregate_to_gene=True)
    hip_obj.meta.insert(3, 'batch', hip_obj.batch_id)
    # hip_obj.meta.insert(3, 'batch', 'HipSci')

    # reduce the number in a (repeatably) random fashion
    rs = np.random.RandomState(42)  # set the seed so we always get the same samples
    keep = np.zeros(hip_obj.meta.shape[0]).astype(bool)
    idx = hip_obj.meta.index.tolist()
    rs.shuffle(idx)
    idx = idx[:n_hipsci]
    hip_obj.meta = hip_obj.meta.loc[idx]
    hip_obj.data = hip_obj.data.loc[:, idx]

    # References without NSC
    ref_obj = load_refs(ref_dict, **load_kwds)

    # discard irrelevant samples
    for td in to_discard:
        the_idx = ~ref_obj.data.columns.str.contains(td)
        ref_obj.filter_samples(the_idx)

    for srch, repl in to_aggr:
        ref_obj.aggregate_by_pattern(srch, repl)

    # fill in missing cell types
    ref_obj.meta.loc[ref_obj.meta.type.isnull(), 'type'] = ref_obj.meta.loc[ref_obj.meta.type.isnull(), 'cell type']
    ref_obj.rename_with_attributes(existing_attr='batch')

    # assemble the data we need for comparison:
    # - ESC (2 studies)
    # - iPSC (ours, GSE97265)
    # - FB (ours, GSE97265)
    ref_labels = ['cacchiarelli', 'encode']

    ix = obj.meta.type.isin(['iPSC', 'FB'])
    dat1 = obj.data.loc[:, ix]

    ix = ref_obj.meta.batch.str.contains('Kogut') & ref_obj.meta.type.isin(['iPSC', 'FB'])
    dat2 = ref_obj.data.loc[:, ix]

    ix = ref_obj.meta.index.str.contains('Cacchiarelli') & (ref_obj.meta.type == 'ESC')
    dat_esc_1 = ref_obj.data.loc[:, ix]

    ix = ref_obj.meta.index.str.contains('ENCODE') & (ref_obj.meta.type == 'ESC')
    dat_esc_2 = ref_obj.data.loc[:, ix]

    the_dat = pd.concat((dat1, dat2, dat_esc_1, dat_esc_2), axis=1)
    the_groups = pd.Series(index=the_dat.columns)
    t = ("ESC_%s" % ref_labels[0], "ESC_%s" % ref_labels[1])
    the_comparisons = {
        t: "%s - %s" % t
    }

    for pid in pids:
        k_fb = 'FB_%s_ours' % pid
        k_ipsc = 'iPSC_%s_ours' % pid
        ix_fb = the_groups.index.str.contains('DURA%s_FB' % pid)
        ix_ipsc = the_groups.index.str.contains('DURA%s_IPSC' % pid)
        the_groups.loc[ix_fb] = k_fb
        the_groups.loc[ix_ipsc] = k_ipsc
        if (ix_fb.sum() > 0) & (ix_ipsc.sum() > 0):
            the_comparisons[(k_ipsc, k_fb)] = "%s - %s" % (k_ipsc, k_fb)
            for r in ref_labels:
                the_comparisons[(k_ipsc, 'ESC_%s' % r)] = "%s - ESC_%s" % (k_ipsc, r)

    for s in ['N2', '50']:
        k_fb = 'FB_%s_kogut' % s
        k_ipsc = 'iPSC_%s_kogut' % s
        ix_fb = the_groups.index.str.contains('F%s' % s) & the_groups.index.str.contains('Kogut')
        ix_ipsc = the_groups.index.str.contains('I%s' % s) & the_groups.index.str.contains('Kogut')
        the_groups.loc[ix_fb] = k_fb
        the_groups.loc[ix_ipsc] = k_ipsc
        if (ix_fb.sum() > 0) & (ix_ipsc.sum() > 0):
            the_comparisons[(k_ipsc, k_fb)] = "%s - %s" % (k_ipsc, k_fb)
            for r in ref_labels:
                the_comparisons[(k_ipsc, 'ESC_%s' % r)] = "%s - ESC_%s" % (k_ipsc, r)

    the_groups.loc[the_groups.index.str.contains('Cacchiarelli')] = 'ESC_cacchiarelli'
    the_groups.loc[the_groups.index.str.contains('ENCODE')] = 'ESC_encode'

    de_res_full = de_grouped_dispersion(
        the_dat,
        the_groups,
        the_comparisons,
        min_cpm=min_cpm,
        return_full=True,
        **de_params
    )
    # rename the keys to simplify
    # de_res_full_s1 = dict([(pid, de_res_full[("GBM%s" % pid, "iNSC%s" % pid)]) for pid in pids])

    # extract only significant DE genes
    de_res_sign = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full.items()])

    # report baseline variation between ESC references
    t = ("ESC_%s" % ref_labels[0], "ESC_%s" % ref_labels[1])
    print "Between the two ESC references, there are %d DE genes. %d up, %d down." % (
        len(de_res_sign[t]),
        (de_res_sign[t].Direction == 'down').sum(),
        (de_res_sign[t].Direction == 'up').sum(),
    )

    # Venn diagrams (two ESC studies)
    s1 = ['_'.join(t) for t in zip(pids, ['ours'] * len(pids))]
    s2 = ['_'.join(t) for t in zip(['N2', '50'], ['kogut'] * 2)]
    fig, axs = plt.subplots(ncols=3, nrows=3)
    i = 0
    for s in s1 + s2:
        this_arr = []
        for r in ref_labels:
            k = ('iPSC_%s' % s, 'ESC_%s' % r)
            if k in de_res_sign:
                this_arr.append(de_res_sign[k])
        if len(this_arr):
            print "Found comparison %s" % s
            venn.venn_diagram(*[t.index for t in this_arr], ax=axs.flat[i], set_labels=ref_labels)
            axs.flat[i].set_title(s)
            i += 1
        else:
            print "No comparison %s" % s

    for j in range(i, axs.size):
        axs.flat[i].axis('off')

    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.05, top=0.95)
    fig.savefig(os.path.join(outdir, "ipsc_esc_venn.png"), dpi=200)

    # core DE genes
    ipsc_vs_esc = dict([
        (k, v) for k, v in de_res_sign.items() if re.search('|'.join(ref_labels), k[1])
                                                  and not re.search('|'.join(ref_labels), k[0])
    ])
    ipsc_vs_esc_core = core_de(ipsc_vs_esc)

    # can check the number of results with opposing signs using this (and the Venn figure):
    # [(k, len(v)) for k, v in ipsc_vs_esc_core.items()]

    # matched FB results
    ipsc_vs_fb_matched = dict([
        (k[0], v) for k, v in de_res_sign.items() if k[1][:2] == 'FB'
    ])

    # classify (and add gene symbols back in)
    ipsc_esc_fb = classify_de_residual_denovo(ipsc_vs_esc_core, ipsc_vs_fb_matched)
    for v in ipsc_esc_fb.values():
        general.add_gene_symbols_to_ensembl_data(v)
    core_de_classified_count = pd.DataFrame(dict(
        [(k, v.classification.value_counts()) for k, v in ipsc_esc_fb.items()]
    ))

    # residual / de novo bar chart
    # combine these results to generate bar charts
    set_colours_hypo = ['#b2df8a', '#33a02c']
    set_colours_hyper = ['#fb9a99', '#e31a1c']
    plot_colours = {'down': set_colours_hypo[::-1], 'up': set_colours_hyper[::-1]}

    df = core_de_classified_count.transpose()

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.5, 5.))
    for i, typ in enumerate(['up', 'down']):
        ax = axs[i]
        this = df.loc[:, df.columns.str.contains(typ)].transpose()
        this.index = ['De novo', 'Residual']
        new_cols = this.columns.str.replace('_ours', '').str.replace('_', '').str.replace('kogut', ' (Kogut et al.)')
        this.columns = new_cols
        bar.stacked_bar_chart(this, colours=plot_colours[typ], ax=ax, width=0.8, ec='k', lw=1.0)
        # df.loc[:, df.columns.str.contains(typ)].plot.bar(stacked=True, colors=plot_colours[typ], ax=ax, width=0.9)
        ax.set_ylabel('%s DE genes' % ('Downregulated' if typ == 'down' else 'Upregulated'))
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
        ax.get_legend().set_frame_on(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "number_de_residual_denovo.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "number_de_residual_denovo.tiff"), dpi=200)

    # check in DMRs
    dmrs_classified = {}
    dedmr_results = {}
    both_genes = {}

    for pid in pids:
        fn = os.path.join(dmr_indir, "iPSC%s_classified_dmrs.csv" % pid)
        if os.path.exists(fn):
            this_dmr = pd.read_csv(fn, header=0, index_col=0)
            this_dmr.loc[:, 'genes'] = this_dmr.genes.apply(make_tuple)
            dmrs_classified[pid] = this_dmr
            if "iPSC_%s_ours" % pid in ipsc_esc_fb:
                this_de_res = ipsc_esc_fb["iPSC_%s_ours" % pid]

                de_genes = this_de_res.loc[:, 'Gene Symbol'].dropna()
                dmr_genes = this_dmr.loc[:, 'genes'].values
                dmr_genes = setops.reduce_union(*dmr_genes) if len(dmr_genes) else []
                both_genes[pid] = set(de_genes).intersection(dmr_genes)

                if len(both_genes[pid]):
                    # DE
                    the_de_res = this_de_res.loc[this_de_res['Gene Symbol'].isin(both_genes[pid])]
                    the_de_res = the_de_res.loc[:, ~the_de_res.columns.str.contains('Direction')]
                    the_de_res.set_index('Gene Symbol', inplace=True)  # TODO: may break if there are duplicates?

                    # DMR
                    the_dmr_res = this_dmr.loc[this_dmr.genes.astype(str).str.contains('|'.join(both_genes[pid]))]
                    # 'unfold' by gene, group by gene, take mean
                    ix = reduce(lambda x, y: x + y, [[(i, t) for t in row.genes] for i, row in the_dmr_res.iterrows()])
                    the_dmr_res_unfold = {}
                    for i, row in the_dmr_res.iterrows():
                        for g in row.genes:

                            the_dmr_res_unfold["%d_%s" % (i, g)] = row.copy()
                            the_dmr_res_unfold["%d_%s" % (i, g)].loc['G'] = g
                    the_dmr_res_unfold = pd.DataFrame(the_dmr_res_unfold).transpose()
                    # reduce to numeric types and aggregate by gene
                    the_dmr_res_unfold_num = the_dmr_res_unfold.loc[
                                             :,
                                             the_dmr_res_unfold.columns.str.contains(r'(median_delta|padj)')
                                             ].astype(float)
                    the_dmr_res_by_gene = the_dmr_res_unfold_num.groupby(the_dmr_res_unfold.G).mean()
                    # combine classifications: keep consistent classifications, otherwise drop them
                    gb = the_dmr_res_unfold.classification.groupby(the_dmr_res_unfold.G)
                    the_dmr_res_class_by_gene = gb.apply(lambda x: x.unique()[0] if len(x.unique()) == 1 else "N/A")
                    the_dmr_res_by_gene.insert(0, 'classification', the_dmr_res_class_by_gene)

                    # combine results
                    dmr_comparisons = [('encode', 'ESC line 1'), ('e6194', 'ESC line 2')]
                    de_comparisons = [('encode', 'ESC line 1'), ('cacchiarelli', 'ESC line 2')]
                    this = {}
                    for g, de_row in the_de_res.iterrows():
                        dm_row = the_dmr_res_by_gene.loc[g]
                        the_row = [de_row['classification']] + \
                                  [de_row['FDR_ESC_%s' % t[0]] for t in de_comparisons] + \
                                  [de_row['logFC_ESC_%s' % t[0]] for t in de_comparisons]
                        the_row += [dm_row['classification']] + \
                                   [dm_row['padj_%s' % t[0]] for t in dmr_comparisons] + \
                                   [dm_row['median_delta_%s' % t[0]] for t in dmr_comparisons]
                        this[g] = the_row
                    dedmr_results[pid] = pd.DataFrame(this).transpose()
                    dedmr_results[pid].columns = pd.MultiIndex.from_arrays([
                        ['DE'] * 5 + ['DMR'] * 5,
                        ['Classification', 'FDR ESC line 1', 'FDR ESC line 2', 'logFC ESC line 1', 'logFC ESC line 2'] \
                        + ['Classification', 'FDR ESC line 1', 'FDR ESC line 2', 'Median delta M ESC line 1', 'Median delta M ESC line 2']
                    ])

                print "Patient %s: %d in DE and DMR." % (pid, len(both_genes[pid]))

    # report: how many genes are DE and DMR across all lines?
    logger.info(
        "In total, %d unique genes are DE and DMR across all %d lines. "
        "Of these, %d are shared across all lines.",
        len(setops.reduce_union(*both_genes.values())),
        len(both_genes),
        len(setops.reduce_intersection(*both_genes.values()))
    )

    # output these in tabular form
    # add a sheet of explanation
    to_export = collections.OrderedDict()
    to_export['Explanation'] = pd.Series(collections.OrderedDict(
        [
            ('DE ESC line 1', 'ENCODE H1 ESC (%d replicates)' % the_groups.value_counts()['ESC_encode']),
            ('DE ESC line 2', 'Cacchiarelli et al. (%d lines)' % the_groups.value_counts()['ESC_cacchiarelli']),
            ('DMR ESC line 1', 'ENCODE H7 ESC (no replicates)'),
            ('DMR ESC line 2', 'Weltner et al. H9 (3 replicates)'),
        ],
    ), name='All comparisons are stated as iPSC - ESC.')
    to_export.update(collections.OrderedDict([(pid, dedmr_results[pid]) for pid in pids if pid in dedmr_results]))

    excel.pandas_to_excel(to_export, os.path.join(outdir, "DE_DMR_results_combined.xlsx"))


    # more nuanced approach: check the correlation between DMRs and corresponding genes at the level of median delta M
    # and logFC.

    de_logfc = {}
    dmr_median_delta = {}

    for pid in pids:
        fn = os.path.join(dmr_indir, "iPSC%s_classified_dmrs.csv" % pid)
        if (pid in dmrs_classified) and ("iPSC_%s_ours" % pid in ipsc_esc_fb):
            this_dmr = dmrs_classified[pid]
            dmr_genes = this_dmr.loc[:, 'genes'].values
            dmr_genes = setops.reduce_union(*dmr_genes) if len(dmr_genes) else []
            this_de_sign = ipsc_esc_fb["iPSC_%s_ours" % pid]
            this_de_full = {}
            for r in ref_labels:
                tt = de_res_full[("iPSC_%s_ours" % pid, "ESC_%s" % r)]
                tt = tt.loc[tt['Gene Symbol'].isin(dmr_genes)]
                this_de_full[r] = tt

            common_ix = sorted(setops.reduce_intersection(*[t.index for t in this_de_full.values()]))
            common_gs = this_de_full.values()[0].loc[common_ix, 'Gene Symbol']

            this_dmr_dm = {}

            for g in common_gs:
                this_ix = this_dmr.index[this_dmr.genes.apply(lambda x: g in x)]
                this_dmr_dm[g] = this_dmr.loc[this_ix, this_dmr.columns.str.contains('median_delta')].mean(axis=0)
            dmr_median_delta[pid] = pd.DataFrame(this_dmr_dm).transpose().sort_index()

            this_logfc = pd.concat(
                (this_de_full[r].loc[common_ix, 'logFC'] for r in ref_labels),
                axis=1
            )
            this_logfc.columns = ref_labels
            this_logfc.index = common_gs

            de_logfc[pid] = this_logfc.sort_index()

    # scatterplot (one per patient)
    import statsmodels.formula.api as sm

    for pid in pids:
        if pid in de_logfc:
            fig = plt.figure(figsize=(5, 4))
            ax = fig.add_subplot(111)
            this_de = de_logfc[pid].mean(axis=1)
            this_dm = dmr_median_delta[pid].mean(axis=1)
            ax.scatter(this_dm, this_de, edgecolor='k', facecolor='royalblue', zorder=2)
            ax.errorbar(
                this_dm,
                this_de,
                xerr=dmr_median_delta[pid].diff(axis=1).dropna(axis=1).values / 2.,
                yerr=de_logfc[pid].diff(axis=1).dropna(axis=1).values / 2.,
                ls='none',
                color='royalblue',
                zorder=1
            )
            ws = pd.DataFrame({'x': this_dm, 'y': this_de})
            xerr = dmr_median_delta[pid].diff(axis=1).dropna(axis=1).values
            yerr = de_logfc[pid].diff(axis=1).dropna(axis=1).values
            weights = pd.Series(((xerr ** 2 + yerr ** 2) ** -.5).flatten(), index=ws.index)
            wls_fit = sm.wls('y ~ x', data=ws, weights=weights).fit()
            print "Patient %s" % pid
            print wls_fit.summary()
            ax.plot(ws.x, wls_fit.predict(), color='darkgray', linewidth=1.5)
            ax.axhline(0, linestyle='--', color='k', linewidth=1.)
            ax.axvline(0, linestyle='--', color='k', linewidth=1.)
            ax.set_xlabel(r'Methylation median $\Delta M$')
            ax.set_ylabel(r'Gene expression $\log_2$ fold change')
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, "de_logfc_vs_dmr_delta_m_%s.png" % pid), dpi=200)
            fig.savefig(os.path.join(outdir, "de_logfc_vs_dmr_delta_m_%s.tiff" % pid), dpi=200)


    # test: are these DMRs significantly more concordant with DE than randomly chosen ones??
    # statistic to use: rsq of weighted linear regression
    # This will replace the statistics reported by WLS (to some extent)

    # to make this work. we need the FULL list of regions, whether DM or not
    dmr_full_indir = os.path.join(output.OUTPUT_DIR, "assess_reprog_alt1_apocrita", "results")


    n_perm = 1000
    rsq_true = {}
    rsq_permute = {}

    for pid in de_logfc:
        this_de = de_logfc[pid].mean(axis=1)
        this_dm = dmr_median_delta[pid].mean(axis=1)
        ws = pd.DataFrame({'x': this_dm, 'y': this_de})
        xerr = dmr_median_delta[pid].diff(axis=1).dropna(axis=1).values
        yerr = de_logfc[pid].diff(axis=1).dropna(axis=1).values
        weights = pd.Series(((xerr ** 2 + yerr ** 2) ** -.5).flatten(), index=ws.index)
        wls_fit = sm.wls('y ~ x', data=ws, weights=weights).fit()
        rsq_true[pid] = rsq_true

        # now repeat for randomly-selected genes (no need for significant DM, just use whatever the delta M value is)



