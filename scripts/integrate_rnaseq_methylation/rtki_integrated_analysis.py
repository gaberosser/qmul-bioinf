from methylation import dmr, process, plots
from load_data import methylation_array
from settings import GIT_LFS_DATA_DIR
from utils.output import unique_output_dir
import collections
import os
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing as mp
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import venn
from utils import setops
import json


patient_pairs = collections.OrderedDict([
    ('018', (
        ('GBM018_P10', 'GBM018_P12'), ('DURA018_NSC_N4_P4', 'DURA018_NSC_N2_P6')
    )),
    ('019', ('GBM019_P4', 'DURA019_NSC_N8C_P2')),
    ('030', ('GBM030_P5', 'DURA030_NSC_N16B6_P1')),
    ('031', ('GBM031_P4', 'DURA031_NSC_N44B_P2')),
    ('all', (
        ('GBM018_P10', 'GBM018_P12', 'GBM019_P4', 'GBM030_P5', 'GBM031_P4'),
        ('DURA018_NSC_N4_P4', 'DURA018_NSC_N2_P6', 'DURA019_NSC_N8C_P2', 'DURA030_NSC_N16B6_P1', 'DURA031_NSC_N44B_P2')
    )
     ),
])
ref_name = 'GIBCONSC_P4'


class TestResultEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, set):
            return tuple(o)
        elif isinstance(o, bool) or isinstance(o, np.bool_):
            return int(o)
        return super(TestResultEncoder, self).default(o)



def construct_contingency(x, y):
    return np.array([
        [((x < 0) & (y < 0)).sum(), ((x > 0) & (y < 0)).sum()],
        [((x < 0) & (y > 0)).sum(), ((x > 0) & (y > 0)).sum()],
    ])


def count_dmr_genes(res):
    the_genes = reduce(
        lambda x, y: x.union(y),
        [t[1]['genes'] for t in dmr.dict_iterator(res, n_level=3)],
        set()
    )
    return len(the_genes)


def n_overlap_datum(sid, joint_dmr_de, this_de, this_dmr):
    n_de = this_de[sid].shape[0]
    n_dmr = len(list(dmr.dict_iterator(this_dmr[sid], n_level=3)))
    n_dmr_genes = count_dmr_genes(this_dmr[sid])
    n_overlaps = joint_dmr_de[sid]['all'].shape[0]
    n_overlaps_unique = joint_dmr_de[sid]['all'].me_genes.unique().shape[0]
    this_datum = [n_de, n_dmr, n_dmr_genes, n_overlaps, n_overlaps_unique]
    for cls in dmr.CLASSES:
        this_datum.append(joint_dmr_de[sid][cls].shape[0])
        this_datum.append(joint_dmr_de[sid][cls].me_genes.unique().shape[0])
    return this_datum


def compute_joint_de_dmr(dmr_results, de_results):
    res = {}

    for sid in dmr_results:
        print sid
        res[sid] = {}

        de_cols = ['genes', 'logFC', 'ensembl', 'direction', 'FDR', 'logCPM']
        meth_cols = ['me_genes', 'chr', 'me_cid', 'me_mediandelta', 'me_median1', 'me_median2', 'me_fdr']
        meth_attrs = ['median_change', 'median1', 'median2']

        for (chr, cls, cid), attrs in dmr.dict_iterator(dmr_results[sid], n_level=3):
            res[sid].setdefault(cls, pd.DataFrame(columns=de_cols + meth_cols))

            if len(attrs['genes']) == 0:
                continue

            try:
                # matching entry in DE (by gene name)
                de_match = de_results[sid].loc[de_results[sid].loc[:, 'genes'].isin(attrs['genes'])]

                if de_match.shape[0] > 0:
                    # form the DMR data block by repeating the same row
                    me_data = np.tile(
                        [chr, cid] + [attrs[k] for k in meth_attrs] + [attrs['padj']],
                        (de_match.shape[0], 1)
                    )
                    me_data = np.concatenate(
                        (
                            np.reshape(de_match.genes.values, (de_match.shape[0], 1)),
                            me_data
                        ),
                        axis=1
                    )
                    me_match = pd.DataFrame(data=me_data, columns=meth_cols, index=de_match.index)

                    this_match = pd.concat((de_match, me_match), axis=1)
                    res[sid][cls] = pd.concat(
                        (res[sid][cls], this_match), axis=0, ignore_index=True
                    )
            except Exception as exc:
                print "Failed to add data: (%s, %s, %d)" % (chr, cls, cid)
                print repr(exc)
                continue

        # combine all methylation cluster classes
        res[sid]['all'] = pd.concat(res[sid].values(), axis=0, ignore_index=True)

    return res


def load_methylation_data(ref='gibco', norm_method='swan'):

    anno = methylation_array.load_illumina_methylationepic_annotation()
    b, me_meta = methylation_array.gbm_rtk1_and_paired_nsc(norm_method='swan', ref='gibco')
    b.dropna(inplace=True)
    m = process.m_from_beta(b)

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(b.index)
    anno = anno.loc[common_probes]
    b = b.loc[common_probes]
    m = m.loc[common_probes]

    # split gene symbols and store as a set
    # this should now be carried out at the point of loading
    # anno.loc[:, 'UCSC_RefGene_Name'] = \
    #     anno.UCSC_RefGene_Name.str.split(';').apply(lambda x: set(x) if isinstance(x, list) else None)

    return b, m, anno


def load_de_data(indir_de=None):
    ncol_per_de_block = 6
    if indir_de is None:
        indir_de = os.path.join(GIT_LFS_DATA_DIR, 'rnaseq_de', 'rtk1', 'insc_gibco')

    de = {}

    # blank string corresponds to full set
    for lbl in patient_pairs.keys():
        de.setdefault(lbl, {})
        # label used in DE structure
        if lbl == 'all':
            # blank string for file naming purposes
            p = ''
        else:
            p = lbl

        fn = os.path.join(indir_de, 'GBM{0}.vs.iNSC{0}-GBM{0}.vs.refNSC.csv'.format(p))

        # Data are in two 'blocks', each comprising ncol_per_block columns
        # The first block corresponds to GBM vs iNSC
        # The second block corresponds to GBM vs reference
        # Missing values mean that the gene in question features in only one comparison
        this_de = pd.read_csv(fn, header=0, index_col=None)
        in_insc = ~this_de.iloc[:, 1].isnull()
        in_ensc = ~this_de.iloc[:, ncol_per_de_block + 1].isnull()

        de[lbl]['gbm_insc'] = this_de.loc[in_insc].iloc[:, :ncol_per_de_block]
        de[lbl]['gbm_insc_only'] = this_de.loc[in_insc & ~in_ensc].iloc[:, :ncol_per_de_block]
        de[lbl]['gbm_insc_and_ref'] = this_de.loc[in_insc & in_ensc].iloc[:, :ncol_per_de_block]

    return de


def compute_dmr(mvals, anno, n_jobs=None, outdir=None, **params):
    n_jobs = n_jobs or mp.cpu_count()

    dmr.add_merged_probe_classes(anno)
    clusters = dmr.identify_clusters(anno, n_min=params['n_min'], d_max=params['d_max'], n_jobs=n_jobs)

    test_results = {}
    test_results_relevant = {}
    test_results_significant = {}

    for sid, samples in patient_pairs.items():

        test_results.setdefault(sid, {})
        test_results_relevant.setdefault(sid, {})
        test_results_significant.setdefault(sid, {})

        test_results[sid]['gbm_insc'] = dmr.test_clusters(
            clusters,
            mvals,
            samples=samples,
            min_median_change=params['delta_m_min'],
            n_jobs=n_jobs,
            method=params['dmr_test_method'],
            test_kwargs=params['test_kwargs']
        )

        # for insc and ref: use a replacement `samples`
        samples_ref = (samples[0], (ref_name,))
        test_results[sid]['gbm_ref'] = dmr.test_clusters(
            clusters,
            mvals,
            samples=samples_ref,
            min_median_change=params['delta_m_min'],
            n_jobs=n_jobs,
            method=params['dmr_test_method'],
            test_kwargs=params['test_kwargs']
        )

        for typ in comparisons:
            ## FIXME: the function  mht_correction now modifies in place; this will error
            test_results_relevant[sid][typ] = dmr.mht_correction(
                test_results[sid][typ],
                alpha=params['fdr']
            )
            test_results_significant[sid][typ] = dmr.filter_dictionary(
                test_results_relevant[sid][typ],
                filt=lambda x: x['rej_h0'],
                n_level=3
            )

    # add list of annotated genes to all clusters
    for sid in test_results:
        for (typ, chr, cls, cid), attrs in dmr.dict_iterator(test_results[sid], n_level=4):
            pids = attrs['probes']
            genes = anno.loc[attrs['probes']].UCSC_RefGene_Name.dropna()
            geneset = reduce(lambda x, y: x.union(y), genes, set())
            attrs['genes'] = geneset

    # write full set of results to disk (json format)
    # NB this converts gene sets to lists and cluster ID to a string
    if outdir is not None:
        fout = os.path.join(outdir, "dmr_results.json")
        with open(fout, 'wb') as f:
            json.dump(test_results, f, cls=TestResultEncoder)
        print "Saved DMR results to %s" % fout

    return {
        'clusters': clusters,
        'test_results': test_results,
        'test_results_relevant': test_results_relevant,
        'test_results_significant': test_results_significant,
    }


def venn_diagram_and_core_genes(meth_de, text_file, fig_file, min_overlap=4, fig_title=None, plot_genes=None):
    """

    :param meth_de:
    :param text_file:
    :param fig_file:
    :param min_overlap:
    :param fig_title:
    :param plot_genes: None, 'scatter', 'colours'
    :return:
    """
    n_sample = len(meth_de)
    # one of the samples is 'all', which we ignore
    bns = list(setops.binary_combinations_sum_gte(n_sample - 1, min_overlap))

    fig_kws = {'figsize': (8, 3.2)}
    # if fig_title is not None:
    #     fig_kws['num'] = fig_title
    fig, axs = plt.subplots(ncols=3, num=fig_title, **fig_kws)

    f = open(text_file, 'wb')
    venn_counts = {}
    venn_sets = {}
    core_sets = {}

    if plot_genes:
        fig_gene_scatter, axs_gene_scatter = plt.subplots(nrows=1, ncols=3, **fig_kws)
        fig_gene_colours, axs_gene_colours = plt.subplots(nrows=1, ncols=3, figsize=(5, 8))
        max_y = 0
        max_de = 0.
        max_dmr = 0.

    # one Venn diagram per class
    # Venn components are patients
    # we therefore build a structure containing patient-wise gene sets for each class
    for i, cls in enumerate(dmr.CLASSES):
        ax = axs[i]
        # reordered dict by sublevel
        tmp = dmr.dict_by_sublevel(meth_de, 2, cls)
        # build labels and gene lists to control the order (otherwise we're relying on dict iteration order)
        # we skip over 'all' here
        set_labels = []
        gene_sets = []
        for k, t in tmp.items():
            if k != 'all':
                gene_sets.append(set(t.genes))
                set_labels.append(k)
        venn_res, venn_sets[cls], venn_counts[cls] = venn.venn_diagram(
            *gene_sets,
            set_labels=set_labels,
            ax=ax
        )
        ax.set_title(cls)

        # Core genes are defined as those supported by min_overlap samples
        core_genes = set()
        for bn in bns:
            core_genes = core_genes.union(venn_sets[cls][bn])
        core_sets[cls] = core_genes
        core_genes = sorted(core_genes)
        print "%s core genes [%d]: %s" % (cls, len(core_genes), ', '.join(core_genes))
        f.write("%s core genes [%d]: %s\n" % (cls, len(core_genes), ', '.join(core_genes)))

        if plot_genes:
            j = 0
            de_vals = []
            dmr_vals = []
            for cg in core_genes:
                val_de = 0.
                val_dmr = 0.
                n = 0.
                # loop over patients
                for t in tmp.values():
                    if (t.genes == cg).any():
                        n += (t.genes == cg).sum()
                        the_rows = t.loc[t.genes == cg]
                        val_de += the_rows.logFC.astype(float).sum()
                        val_dmr += the_rows.me_mediandelta.astype(float).sum()
                x = val_de / n
                y = val_dmr / n
                de_vals.append(x)
                dmr_vals.append(y)
                axs_gene_scatter[i].text(val_de / n, val_dmr / n, cg)

                if x > 0 and y > 0:
                    cmap = plt.cm.Blues
                elif x <0 and y < 0:
                    cmap = plt.cm.Purples
                elif x > 0 and y < 0:
                    cmap = plt.cm.Reds
                elif x < 0 and y > 0:
                    cmap = plt.cm.Greens
                else:
                    cmap = plt.cm.Greys

                # scaling
                ## FIXME: this would be better applied afterwards, like the axis scaling
                xn = x / 12.
                yn = y / 6.
                rn = (xn ** 2 + yn ** 2) ** .5
                rn = min(max(rn, 0.3), 0.9)
                c = cmap(rn)
                axs_gene_colours[i].text(0, j, cg, color=c)
                j += 1.
                axs_gene_colours[i].axis('off')
                axs_gene_colours[i].set_title(cls)

            axs_gene_scatter[i].scatter(de_vals, dmr_vals)
            axs_gene_scatter[i].axhline(0.)
            axs_gene_scatter[i].axvline(0.)
            axs_gene_scatter[i].set_xlabel('RNASeq DE logFC')
            axs_gene_scatter[i].set_ylabel('EPIC DMR median delta M')

            max_de = max(max_de, np.abs(de_vals).max())
            max_dmr = max(max_dmr, np.abs(dmr_vals).max())
            max_y = max(max_y, j)

    if plot_genes:
        for i in range(len(dmr.CLASSES)):
            axs_gene_colours[i].set_ylim([0, max_y])
            axs_gene_scatter[i].set_xlim([-max_de * 1.2, max_de * 1.2])
            axs_gene_scatter[i].set_ylim([-max_dmr * 1.2, max_dmr * 1.2])


    # intersection of all
    core_all = sorted(reduce(set.intersection, core_sets.values()))
    print "Core genes shared across all classes [%d]: %s" % (len(core_all), ', '.join(core_all))
    f.write("Core genes shared across all classes [%d]: %s\n" % (len(core_all), ', '.join(core_all)))

    # union of all
    core_union = sorted(reduce(set.union, core_sets.values()))
    print "Core genes in >=1 class [%d]: %s" % (len(core_union), ', '.join(core_union))
    f.write("Core genes in >= 1 class [%d]: %s\n" % (len(core_union), ', '.join(core_union)))
    f.close()

    fig.tight_layout()
    fig.savefig("%s.png" % fig_file, dpi=200)
    fig.savefig("%s.pdf" % fig_file)

    if plot_genes:
        fig_gene_scatter.tight_layout()
        fig_gene_scatter.savefig("%s_genescatter.png" % fig_file, dpi=200)
        fig_gene_scatter.savefig("%s_genescatter.pdf" % fig_file)

        fig_gene_colours.tight_layout()
        fig_gene_colours.savefig("%s_genecolours.png" % fig_file, dpi=200)
        fig_gene_colours.savefig("%s_genecolours.pdf" % fig_file)

    return venn_sets


def cohort_core_genes(meth_de, text_file):
    core_sets = {}
    with open(text_file, 'wb') as f:
        for i, cls in enumerate(dmr.CLASSES):
            core_sets[cls] = sorted(set(meth_de['all'][cls].genes))
            print "%s core genes [%d]: %s" % (cls, len(core_sets[cls]), ', '.join(core_sets[cls]))
            f.write("%s core genes [%d]: %s\n" % (cls, len(core_sets[cls]), ', '.join(core_sets[cls])))

        # intersection of all
        core_all = sorted(reduce(lambda x, y: set(x).intersection(y), core_sets.values()))
        print "Core genes shared across all classes [%d]: %s" % (len(core_all), ', '.join(core_all))
        f.write("Core genes shared across all classes [%d]: %s\n" % (len(core_all), ', '.join(core_all)))

        # union of all
        core_union = sorted(reduce(lambda x, y: set(x).union(y), core_sets.values()))
        print "Core genes in >=1 class [%d]: %s" % (len(core_union), ', '.join(core_union))
        f.write("Core genes in >= 1 class [%d]: %s\n" % (len(core_union), ', '.join(core_union)))


if __name__ == '__main__':

    params = {
        'core_min_sample_overlap': 3,  # 3 / 4 samples must match
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'fdr': 0.05,
        # 'dmr_test_method': 'mwu_permute',  # 'mwu', 'mwu_permute'
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        # 'test_kwargs': {'n_max': 1999},
        'test_kwargs': {},
    }

    # set outdir to use pre-computed results
    outdir = unique_output_dir(
        "rtk1_de_dmr.%s.gibco_reference" % params['dmr_test_method'],
        reuse_empty=True
    )
    indir_de = os.path.join(GIT_LFS_DATA_DIR, 'rnaseq_de', 'rtk1', 'insc_gibco')  # can also use H9

    # write params to a file
    with open(os.path.join(outdir, 'parameters.json'), 'wb') as f:
        json.dump(params, f)

    ref_name = 'GIBCONSC_P4'
    comparisons = ['gbm_insc', 'gbm_ref']
    n_jobs = mp.cpu_count()

    de = load_de_data(indir_de)

    b, m, anno = load_methylation_data(ref='gibco')

    # compute or load DMR
    res = compute_dmr(m, anno, outdir=outdir, **params)
    clusters = res['clusters']
    test_results = res['test_results']
    test_results_relevant = res['test_results_relevant']
    test_results_significant = res['test_results_significant']

    # create a table of the numbers of DMRs
    cols = (
        'sample', 'clusters_proposed', 'clusters_relevant', 'clusters_significant',
        'genes_proposed', 'genes_relevant', 'genes_significant'
    )
    table_cluster_numbers = {}

    # the number of proposed clusters and proposed genes is constant across all results, so just use the first
    k0 = patient_pairs.keys()[0]
    k1 = test_results.values()[0].keys()[0]
    ncl = len(list(
        dmr.dict_iterator(test_results[k0][k1], n_level=3)
    ))
    ng = count_dmr_genes(test_results[k0][k1])

    for sid in test_results:
        for typ in comparisons:
            table_cluster_numbers.setdefault(typ, pd.DataFrame(columns=cols))
            ncl_re = len(list(
                dmr.dict_iterator(test_results_relevant[sid][typ], n_level=3)
            ))
            ncl_si = len(list(
                dmr.dict_iterator(test_results_significant[sid][typ], n_level=3)
            ))
            ng_re = count_dmr_genes(test_results_relevant[sid][typ])
            ng_si = count_dmr_genes(test_results_significant[sid][typ])
            this_row = pd.Series({
                'sample': sid,
                'clusters_proposed': ncl,
                'clusters_relevant': ncl_re,
                'clusters_significant': ncl_si,
                'genes_proposed': ng,
                'genes_relevant': ng_re,
                'genes_significant': ng_si
            })
            table_cluster_numbers[typ] = table_cluster_numbers[typ].append(this_row, ignore_index=True)

    for typ in comparisons:
        table_cluster_numbers[typ].to_csv(os.path.join(outdir, "cluster_numbers.%s.csv" % typ))

    # split into (insc not ref) and (insc and ref) groups, like DE

    this_dmr_insc = dmr.dict_by_sublevel(test_results_significant, 2, 'gbm_insc')
    this_dmr_ref = dmr.dict_by_sublevel(test_results_significant, 2, 'gbm_ref')

    this_dmr_insc_only = {}
    this_dmr_insc_ref = {}

    for sid in test_results_significant:

        this_dmr_insc_only.setdefault(sid, {})
        this_dmr_insc_ref.setdefault(sid, {})

        for chr in this_dmr_insc[sid]:
            this_dmr_insc_only[sid].setdefault(chr, {})
            this_dmr_insc_ref[sid].setdefault(chr, {})
            for cls in dmr.CLASSES:
                this_dmr_insc_only[sid][chr].setdefault(cls, {})
                this_dmr_insc_ref[sid][chr].setdefault(cls, {})

                in_insc = set(this_dmr_insc[sid][chr][cls].keys())
                in_ref = set(this_dmr_ref[sid][chr][cls].keys())
                in_insc_only = in_insc.difference(in_ref)
                in_insc_and_ref = in_insc.intersection(in_ref)

                # use the iNSC results in both cases, as that is what we are interested in
                this_dmr_insc_only[sid][chr][cls] = dict([
                    (cid, this_dmr_insc[sid][chr][cls][cid]) for cid in in_insc_only
                ])
                this_dmr_insc_ref[sid][chr][cls] = dict([
                    (cid, this_dmr_insc[sid][chr][cls][cid]) for cid in in_insc_and_ref
                ])

    # counts of clusters and genes in DMR
    cols += (
        'clusters_significant_gbm_insc',
        'genes_significant_gbm_insc',
        'clusters_significant_gbm_ref',
        'genes_significant_gbm_ref',
    )
    table_clusters_ref_insc = pd.DataFrame(columns=cols)

    for sid in this_dmr_insc:
        old_row = table_cluster_numbers['gbm_insc'].loc[
            table_cluster_numbers['gbm_insc'].loc[:, 'sample'] == sid
        ]
        this_row = pd.Series(old_row.iloc[0], copy=True)
        this_row.loc['clusters_significant_gbm_insc'] = len(list(
            dmr.dict_iterator(this_dmr_insc_only[sid], n_level=3)
        ))
        this_row.loc['clusters_significant_gbm_ref'] = len(list(
            dmr.dict_iterator(this_dmr_insc_ref[sid], n_level=3)
        ))
        this_row.loc['genes_significant_gbm_insc'] = count_dmr_genes(this_dmr_insc_only[sid])
        this_row.loc['genes_significant_gbm_ref'] = count_dmr_genes(this_dmr_insc_ref[sid])
    table_clusters_ref_insc.to_csv(os.path.join(outdir, "cluster_numbers_gbm_ref.csv"))

    """
    Before we start combining DMR and DE, I want to check the overlap between DMR on individuals and
    DMR on the pooled samples.
    We use the chr, class and cluster ID as unique IDs here.
    """
    print "(iNSC vs GBM) AND NOT (iNSC vs reference)"
    for cls in dmr.CLASSES:
        a = dmr.dict_by_sublevel(
            this_dmr_insc_only['all'],
            2,
            cls
        )
        ga = dmr.dict_iterator(a, n_level=2)
        a_keys = set([tuple(t[0]) for t in ga])

        b_keys = set()
        for sid in [t for t in patient_pairs if t != 'all']:
            gb = dmr.dict_iterator(
                dmr.dict_by_sublevel(this_dmr_insc_only[sid], 2, cls),
                n_level=2
            )
            b_keys.update(set([tuple(t[0]) for t in gb]))

        print "***"
        print "Individual comparisons: %d DMRs in class %s" % (len(a_keys), cls)
        print "Lumped comparisons: %d DMRs in class %s" % (len(b_keys), cls)
        print "Intersection of those: %d DMRs in class %s" % (len(b_keys.intersection(a_keys)), cls)
    print "***"

    print "(iNSC vs GBM) AND (iNSC vs reference)"
    for cls in dmr.CLASSES:
        a = dmr.dict_by_sublevel(
            this_dmr_insc_ref['all'],
            2,
            cls
        )
        ga = dmr.dict_iterator(a, n_level=2)
        a_keys = set([tuple(t[0]) for t in ga])

        b_keys = set()
        for sid in [t for t in patient_pairs if t != 'all']:
            gb = dmr.dict_iterator(
                dmr.dict_by_sublevel(this_dmr_insc_ref[sid], 2, cls),
                n_level=2
            )
            b_keys.update(set([tuple(t[0]) for t in gb]))

        print "***"
        print "Individual comparisons: %d DMRs in class %s" % (len(a_keys), cls)
        print "Lumped comparisons: %d DMRs in class %s" % (len(b_keys), cls)
        print "Intersection of those: %d DMRs in class %s" % (len(b_keys.intersection(a_keys)), cls)
    print "***"

    # 1: What is the joint distribution of methylation / mRNA fold change?
    # Get methylation level and DE fold change for linked genes (pairwise only)

    this_de_insc = dmr.dict_by_sublevel(de, 2, 'gbm_insc')
    this_de_insc_only = dmr.dict_by_sublevel(de, 2, 'gbm_insc_only')
    this_de_insc_ref = dmr.dict_by_sublevel(de, 2, 'gbm_insc_and_ref')

    joint_de_dmr = {
        ('insc', 'insc'): compute_joint_de_dmr(this_dmr_insc, this_de_insc),
        ('insc_only', 'insc_only'): compute_joint_de_dmr(this_dmr_insc_only, this_de_insc_only),
        ('insc_ref', 'insc_ref'): compute_joint_de_dmr(this_dmr_insc_ref, this_de_insc_ref),
        ('insc_only', 'insc'): compute_joint_de_dmr(this_dmr_insc_only, this_de_insc),
        ('insc_ref', 'insc'): compute_joint_de_dmr(this_dmr_insc_ref, this_de_insc),
        ('insc', 'insc_only'): compute_joint_de_dmr(this_dmr_insc, this_de_insc_only),
        ('insc', 'insc_ref'): compute_joint_de_dmr(this_dmr_insc, this_de_insc_ref),
    }


    joint_de_insc_dmr_insc = compute_joint_de_dmr(this_dmr_insc, this_de_insc)

    joint_de_insc_only_dmr_insc_only = compute_joint_de_dmr(this_dmr_insc_only, this_de_insc_only)
    joint_de_insc_ref_dmr_insc_ref = compute_joint_de_dmr(this_dmr_insc_ref, this_de_insc_ref)

    # compute the overlap with split DE and unsplit DMR
    joint_de_insc_only_dmr_insc = compute_joint_de_dmr(this_dmr_insc, this_de_insc_only)
    joint_de_insc_ref_dmr_insc = compute_joint_de_dmr(this_dmr_insc, this_de_insc_ref)

    # compute the overlap with unsplit DE and split DMR
    joint_de_insc_dmr_insc_only = compute_joint_de_dmr(this_dmr_insc_only, this_de_insc)
    joint_de_insc_dmr_insc_ref = compute_joint_de_dmr(this_dmr_insc_ref, this_de_insc)

    # Generate table giving the number of overlaps in each patient and cluster class
    # this includes the number of absolute overlaps AND the number of unique overlaps
    the_cols = ['DE genes', 'DMR', 'DMR genes', 'overlaps', 'unique overlaps']
    the_cols += reduce(lambda x, y: x + y, [['%s' % t, '%s_unique' % t] for t in dmr.CLASSES], [])
    count_overlap_de_insc_only_dmr_insc_only = pd.DataFrame(
        columns=the_cols,
        index=pd.Index(patient_pairs, name='patient'),
    )
    count_overlap_de_insc_ref_dmr_insc_ref = pd.DataFrame.copy(count_overlap_de_insc_only_dmr_insc_only)
    count_overlap_de_insc_dmr_insc_only = pd.DataFrame.copy(count_overlap_de_insc_only_dmr_insc_only)
    count_overlap_de_insc_dmr_insc_ref = pd.DataFrame.copy(count_overlap_de_insc_only_dmr_insc_only)
    count_overlap_de_insc_only_dmr_insc = pd.DataFrame.copy(count_overlap_de_insc_only_dmr_insc_only)
    count_overlap_de_insc_ref_dmr_insc = pd.DataFrame.copy(count_overlap_de_insc_only_dmr_insc_only)
    count_overlap_de_insc_dmr_insc = pd.DataFrame.copy(count_overlap_de_insc_only_dmr_insc_only)

    for sid in test_results:
        count_overlap_de_insc_only_dmr_insc_only.loc[sid] = n_overlap_datum(sid, joint_de_insc_only_dmr_insc_only, this_de_insc_only, this_dmr_insc_only)
        count_overlap_de_insc_ref_dmr_insc_ref.loc[sid] = n_overlap_datum(sid, joint_de_insc_ref_dmr_insc_ref, this_de_insc_ref, this_dmr_insc_ref)
        count_overlap_de_insc_dmr_insc_only.loc[sid] = n_overlap_datum(sid, joint_de_insc_dmr_insc_only, this_de_insc, this_dmr_insc_only)
        count_overlap_de_insc_dmr_insc_ref.loc[sid] = n_overlap_datum(sid, joint_de_insc_dmr_insc_ref, this_de_insc, this_dmr_insc_ref)
        count_overlap_de_insc_only_dmr_insc.loc[sid] = n_overlap_datum(sid, joint_de_insc_only_dmr_insc, this_de_insc_only, this_dmr_insc)
        count_overlap_de_insc_ref_dmr_insc.loc[sid] = n_overlap_datum(sid, joint_de_insc_ref_dmr_insc, this_de_insc_ref, this_dmr_insc)
        count_overlap_de_insc_dmr_insc.loc[sid] = n_overlap_datum(sid, joint_de_insc_dmr_insc, this_de_insc, this_dmr_insc)

    def scatter_plot_dmr_de(meth_de, fig_filestem, fig_titlestem=''):
        for sid in meth_de:
            fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, num="%s %s" % (fig_titlestem, sid))
            for i, cls in enumerate(['all', 'tss', 'gene', 'island']):
                ax = axs.flat[i]

                # get values for ALL DMR clusters of this class
                x = meth_de[sid][cls].loc[:, 'logFC'].astype(float)
                y = meth_de[sid][cls].loc[:, 'me_mediandelta'].astype(float)

                # contingency table for Fisher's exact test
                conting = construct_contingency(x, y)
                logodds, fisherp = stats.fisher_exact(conting)

                # get values for DMR clusters that are ONLY in this class (no overlaps)
                if cls == 'all':
                    print "%s - %s p = %.3e" % (
                        sid, cls, fisherp
                    )
                    ax.scatter(x, y, c='k')
                    ax.axhline(0, c=0.4 * np.ones(3))
                    ax.axvline(0, c=0.4 * np.ones(3))
                    ttl = "%s; p={0}" % cls
                    if fisherp < 0.001:
                        ttl = ttl.format('%.2e') % fisherp
                    else:
                        ttl = ttl.format('%.3f') % fisherp

                else:
                    cid_other = set(
                        np.unique(
                            np.concatenate([meth_de[sid][t].me_cid.values for t in dmr.CLASSES.difference({cls, })])
                        )
                    )
                    xu = x.loc[~meth_de[sid][cls].loc[:, 'me_cid'].isin(cid_other)].astype(float)
                    yu = y.loc[~meth_de[sid][cls].loc[:, 'me_cid'].isin(cid_other)].astype(float)
                    contingu = construct_contingency(xu, yu)
                    logoddsu, fisherpu = stats.fisher_exact(contingu)

                    print "%s - %s p = %.3e (incl overlaps), p = %.3e (excl overlaps)" % (
                        sid, cls, fisherp, fisherpu
                    )

                    ax.scatter(x, y, c='k')
                    ax.scatter(xu, yu, c='b')
                    ax.axhline(0, c=0.4 * np.ones(3))
                    ax.axvline(0, c=0.4 * np.ones(3))

                    ttl = "%s; p={0}; unique p={0}" % cls

                    if fisherp < 0.001:
                        ttl = ttl.format('%.2e') % (fisherp, fisherpu)
                    else:
                        ttl = ttl.format('%.3f') % (fisherp, fisherpu)

                ax.set_title(ttl)

            fig.text(0.5, 0.04, 'RNASeq DE logFC', ha='center')
            fig.text(0.04, 0.5, 'EPIC DMR median delta M', va='center', rotation='vertical')
            fig.savefig("%s_%s.png" % (fig_filestem, sid), dpi=200)
            fig.savefig("%s_%s.png" % (fig_filestem, sid), dpi=200)

    print "*** DE: iNSC and not ref; DMR: iNSC and not ref ***"
    scatter_plot_dmr_de(joint_de_insc_only_dmr_insc_only, os.path.join(outdir, "de_insc_only_dmr_insc_only"), "DE: iNSC not ref; DMR: iNSC not ref")
    print "*** DE: iNSC and ref; DMR: iNSC and ref ***"
    scatter_plot_dmr_de(joint_de_insc_ref_dmr_insc_ref, os.path.join(outdir, "de_insc_ref_dmr_insc_ref"), "DE: iNSC and ref; DMR: iNSC and ref")

    print "*** DE: iNSC; DMR: iNSC not ref ***"
    scatter_plot_dmr_de(joint_de_insc_dmr_insc_only, os.path.join(outdir, "de_insc_dmr_insc_only"), "DE: iNSC; DMR: iNSC not ref")
    print "*** DE: iNSC; DMR: iNSC and ref ***"
    scatter_plot_dmr_de(joint_de_insc_dmr_insc_ref, os.path.join(outdir, "de_insc_dmr_insc_ref"), "DE: iNSC; DMR: iNSC and ref")

    print "*** DE: iNSC not ref; DMR: iNSC ***"
    scatter_plot_dmr_de(joint_de_insc_only_dmr_insc, os.path.join(outdir, "de_insc_only_dmr_insc"), "DE: iNSC not ref; DMR: iNSC")
    print "*** DE: iNSC and ref; DMR: iNSC ***"
    scatter_plot_dmr_de(joint_de_insc_ref_dmr_insc, os.path.join(outdir, "de_insc_ref_dmr_insc"), "DE: iNSC and ref; DMR: iNSC")

    print "*** DE: iNSC; DMR: iNSC ***"
    scatter_plot_dmr_de(joint_de_insc_dmr_insc, os.path.join(outdir, "de_insc_dmr_insc"), "DE: iNSC; DMR: iNSC")

    # 2: Venn diagrams: to what extent do the same genes appear in all RTK 1 samples?

    # since a 5-way Venn is not supported (and is horrible to look at even if you can draw it...), remove the
    # 'all' results here and report them separately

    venn_results = {}
    bns = list(setops.binary_combinations_sum_gte(len(patient_pairs) - 1, len(patient_pairs) - 2))
    for ii in ['insc_only', 'insc_ref', 'insc']:
        for jj in ['insc_only', 'insc_ref', 'insc']:
            if (ii, jj) not in joint_de_dmr:
                continue
            print "*** DE: %s; DMR: %s ***" % (ii, jj)
            venn_results[(ii, jj)] = venn_diagram_and_core_genes(
                joint_de_dmr[(ii, jj)],
                os.path.join(outdir, "core_genes_joint_de_%s_dmr_%s.by_patient.txt" % (ii, jj)),
                os.path.join(outdir, "venn_joint_de_%s_dmr_%s" % (ii, jj)),
                min_overlap=params['core_min_sample_overlap'],
                fig_title='DE: %s; DMR: %s' % (ii, jj),
                plot_genes=True,
            )


    print "*** DE: iNSC and not ref; DMR: iNSC and not ref - individual patient ***"
    venn_results[('insc_only', 'insc_only')] = venn_diagram_and_core_genes(
        joint_de_insc_only_dmr_insc_only,
        os.path.join(outdir, "core_genes_joint_de_insc_only_dmr_insc_only.by_patient.txt"),
        os.path.join(outdir, "venn_joint_de_insc_only_dmr_insc_only"),
        min_overlap=params['core_min_sample_overlap'],
        fig_title='DE: iNSC and not ref; DMR: iNSC and not ref'
    )

    print "*** DE: iNSC and ref; DMR: iNSC and ref - individual patient ***"
    venn_results[('insc_ref', 'insc_ref')] = venn_diagram_and_core_genes(
        joint_de_insc_ref_dmr_insc_ref,
        os.path.join(outdir, "core_genes_joint_de_insc_ref_dmr_insc_ref.by_patient.txt"),
        os.path.join(outdir, "venn_joint_de_insc_ref_dmr_insc_ref"),
        min_overlap=params['core_min_sample_overlap'],
        fig_title='DE: iNSC and ref; DMR: iNSC and ref'
    )

    print "*** DE: iNSC; DMR: iNSC not ref - individual patient ***"
    venn_results[('insc', 'insc_only')] = venn_diagram_and_core_genes(
        joint_de_insc_dmr_insc_only,
        os.path.join(outdir, "core_genes_joint_de_insc_dmr_insc_only.by_patient.txt"),
        os.path.join(outdir, "venn_joint_de_insc_dmr_insc_only"),
        min_overlap=params['core_min_sample_overlap'],
        fig_title='DE: iNSC; DMR: iNSC not ref'
    )
    print "*** DE: iNSC; DMR: iNSC and ref - individual patient ***"
    venn_results[('insc', 'insc_ref')] = venn_diagram_and_core_genes(
        joint_de_insc_dmr_insc_ref,
        os.path.join(outdir, "core_genes_joint_de_insc_dmr_insc_ref.by_patient.txt"),
        os.path.join(outdir, "venn_joint_de_insc_dmr_insc_ref"),
        min_overlap=params['core_min_sample_overlap'],
        fig_title='DE: iNSC; DMR: iNSC and ref'
    )

    print "*** DE: iNSC not ref; DMR: iNSC - individual patient ***"
    venn_results[('insc_only', 'insc')] = venn_diagram_and_core_genes(
        joint_de_insc_only_dmr_insc,
        os.path.join(outdir, "core_genes_joint_de_insc_only_dmr_insc.by_patient.txt"),
        os.path.join(outdir, "venn_joint_de_insc_only_dmr_insc"),
        min_overlap=params['core_min_sample_overlap'],
        fig_title='DE: iNSC not ref; DMR: iNSC'
    )
    print "*** DE: iNSC and ref; DMR: iNSC - individual patient ***"
    venn_results[('insc_ref', 'insc')] = venn_diagram_and_core_genes(
        joint_de_insc_ref_dmr_insc,
        os.path.join(outdir, "core_genes_joint_de_insc_ref_dmr_insc.by_patient.txt"),
        os.path.join(outdir, "venn_joint_de_insc_ref_dmr_insc"),
        min_overlap=params['core_min_sample_overlap'],
        fig_title='DE: iNSC and ref; DMR: iNSC'
    )

    print "*** DE: iNSC; DMR: iNSC - individual patient ***"
    venn_results[('insc', 'insc')] = venn_diagram_and_core_genes(
        joint_de_insc_dmr_insc,
        os.path.join(outdir, "core_genes_joint_de_insc_dmr_insc.by_patient.txt"),
        os.path.join(outdir, "venn_joint_de_insc_dmr_insc"),
        min_overlap=params['core_min_sample_overlap'],
        fig_title='DE: iNSC; DMR: iNSC'
    )

    # repeat with whole cohort (core genes only, no Venn)

    print "*** DE: iNSC and not ref; DMR: iNSC and not ref - whole cohort ***"
    cohort_core_genes(
        joint_de_insc_only_dmr_insc_only,
        os.path.join(outdir, "core_genes_joint_de_insc_only_dmr_insc_only.by_cohort.txt"),
    )
    print "*** DE: iNSC and ref; DMR: iNSC and ref - whole cohort ***"
    cohort_core_genes(
        joint_de_insc_ref_dmr_insc_ref,
        os.path.join(outdir, "core_genes_joint_de_insc_ref_dmr_insc_ref.by_cohort.txt"),
    )

    print "*** DE: iNSC; DMR: iNSC not ref - whole cohort ***"
    cohort_core_genes(
        joint_de_insc_dmr_insc_only,
        os.path.join(outdir, "core_genes_joint_de_insc_dmr_insc_only.by_cohort.txt"),
    )
    print "*** DE: iNSC; DMR: iNSC and ref - whole cohort ***"
    cohort_core_genes(
        joint_de_insc_dmr_insc_ref,
        os.path.join(outdir, "core_genes_joint_de_insc_dmr_insc_ref.by_cohort.txt"),
    )

    print "*** DE: iNSC not ref; DMR: iNSC - whole cohort ***"
    cohort_core_genes(
        joint_de_insc_only_dmr_insc,
        os.path.join(outdir, "core_genes_joint_de_insc_only_dmr_insc.by_cohort.txt"),
    )
    print "*** DE: iNSC and ref; DMR: iNSC - whole cohort ***"
    cohort_core_genes(
        joint_de_insc_ref_dmr_insc,
        os.path.join(outdir, "core_genes_joint_de_insc_ref_dmr_insc.by_cohort.txt"),
    )

    print "*** DE: iNSC; DMR: iNSC - whole cohort ***"
    cohort_core_genes(
        joint_de_insc_dmr_insc,
        os.path.join(outdir, "core_genes_joint_de_insc_dmr_insc.by_cohort.txt"),
    )

    # represent these results on a plot




    # check the direction of change in each sample - these must be consistent for the gene to be considered core
    de_values = {}
    dmr_values = {}
    for i, cls in enumerate(dmr.CLASSES):
        # meth_de = dmr.dict_by_sublevel(meth_de_joint_insc_ref, 2, cls)
        meth_de = dmr.dict_by_sublevel(joint_de_insc_only_dmr_insc_only, 2, cls)

        per_sample_de = dict([(k, set(x.genes)) for k, x in meth_de.items()])
        all_genes = reduce(set.union, per_sample_de.values())
        de_values[cls] = pd.DataFrame(index=sorted(all_genes), columns=patient_pairs.keys())

        # Step 1: check each gene for DE compatibility
        for g in all_genes:
            this_row = []
            for sid, this_dat in meth_de.items():
                this_genes = this_dat.genes
                if g in this_genes.values:
                    de_values[cls].loc[g, sid] = this_dat.loc[this_genes == g, 'logFC'].values[0]

        # find those genes shared by >1 sample
        role_call = (~de_values[cls].isnull()).sum(axis=1)
        shared_de = de_values[cls].loc[
            role_call > 1
        ].fillna(0)
        obs = np.sign(shared_de).astype(int).sum(axis=1).abs()
        expctd = shared_de.astype(bool).astype(int).sum(axis=1)
        if not (obs == expctd).all():
            to_remove = obs.loc[obs != expctd].index
            print "Not all DE genes had the same direction change. Removing %s." % ', '.join(to_remove)
            ## TODO: remove
            raise NotImplementedError

        # Step 2: check each matching DMR for compatibility
        per_sample_dmr = dict(
            [(k, set([tuple(t) for t in x.loc[:, ['chr', 'me_cid']].values])) for k, x in meth_de.items()]
        )
        all_dmr = reduce(set.union, per_sample_dmr.values())

        dmr_values[cls] = pd.DataFrame(index=sorted(all_dmr), columns=patient_pairs.keys())
        for chr, cid in all_dmr:
            this_row = []
            for sid, this_dat in meth_de.items():
                idx = ((this_dat.loc[:, 'chr'] == chr) & (this_dat.loc[:, 'me_cid'] == cid))
                if idx.any():
                    dmr_values[cls].loc[(chr, cid), sid] = this_dat.loc[idx, 'me_mediandelta'].values[0]
        role_call = (~dmr_values[cls].isnull()).sum(axis=1)
        shared_dmr = dmr_values[cls].loc[
            role_call > 1
        ].fillna(0)
        obs = np.sign(shared_dmr).astype(int).sum(axis=1).abs()
        expctd = shared_dmr.astype(bool).astype(int).sum(axis=1)
        if not (obs == expctd).all():
            to_remove = obs.loc[obs != expctd].index
            print "Not all DE genes had the same direction change. Removing %s." % ', '.join(to_remove)
            ## TODO: remove
            raise NotImplementedError


    # OK!




