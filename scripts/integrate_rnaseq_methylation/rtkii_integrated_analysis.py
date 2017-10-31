from load_data import rnaseq_data, methylation_array
from rnaseq import differential_expression
from rnaseq.filter import filter_by_cpm
from methylation import process, dmr
from methylation.plots import venn_dmr_counts, dmr_overlap
from integrator.rnaseq_methylationarray import compute_joint_de_dmr
from integrator import plots
import pandas as pd
import numpy as np
from scipy import stats
import references
import os
import itertools
import csv
import datetime
import json
from utils import output, setops
from matplotlib import pyplot as plt
import seaborn as sns


class TestResultEncoder(json.JSONEncoder):
    """
    Handles saving results that may contain Bool or sets
    """
    def default(self, o):
        if isinstance(o, set):
            return tuple(o)
        elif isinstance(o, bool) or isinstance(o, np.bool_):
            return int(o)
        return super(TestResultEncoder, self).default(o)


def add_gene_symbols(df):
    """
    Add gene symbols to the DataFrame df which is indexed by Ensembl IDs
    """
    gs = references.ensembl_to_gene_symbol(df.index)
    # resolve any duplicates arbitrarily (these should be rare)
    gs = gs.loc[~gs.index.duplicated()]
    df.insert(0, 'Gene Symbol', gs)


def add_fc_direction(df):
    direction = pd.Series(index=df.index, name='Direction')
    direction.loc[df.logFC < 0] = 'down'
    direction.loc[df.logFC > 0] = 'up'
    df.insert(df.shape[1], 'Direction', direction)


def compute_de(rnaseq_obj, pids, lfc=1, fdr=0.01, method='QLGLM'):
    if method not in {'QLGLM', 'GLM', 'exact'}:
        raise NotImplementedError("Unsupported method.")
    de_matched = {}
    de_gibco = {}
    de = {}

    # individual DE

    for pid in pids:
        try:
            the_idx = rnaseq_obj.meta.index.str.contains(pid)
            the_data = rnaseq_obj.data.loc[:, the_idx]

            the_data = filter_by_cpm(the_data, min_n_samples=1)
            the_genes = the_data.index

            the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
            the_contrast = "GBM - iNSC"
            the_pair = ['iNSC', 'GBM']

            if method == 'QLGLM':
                de_matched[pid] = differential_expression.edger_glmqlfit(the_data, the_groups, the_contrast, lfc=lfc, fdr=fdr)
            elif method == 'GLM':
                de_matched[pid] = differential_expression.edger_glmfit(the_data, the_groups, the_contrast, lfc=lfc, fdr=fdr)
            elif method == 'exact':
                de_matched[pid] = differential_expression.edger_exacttest(the_data, the_groups, pair=the_pair, lfc=lfc, fdr=fdr)

            add_gene_symbols(de_matched[pid])
            add_fc_direction(de_matched[pid])

            # repeat with gibco reference
            # use the same genes, rather than filtering again
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.type == 'GBM')) | rnaseq_obj.meta.index.str.contains('GIBCO')
            the_data = rnaseq_obj.data.loc[the_genes, the_idx]
            the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
            the_contrast = "GBM - NSC"
            the_pair = ['NSC', 'GBM']

            if method == 'QLGLM':
                de_gibco[pid] = differential_expression.edger_glmqlfit(the_data, the_groups, the_contrast, lfc=lfc, fdr=fdr)
            elif method == 'GLM':
                de_gibco[pid] = differential_expression.edger_glmfit(the_data, the_groups, the_contrast, lfc=lfc, fdr=fdr)
            elif method == 'exact':
                de_gibco[pid] = differential_expression.edger_exacttest(the_data, the_groups, pair=the_pair, lfc=lfc, fdr=fdr)

            add_gene_symbols(de_gibco[pid])
            add_fc_direction(de_gibco[pid])

            # Separate into sets
            de[pid], _ = setops.venn_from_arrays(de_matched[pid].index, de_gibco[pid].index)

            # TODO: finally, modify any entries that are discordant in direction (up/down) in the intersection set?

        except Exception as exc:
            print "Patient %s failed." % pid
            print repr(exc)

    # TODO: add a whole-group comparison?

    return {
        'de': de,
        'de_matched': de_matched,
        'de_gibco': de_gibco
    }


def compute_dmr(me_data, me_meta, anno, pids, dmr_params, ref_name_contains='GIBCO'):

    clusters = dmr.identify_clusters(anno, **dmr_params)

    test_results = {}
    test_results_relevant = {}
    test_results_significant = {}

    for pid in pids:

        test_results.setdefault(pid, {})
        test_results_relevant.setdefault(pid, {})
        test_results_significant.setdefault(pid, {})

        if pid == 'all':
            the_idx = ~me_meta.index.str.contains(ref_name_contains)
        else:
            the_idx = me_meta.index.str.contains(pid)
        the_groups = me_meta.loc[the_idx, 'type'].values
        the_samples = me_meta.index[the_idx].groupby(the_groups).values()

        # matched comparison
        test_results[pid]['matched'] = dmr.test_clusters(
            clusters,
            me_data,
            samples=the_samples,
            n_jobs=dmr_params['n_jobs'],
            min_median_change=dmr_params['delta_m_min'],
            method=dmr_params['dmr_test_method'],
            **dmr_params['test_kwargs']
        )

        # reference comparison
        if pid == 'all':
            the_idx = (me_meta.type != 'iNSC')
        else:
            the_idx = (me_meta.index.str.contains(pid) & (me_meta.type == 'GBM')) | me_meta.index.str.contains(ref_name_contains)
        the_groups = me_meta.loc[the_idx, 'type'].values
        the_samples = me_meta.index[the_idx].groupby(the_groups).values()
        test_results[pid]['gibco'] = dmr.test_clusters(
            clusters,
            me_data,
            samples=the_samples,
            n_jobs=dmr_params['n_jobs'],
            min_median_change=dmr_params['delta_m_min'],
            method=dmr_params['dmr_test_method'],
            **dmr_params['test_kwargs']
        )

        for typ in ['matched', 'gibco']:
            test_results_relevant[pid][typ] = dmr.mht_correction(
                test_results[pid][typ],
                alpha=dmr_params['fdr']
            )
            test_results_significant[pid][typ] = dmr.filter_dictionary(
                test_results_relevant[pid][typ],
                filt=lambda x: x['rej_h0'],
                n_level=3
            )

    # add list of annotated genes to all clusters
    for pid in test_results:
        for _, attrs in dmr.dict_iterator(test_results[pid], n_level=4):
            genes = anno.loc[attrs['probes']].UCSC_RefGene_Name.dropna()
            geneset = reduce(lambda x, y: x.union(y), genes, set())
            attrs['genes'] = geneset

    return {
        'clusters': clusters,
        'results': test_results,
        'relevant': test_results_relevant,
        'significant': test_results_significant,
    }


def tabulate_dmr_results(test_results, outdir=None, **other_results):
    """
    :param dmr_results: Dictionary giving different representations, typically 'proposed', 'relevant' and
    'significant'. Keys are used to name columns.
    :param outdir: If supplied, tables are written to this directory
    :return:
    """
    pids = test_results.keys()
    types = test_results.values()[0].keys()
    the_keys = ['proposed'] + other_results.keys()
    cols = ['sample'] + ['clusters_%s' % t for t in the_keys] + ['genes_%s' % t for t in the_keys]
    table_cluster_numbers = {}

    # the number of proposed clusters and proposed genes is constant across all results, so just use the first
    ncl = len(list(
        dmr.dict_iterator(test_results[pids[0]][types[0]], n_level=3)
    ))
    ng = dmr.count_dmr_genes(test_results[pids[0]][types[0]])
    row_template = {'clusters_proposed': ncl, 'genes_proposed': ng}

    for pid in pids:
        for typ in types:
            table_cluster_numbers.setdefault(typ, pd.DataFrame(columns=cols))
            this_row = dict(row_template)
            this_row['sample'] = pid
            for suffix in other_results.keys():
                this_kwarg = other_results[suffix]
                this_row['clusters_%s' % suffix] = len(list(
                    dmr.dict_iterator(this_kwarg[pid][typ], n_level=3)
                ))
                this_row['genes_%s' % suffix] = dmr.count_dmr_genes(this_kwarg[pid][typ])
            this_row = pd.Series(this_row)
            table_cluster_numbers[typ] = table_cluster_numbers[typ].append(this_row, ignore_index=True)

    if outdir is not None:
        for typ in types:
            table_cluster_numbers[typ].to_csv(os.path.join(outdir, "cluster_numbers.%s.csv" % typ))

    return table_cluster_numbers


def dmr_venn_sets(test_results):
    test_results_reindex = dmr.dict_by_sublevel(test_results, 2, 'matched')
    test_results_exclusive = {}
    test_results_inclusive = {}

    for (pid, chr, cls), attr in dmr.dict_iterator(test_results_reindex, n_level=3):

        test_results_exclusive.setdefault(pid, {}).setdefault('matched', {}).setdefault(chr, {})
        test_results_exclusive.setdefault(pid, {}).setdefault('gibco', {}).setdefault(chr, {})
        test_results_inclusive.setdefault(pid, {}).setdefault('matched', {}).setdefault(chr, {})
        test_results_inclusive.setdefault(pid, {}).setdefault('gibco', {}).setdefault(chr, {})

        tmp_venn_set, tmp_venn_counts = setops.venn_from_arrays(
            test_results[pid]['matched'][chr][cls],
            test_results[pid]['gibco'][chr][cls],
        )

        test_results_exclusive[pid]['matched'][chr][cls] = dict([
            (t, test_results[pid]['matched'][chr][cls][t]) for t in tmp_venn_set['10']
        ])
        test_results_exclusive[pid]['gibco'][chr][cls] = dict([
            (t, test_results[pid]['gibco'][chr][cls][t]) for t in tmp_venn_set['01']
        ])

        test_results_inclusive[pid]['matched'][chr][cls] = dict([
            (t, test_results[pid]['matched'][chr][cls][t]) for t in tmp_venn_set['11']
        ])
        test_results_inclusive[pid]['gibco'][chr][cls] = dict([
            (t, test_results[pid]['gibco'][chr][cls][t]) for t in tmp_venn_set['11']
        ])

    return {
        'exclusive': test_results_exclusive,
        'inclusive': test_results_inclusive,
    }


def n_overlap_datum(joint_dmr_de, this_de, this_dmr):
    n_de = this_de.shape[0]
    n_dmr = len(list(dmr.dict_iterator(this_dmr, n_level=3)))
    n_dmr_genes = dmr.count_dmr_genes(this_dmr)
    n_overlaps = joint_dmr_de['all'].shape[0]
    n_overlaps_unique = joint_dmr_de['all'].me_genes.unique().shape[0]
    this_datum = [n_de, n_dmr, n_dmr_genes, n_overlaps, n_overlaps_unique]
    for cls in dmr.CLASSES:
        this_datum.append(joint_dmr_de[cls].shape[0])
        this_datum.append(joint_dmr_de[cls].me_genes.unique().shape[0])
    return this_datum


def joint_de_dmr_counts(joint_de_dmr_result, de_result, dmr_result):
    the_cols = ['DE genes', 'DMR', 'DMR genes', 'overlaps', 'unique overlaps']
    the_cols += reduce(lambda x, y: x + y, [['%s' % t, '%s_unique' % t] for t in dmr.CLASSES], [])
    joint_count_table = pd.DataFrame(
        columns=the_cols,
        index=pd.Index(joint_de_dmr_result.keys(), name='patient'),
    )
    for pid in joint_de_dmr_result:
        joint_count_table.loc[pid] = n_overlap_datum(joint_de_dmr_result[pid], de_result[pid], dmr_result[pid])
    return joint_count_table


def joint_de_dmr_table(joint_de_dmr_result):

    de_cols = [
        ('Gene Symbol', 'DE gene symbol'),
        ('logFC', 'DE logFC'),
        ('logCPM', 'DE logCPM'),
        ('PValue', 'DE pvalue'),
        ('FDR', 'DE FDR'),
    ]

    dmr_cols = [
        'DMR clusterID',
        'DMR class tss',
        'DMR class gene',
        'DMR class island',
        'DMR chr',
        'DMR median delta M',
        'DMR padj',
    ]

    dmr_col_map = [
        ('DMR clusterID', 'me_cid'),
        ('DMR chr', 'chr'),
        ('DMR median delta M', 'me_mediandelta'),
        ('DMR padj', 'me_fdr'),
    ]
    dmr_class_series = pd.Series(0, index=['DMR class %s' % t for t in dmr.CLASSES])
    direction_series = pd.Series('-', index=['DE direction, DMR direction'])

    all_cols = [t[1] for t in de_cols] + dmr_cols + ['DE direction', 'DMR direction']

    out = {}

    for pid in joint_de_dmr_result:
        df = pd.DataFrame(columns=all_cols)
        row_lookup = {}
        i = 0
        for cls in dmr.CLASSES:
            this_tbl = joint_de_dmr_result[pid][cls]
            for _, row in this_tbl.iterrows():
                hsh = (row.loc['Gene Symbol'], row.loc['chr'], row.loc['me_cid'])
                if hsh in row_lookup:
                    # already seen this precise DE / DMR combination
                    # we only need to update the class types
                    df.loc[row_lookup[hsh], 'DMR class %s' % cls] = 1
                else:
                    new_row_de = pd.Series(
                        row.loc[list(zip(*de_cols)[0])].values,
                        index=list(zip(*de_cols)[1])
                    )

                    new_row_dmr = pd.Series(
                        row.loc[list(zip(*dmr_col_map)[1])].values,
                        index=list(zip(*dmr_col_map)[0])
                    )
                    new_row_dmr = new_row_dmr.append(dmr_class_series)
                    new_row_dmr.loc['DMR class %s' % cls] = 1

                    new_row = pd.concat((new_row_de, new_row_dmr))
                    new_row = new_row.append(direction_series)
                    new_row.loc['DE direction'] = 'U' if float(row.loc['logFC']) >= 0 else 'D'
                    new_row.loc['DMR direction'] = 'U' if float(row.loc['me_mediandelta']) >= 0 else 'D'
                    df.loc[i] = new_row

                    row_lookup[hsh] = i
                    i += 1

        # sort by abs DE logFC
        idx = np.abs(df.loc[:, 'DE logFC']).sort_values(ascending=False).index
        out[pid] = df.loc[idx]

    return out


def results_to_excel(blocks, fn):
    """

    :param blocks: Dictionary containing the different comparisons to save. Values are pandas dataframes.
    :param fn: Output file
    :return:
    """
    xl_writer = pd.ExcelWriter(outfile)
    # sort the keys for a more sensible order
    keys = sorted(blocks.keys())
    for k in keys:
        bl = blocks[k]
        bl.to_excel(xl_writer, k)
    xl_writer.save()


def results_to_ipa_format(blocks, outdir, incl_cols=('logFC', 'FDR')):
    incl_cols = list(incl_cols)

    for k, bl in blocks.iteritems():
        fn = os.path.join(outdir, "%s.txt" % k)
        header = [
            ['Key', 'Value'],
            ['identifier_types', 'Ensembl'],
            ['observation_name', k],
            ['date_created', datetime.datetime.now().isoformat()]
        ]
        with open(fn, 'wb') as f:
            c = csv.writer(f, delimiter='\t')
            # meta header
            c.writerows(header)
            c.writerow(['Data_begins_here'])
            # data column header
            c.writerow(['ID'] + incl_cols)
            # reduced block
            reduced_block = bl.loc[:, incl_cols]
            c.writerows(reduced_block.itertuples())



if __name__ == "__main__":
    # if this is specified, we load the DMR results from a JSON rather than recomputing them to save time
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'integrate_rnaseq_methylation')

    outdir = output.unique_output_dir("integrate_rnaseq_methylation")
    ref_name = 'GIBCONSC_P4'
    pids = ['017', '050', '054', '061']

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'GLM'
    }

    dmr_params = {
        'core_min_sample_overlap': 3,  # 3 / 4 samples must match
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'fdr': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': 8,
    }

    # Load RNA-Seq
    rnaseq_obj = rnaseq_data.load_by_patient(pids, annotate_by='Ensembl Gene ID')
    # discard unmapped, etc
    rnaseq_obj.data = rnaseq_obj.data.loc[rnaseq_obj.data.index.str.contains('ENSG')]
    # filter to remove genes that are not expressed sufficiently
    # this is only used in the group comparison - individual samples are filtered separately in batches
    # rnaseq_dat_filt = filter_by_cpm(rnaseq_obj.data, min_n_samples=2)

    de_res = compute_de(rnaseq_obj, pids, **de_params)
    de = de_res['de']
    de_matched = de_res['de_matched']
    de_ref = de_res['de_gibco']

    # split results by Venn group
    de_matched_only = dict([
        (pid, de_matched[pid].loc[de[pid]['10']]) for pid in pids
    ])
    de_ref_only = dict([
                           (pid, de_ref[pid].loc[de[pid]['01']]) for pid in pids
                           ])
    # DE intersection: use the matched comparison values
    de_intersection_matched = dict([
        (pid, de_matched[pid].loc[de[pid]['11']]) for pid in pids
    ])
    # DE intersection: use the ref comparison values
    de_intersection_ref = dict([
                                   (pid, de_ref[pid].loc[de[pid]['11']]) for pid in pids
                                   ])

    # write DE results to disk
    blocks = {}
    for pid in pids:
        blocks['GBM%s_pair_all' % pid] = de_matched[pid]
        blocks['GBM%s_pair_only' % pid] = de_matched[pid].loc[de[pid]['10']]
        blocks['GBM%s_ref_all' % pid] = de_ref[pid]
        blocks['GBM%s_ref_only' % pid] = de_ref[pid].loc[de[pid]['01']]
        # intersection: both results
        blocks['GBM%s_pair_and_ref' % pid] = de_matched[pid].loc[de[pid]['11']]
        blocks['GBM%s_ref_and_pair' % pid] = de_ref[pid].loc[de[pid]['11']]

    # single Excel workbook
    outfile = os.path.join(outdir, 'individual_gene_lists_de.xlsx')
    results_to_excel(blocks, outfile)

    # IPA lists, one per file
    ipa_outdir = os.path.join(outdir, 'ipa_de')
    if not os.path.exists(ipa_outdir):
        os.makedirs(ipa_outdir)

    results_to_ipa_format(blocks, ipa_outdir)

    # Load DNA Methylation
    me_data, me_meta = methylation_array.load_by_patient(pids)
    me_data.dropna(inplace=True)
    me_data = process.m_from_beta(me_data)
    anno = methylation_array.load_illumina_methylationepic_annotation()

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(me_data.index)
    anno = anno.loc[common_probes]
    dmr.add_merged_probe_classes(anno)
    me_data = me_data.loc[common_probes]

    # Compute DMR
    loaded = False
    if DMR_LOAD_DIR is not None:
        fn_in = os.path.join(DMR_LOAD_DIR, 'dmr_results.json')
        if os.path.isfile(fn_in):
            with open(fn_in, 'rb') as f:
                test_results = json.load(f)
            loaded = True
        # recompute relevant and significant results
        test_results_relevant = dmr.filter_dictionary(
            test_results,
            lambda x: 'pval' in x,
            n_level=5,
        )
        test_results_significant = dmr.filter_dictionary(
            test_results,
            lambda x: 'pval' in x and x['rej_h0'],
            n_level=5,
        )

    if not loaded:
        dmr_res = compute_dmr(me_data, me_meta, anno, pids, dmr_params)
        test_results = dmr_res['results']
        test_results_relevant = dmr_res['relevant']
        test_results_significant = dmr_res['significant']

        # Save DMR results to disk
        fout = os.path.join(outdir, "dmr_results.json")
        with open(fout, 'wb') as f:
            json.dump(test_results, f, cls=TestResultEncoder)
        print "Saved DMR results to %s" % fout

    tmp_venn_set = dmr_venn_sets(test_results_significant)
    test_results_exclusive = tmp_venn_set['exclusive']
    test_results_inclusive = tmp_venn_set['inclusive']

    # summarise DMR results in a table
    table_cluster_numbers = tabulate_dmr_results(
        test_results,
        relevant=test_results_relevant,
        significant=test_results_significant,
        significant_exclusive=test_results_exclusive,
        significant_inclusive=test_results_inclusive,
        outdir=outdir
    )

    # plot the DMR counts, classified into hyper and hypomethylation
    venn_dmr_counts(test_results_significant, outdir=outdir, figname="dmr_venn_all")
    for cls in dmr.CLASSES:
        venn_dmr_counts(test_results_significant, probe_class=cls, outdir=outdir, figname="dmr_venn_%s" % cls)
    dmr_overlap(test_results_significant, outdir=outdir, figname="dmr_overlap_all")
    for cls in dmr.CLASSES:
        dmr_overlap(test_results_significant, probe_class=cls, outdir=outdir, figname="dmr_overlap_%s" % cls)

    # side story: check overlap between individuals and group
    # skip this for now but useful to know!

    # integrate the two results
    dmr_matched = dmr.dict_by_sublevel(test_results_significant, 2, 'matched')
    dmr_matched_only = dmr.dict_by_sublevel(test_results_exclusive, 2, 'matched')
    dmr_ref = dmr.dict_by_sublevel(test_results_significant, 2, 'gibco')
    dmr_ref_only = dmr.dict_by_sublevel(test_results_exclusive, 2, 'gibco')
    dmr_intersection = dmr.dict_by_sublevel(test_results_inclusive, 2, 'matched')

    joint_de_dmr = {
        'pair_all': compute_joint_de_dmr(dmr_matched, de_matched),
        'pair_only': compute_joint_de_dmr(dmr_matched_only, de_matched_only),
        'pair_and_ref': compute_joint_de_dmr(dmr_intersection, de_intersection_matched),
        'ref_and_pair': compute_joint_de_dmr(dmr_intersection, de_intersection_ref),
        'ref_all': compute_joint_de_dmr(dmr_ref, de_ref),
        'ref_only': compute_joint_de_dmr(dmr_ref_only, de_ref_only),
    }

    joint_de_dmr_counts = {
        'pair_all': joint_de_dmr_counts(joint_de_dmr['pair_all'], de_matched, dmr_matched),
        'pair_only': joint_de_dmr_counts(
            joint_de_dmr['pair_only'],
            de_matched_only,
            dmr_matched_only
        ),
        'pair_and_ref': joint_de_dmr_counts(
            joint_de_dmr['pair_and_ref'],
            de_intersection_matched,
            dmr_intersection
        ),
    }


    ############

    # write DE / DMR results to disk
    # 1. Excel (human-readable)
    blocks = {}

    for k, v in joint_de_dmr.items():
        this_tbl = joint_de_dmr_table(v)
        for pid in this_tbl:
            blocks['GBM%s_%s' % (pid, k)] = this_tbl[pid]

    outfile = os.path.join(outdir, 'individual_gene_lists_de_dmr.xlsx')
    results_to_excel(blocks, outfile)

    # 2. IPA format

    blocks = {}
    for pid in pids:
        for k, v in joint_de_dmr.items():
            for cls in v[pid]:
                blocks["GBM%s_%s_%s" % (pid, k, cls)] = v[pid][cls]

    # IPA lists, one per file
    ipa_outdir = os.path.join(outdir, 'ipa_dmr')
    if not os.path.exists(ipa_outdir):
        os.makedirs(ipa_outdir)

    incl_cols = ['logFC', 'FDR']

    for k, bl in blocks.iteritems():
        fn = os.path.join(outdir, "%s.txt" % k)
        header = [
            ['Key', 'Value'],
            ['identifier_types', 'Ensembl'],
            ['observation_name', k],
            ['date_created', datetime.datetime.now().isoformat()]
        ]
        with open(fn, 'wb') as f:
            c = csv.writer(f, delimiter='\t')
            # meta header
            c.writerows(header)
            c.writerow(['Data_begins_here'])
            # data column header
            c.writerow(['ID'] + incl_cols)
            # reduced block
            reduced_block = bl.loc[:, incl_cols]
            c.writerows(reduced_block.itertuples())





    # results_to_ipa_format(blocks, ipa_outdir)

    #################




    ## TODO: move this elsewhere
    aa = joint_de_dmr['pair_only']['054']['tss'].set_index('Gene Symbol')
    x = aa.logFC
    y = aa.me_mediandelta.astype(float)
    r = (x ** 2 + y ** 2) ** .5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    idx_ul = ((x < 0) & (y > 0))
    ax.scatter(x.loc[idx_ul], y.loc[idx_ul], c='r')
    idx_br = ((x > 0) & (y < 0))
    ax.scatter(x.loc[idx_br], y.loc[idx_br], c='b')
    idx_od = ((x > 0) & (y > 0)) | ((x < 0) & (y < 0))
    ax.scatter(x.loc[idx_od], y.loc[idx_od], c='gray')

    for i in np.where(idx_ul & (r > 8))[0]:
        ax.text(x.iloc[i], y.iloc[i], x.index[i], color='r')
    for i in np.where(idx_br & (r > 5))[0]:
        ax.text(x.iloc[i], y.iloc[i], x.index[i], color='b')
    ax.axhline(0, ls='--', c='k', alpha=0.4)
    ax.axvline(0, ls='--', c='k', alpha=0.4)
    ax.set_xlabel("DE logFC")
    ax.set_ylabel("DMR median delta")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "integrated_scatter_labelled_054_tss_pair_only.png"), dpi=200)

    aa = joint_de_dmr['pair_and_ref']['054']['tss'].set_index('Gene Symbol')
    x = aa.logFC
    y = aa.me_mediandelta.astype(float)
    r = (x ** 2 + y ** 2) ** .5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    idx_ul = ((x < 0) & (y > 0))
    ax.scatter(x.loc[idx_ul], y.loc[idx_ul], c='r')
    idx_br = ((x > 0) & (y < 0))
    ax.scatter(x.loc[idx_br], y.loc[idx_br], c='b')
    idx_od = ((x > 0) & (y > 0)) | ((x < 0) & (y < 0))
    ax.scatter(x.loc[idx_od], y.loc[idx_od], c='gray')

    for i in np.where(idx_ul & (r > 12))[0]:
        ax.text(x.iloc[i], y.iloc[i], x.index[i], color='r')
    for i in np.where(idx_br & (r > 8))[0]:
        ax.text(x.iloc[i], y.iloc[i], x.index[i], color='b')
    ax.axhline(0, ls='--', c='k', alpha=0.4)
    ax.axvline(0, ls='--', c='k', alpha=0.4)
    ax.set_xlabel("DE logFC")
    ax.set_ylabel("DMR median delta")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "integrated_scatter_labelled_054_tss_pair_and_ref.png"), dpi=200)

    # export integrated results to lists
    # tmp_for_xl = {}
    #
    # for k, v in joint_de_dmr.items():
    #     this_tbl = joint_de_dmr_table(v)
    #     for pid in this_tbl:
    #         tmp_for_xl.setdefault(pid, {})[k] = this_tbl[pid]
    #
    # xl_writer = pd.ExcelWriter(os.path.join(outdir, "individual_gene_lists_de_dmr.xlsx"))
    # for pid in tmp_for_xl:
    #     for k in tmp_for_xl[pid]:
    #         sheet_name = 'GBM%s_%s' % (pid, k)
    #         tmp_for_xl[pid][k].to_excel(xl_writer, sheet_name, index=False)
    # xl_writer.save()

    # TODO: move to plots
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for pid in pids:
        for i, cmp in enumerate(joint_de_dmr.keys()):
            for j, cls in enumerate(joint_de_dmr[cmp][pid].keys()):
                # iterate over genes, get quadrant, define colour, plot text
                pass




    # def venn_diagram_and_core_genes(meth_de, text_file, fig_file, min_overlap=4, fig_title=None, plot_genes=None):
    #     """
    #
    #     :param meth_de:
    #     :param text_file:
    #     :param fig_file:
    #     :param min_overlap:
    #     :param fig_title:
    #     :param plot_genes: None, 'scatter', 'colours'
    #     :return:
    #     """
    #     n_sample = len(meth_de)
    #     # one of the samples is 'all', which we ignore
    #     bns = list(setops.binary_combinations_sum_gte(n_sample - 1, min_overlap))
    #
    #     fig_kws = {'figsize': (8, 3.2)}
    #     # if fig_title is not None:
    #     #     fig_kws['num'] = fig_title
    #     fig, axs = plt.subplots(ncols=3, num=fig_title, **fig_kws)
    #
    #     f = open(text_file, 'wb')
    #     venn_counts = {}
    #     venn_sets = {}
    #     core_sets = {}
    #
    #     if plot_genes:
    #         fig_gene_scatter, axs_gene_scatter = plt.subplots(nrows=1, ncols=3, **fig_kws)
    #         fig_gene_colours, axs_gene_colours = plt.subplots(nrows=1, ncols=3, figsize=(5, 8))
    #         max_y = 0
    #         max_de = 0.
    #         max_dmr = 0.
    #
    #     # one Venn diagram per class
    #     # Venn components are patients
    #     # we therefore build a structure containing patient-wise gene sets for each class
    #     for i, cls in enumerate(dmr.CLASSES):
    #         ax = axs[i]
    #         # reordered dict by sublevel
    #         tmp = dmr.dict_by_sublevel(meth_de, 2, cls)
    #         # build labels and gene lists to control the order (otherwise we're relying on dict iteration order)
    #         # we skip over 'all' here
    #         set_labels = []
    #         gene_sets = []
    #         for k, t in tmp.items():
    #             if k != 'all':
    #                 gene_sets.append(set(t.genes))
    #                 set_labels.append(k)
    #         venn_res, venn_sets[cls], venn_counts[cls] = venn.venn_diagram(
    #             *gene_sets,
    #             set_labels=set_labels,
    #             ax=ax
    #         )
    #         ax.set_title(cls)
    #
    #         # Core genes are defined as those supported by min_overlap samples
    #         core_genes = set()
    #         for bn in bns:
    #             core_genes = core_genes.union(venn_sets[cls][bn])
    #         core_sets[cls] = core_genes
    #         core_genes = sorted(core_genes)
    #         print "%s core genes [%d]: %s" % (cls, len(core_genes), ', '.join(core_genes))
    #         f.write("%s core genes [%d]: %s\n" % (cls, len(core_genes), ', '.join(core_genes)))
    #
    #         if plot_genes:
    #             j = 0
    #             de_vals = []
    #             dmr_vals = []
    #             for cg in core_genes:
    #                 val_de = 0.
    #                 val_dmr = 0.
    #                 n = 0.
    #                 # loop over patients
    #                 for t in tmp.values():
    #                     if (t.genes == cg).any():
    #                         n += (t.genes == cg).sum()
    #                         the_rows = t.loc[t.genes == cg]
    #                         val_de += the_rows.logFC.astype(float).sum()
    #                         val_dmr += the_rows.me_mediandelta.astype(float).sum()
    #                 x = val_de / n
    #                 y = val_dmr / n
    #                 de_vals.append(x)
    #                 dmr_vals.append(y)
    #                 axs_gene_scatter[i].text(val_de / n, val_dmr / n, cg)
    #
    #                 if x > 0 and y > 0:
    #                     cmap = plt.cm.Blues
    #                 elif x < 0 and y < 0:
    #                     cmap = plt.cm.Purples
    #                 elif x > 0 and y < 0:
    #                     cmap = plt.cm.Reds
    #                 elif x < 0 and y > 0:
    #                     cmap = plt.cm.Greens
    #                 else:
    #                     cmap = plt.cm.Greys
    #
    #                 # scaling
    #                 ## FIXME: this would be better applied afterwards, like the axis scaling
    #                 xn = x / 12.
    #                 yn = y / 6.
    #                 rn = (xn ** 2 + yn ** 2) ** .5
    #                 rn = min(max(rn, 0.3), 0.9)
    #                 c = cmap(rn)
    #                 axs_gene_colours[i].text(0, j, cg, color=c)
    #                 j += 1.
    #                 axs_gene_colours[i].axis('off')
    #                 axs_gene_colours[i].set_title(cls)
    #
    #             axs_gene_scatter[i].scatter(de_vals, dmr_vals)
    #             axs_gene_scatter[i].axhline(0.)
    #             axs_gene_scatter[i].axvline(0.)
    #             axs_gene_scatter[i].set_xlabel('RNASeq DE logFC')
    #             axs_gene_scatter[i].set_ylabel('EPIC DMR median delta M')
    #
    #             max_de = max(max_de, np.abs(de_vals).max())
    #             max_dmr = max(max_dmr, np.abs(dmr_vals).max())
    #             max_y = max(max_y, j)
    #
    #     if plot_genes:
    #         for i in range(len(dmr.CLASSES)):
    #             axs_gene_colours[i].set_ylim([0, max_y])
    #             axs_gene_scatter[i].set_xlim([-max_de * 1.2, max_de * 1.2])
    #             axs_gene_scatter[i].set_ylim([-max_dmr * 1.2, max_dmr * 1.2])
    #
    #     # intersection of all
    #     core_all = sorted(reduce(set.intersection, core_sets.values()))
    #     print "Core genes shared across all classes [%d]: %s" % (len(core_all), ', '.join(core_all))
    #     f.write("Core genes shared across all classes [%d]: %s\n" % (len(core_all), ', '.join(core_all)))
    #
    #     # union of all
    #     core_union = sorted(reduce(set.union, core_sets.values()))
    #     print "Core genes in >=1 class [%d]: %s" % (len(core_union), ', '.join(core_union))
    #     f.write("Core genes in >= 1 class [%d]: %s\n" % (len(core_union), ', '.join(core_union)))
    #     f.close()
    #
    #     fig.tight_layout()
    #     fig.savefig("%s.png" % fig_file, dpi=200)
    #     fig.savefig("%s.pdf" % fig_file)
    #
    #     if plot_genes:
    #         fig_gene_scatter.tight_layout()
    #         fig_gene_scatter.savefig("%s_genescatter.png" % fig_file, dpi=200)
    #         fig_gene_scatter.savefig("%s_genescatter.pdf" % fig_file)
    #
    #         fig_gene_colours.tight_layout()
    #         fig_gene_colours.savefig("%s_genecolours.png" % fig_file, dpi=200)
    #         fig_gene_colours.savefig("%s_genecolours.pdf" % fig_file)
    #
    #     return venn_sets
