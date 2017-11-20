from load_data import rnaseq_data, methylation_array
from rnaseq import differential_expression
from rnaseq.filter import filter_by_cpm
from methylation import process, dmr
from methylation import plots as me_plots
from integrator.rnaseq_methylationarray import compute_joint_de_dmr
from integrator import plots
import copy
import pandas as pd
import numpy as np
from scipy import stats
import references
import os
import itertools
import csv
import datetime
import json
from utils import output, setops, dictionary
from matplotlib import pyplot as plt
import seaborn as sns


class BasicLogicException(Exception):
    pass


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


def compute_dmr2(me_data, me_meta, anno, pids, dmr_params, ref_name_contains='GIBCO'):

    obj = dmr.DmrResults(anno=anno)
    obj.identify_clusters(**dmr_params)
    res = {}

    for pid in pids:
        if pid == 'all':
            the_idx = ~me_meta.index.str.contains(ref_name_contains)
        else:
            the_idx = me_meta.index.str.contains(pid)
        the_groups = me_meta.loc[the_idx, 'type'].values
        the_samples = me_meta.index[the_idx].groupby(the_groups).values()

        the_obj = obj.copy()
        the_obj.test_clusters(me_data,
                              samples=the_samples,
                              n_jobs=dmr_params['n_jobs'],
                              min_median_change=dmr_params['delta_m_min'],
                              method=dmr_params['dmr_test_method'],
                              **dmr_params['test_kwargs']
                              )
        res.setdefault(pid, {})['matched'] = the_obj

        if pid == 'all':
            the_idx = (me_meta.type != 'iNSC')
        else:
            the_idx = (me_meta.index.str.contains(pid) & (me_meta.type == 'GBM')) | me_meta.index.str.contains(ref_name_contains)
        the_groups = me_meta.loc[the_idx, 'type'].values
        the_samples = me_meta.index[the_idx].groupby(the_groups).values()

        the_obj = obj.copy()
        the_obj.test_clusters(me_data,
                              samples=the_samples,
                              n_jobs=dmr_params['n_jobs'],
                              min_median_change=dmr_params['delta_m_min'],
                              method=dmr_params['dmr_test_method'],
                              **dmr_params['test_kwargs']
                              )
        res.setdefault(pid, {})['ref'] = the_obj

    return dmr.DmrResultCollection(**res)


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
        test_results[pid]['matched'] = dmr.test_clusters_in_place(
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
        test_results[pid]['gibco'] = dmr.test_clusters_in_place(
            clusters,
            me_data,
            samples=the_samples,
            n_jobs=dmr_params['n_jobs'],
            min_median_change=dmr_params['delta_m_min'],
            method=dmr_params['dmr_test_method'],
            **dmr_params['test_kwargs']
        )

        for typ in ['matched', 'gibco']:
            dmr.mht_correction(
                test_results[pid][typ],
                alpha=dmr_params['fdr']
            )
            for k, v in test_results[pid][typ].iteritems():
                test_results_relevant[pid].setdefault(typ, {})[k] = dictionary.filter_dictionary(
                    v,
                    filt=lambda x: 'pval' in x,
                    n_level=1
                )
                test_results_significant[pid].setdefault(typ, {})[k] = dictionary.filter_dictionary(
                    v,
                    filt=lambda x: x.get('rej_h0', False),
                    n_level=1
                )

    return {
        'clusters': clusters,
        'results': test_results,
        'relevant': test_results_relevant,
        'significant': test_results_significant,
    }


def save_dmr_results(dmr_results, fn):
    with open(fn, 'wb') as f:
        json.dump(dmr_results['all'], f, cls=TestResultEncoder)


def load_dmr_results(fn):
    test_results = {}
    with open(fn, 'wb') as f:
        test_results['all'] = json.load(f)

    # for convenience, provide nested dictionaries pointing to specific class types
    test_results = dmr.cluster_class_nests(test_results['all'])

    test_results_relevant = dmr.filter_dictionary(
        test_results,
        lambda x: 'pval' in x,
        n_level=2,
    )
    test_results_significant = dmr.filter_dictionary(
        test_results,
        lambda x: 'pval' in x and x['rej_h0'],
        n_level=2,
    )
    return {
        'test_results': test_results,
        'test_results_relevant': test_results_relevant,
        'test_results_significant': test_results_significant
    }


def tabulate_dmr_results(dmr_res, outdir=None, comparisons = ('matched', 'ref')):
    for_counting = {
        'relevant': dmr_res.clusters_relevant,
        'significant': dmr_res.clusters_significant,
    }

    pids = dmr_res.keys()
    the_keys = ['proposed', 'relevant', 'significant', 'significant_exclusive', 'significant_intersection']
    cols = ['sample'] + ['clusters_%s' % t for t in the_keys] + ['genes_%s' % t for t in the_keys]
    table_cluster_numbers = {}

    # the number of proposed clusters and proposed genes is constant across all results, so just use the first
    ncl = len(dmr_res.clusters)
    ng = dmr.count_dmr_genes(dmr_res.clusters)
    row_template = {'clusters_proposed': ncl, 'genes_proposed': ng}

    for pid in pids:

        tmp, _ = setops.venn_from_arrays(*[dmr_res[pid][typ].clusters_significant.keys() for typ in comparisons])

        for typ in comparisons:
            the_venn_key = '10' if typ == 'matched' else '01'
            po = dict([(k, for_counting['significant'][pid][typ][k]) for k in tmp[the_venn_key]])
            pr = dict([(k, for_counting['significant'][pid][typ][k]) for k in tmp['11']])
            table_cluster_numbers.setdefault(typ, pd.DataFrame(columns=cols))
            this_row = dict(row_template)
            this_row.update({
                'sample': pid,
                'clusters_relevant': len(for_counting['relevant'][pid][typ]),
                'clusters_significant': len(for_counting['significant'][pid][typ]),
                'clusters_significant_exclusive': len(po),
                'clusters_significant_intersection': len(pr),
                'genes_relevant': dmr.count_dmr_genes(for_counting['relevant'][pid][typ]),
                'genes_significant': dmr.count_dmr_genes(for_counting['significant'][pid][typ]),
                'genes_significant_exclusive': dmr.count_dmr_genes(po),
                'genes_significant_intersection': dmr.count_dmr_genes(pr),
            })
            this_row = pd.Series(this_row)
            table_cluster_numbers[typ] = table_cluster_numbers[typ].append(this_row, ignore_index=True)

    if outdir is not None:
        for typ in comparisons:
            table_cluster_numbers[typ].to_csv(os.path.join(outdir, "cluster_numbers.%s.csv" % typ))

    return table_cluster_numbers



def dmr_venn_sets(test_results):
    test_results_exclusive = {}
    test_results_inclusive = {}

    for pid, v in test_results.iteritems():
        test_results_exclusive.setdefault(pid, {}).setdefault('matched', {})
        test_results_exclusive.setdefault(pid, {}).setdefault('gibco', {})
        test_results_inclusive.setdefault(pid, {}).setdefault('matched', {})
        test_results_inclusive.setdefault(pid, {}).setdefault('gibco', {})

        tmp_venn_set, tmp_venn_counts = setops.venn_from_arrays(
            test_results[pid]['matched']['all'],
            test_results[pid]['gibco']['all'],
        )

        test_results_exclusive[pid]['matched'] = dmr.cluster_class_nests(
            dict([
                (t, test_results[pid]['matched']['all'][t]) for t in tmp_venn_set['10']
            ])
        )
        test_results_exclusive[pid]['gibco'] = dmr.cluster_class_nests(
            dict([
                (t, test_results[pid]['gibco']['all'][t]) for t in tmp_venn_set['01']
            ])
        )

        test_results_inclusive[pid]['matched'] = dmr.cluster_class_nests(
            dict([
                (t, test_results[pid]['matched']['all'][t]) for t in tmp_venn_set['11']
            ])
        )
        test_results_inclusive[pid]['gibco'] = dmr.cluster_class_nests(
            dict([
                (t, test_results[pid]['gibco']['all'][t]) for t in tmp_venn_set['11']
            ])
        )

    return {
        'exclusive': test_results_exclusive,
        'inclusive': test_results_inclusive,
    }


def n_overlap_datum(joint_dmr_de, this_de, this_dmr):
    n_de = this_de.shape[0]
    n_dmr = len(list(dictionary.dict_iterator(this_dmr, n_level=3)))
    n_dmr_genes = dmr.count_dmr_genes(this_dmr)
    n_overlaps = joint_dmr_de['all'].shape[0]
    n_overlaps_unique = joint_dmr_de['all'].loc[:, 'Gene Symbol'].unique().shape[0]
    this_datum = [n_de, n_dmr, n_dmr_genes, n_overlaps, n_overlaps_unique]
    for cls in dmr.CLASSES:
        this_datum.append(joint_dmr_de[cls].shape[0])
        this_datum.append(joint_dmr_de[cls].loc[:, 'Gene Symbol'].unique().shape[0])
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


def de_dmr_hash(dat, cols=('Gene Symbol', 'chr', 'me_cid')):
    """
    Generate a list of hash values from the supplied DataFrame
    :param dat:
    :param cols:
    :return:
    """
    return list(dat.loc[:, cols].itertuples(index=False, name=None))


def de_dmr_concordance_status(dat):
    """
    For each entry supplied, compute the concordance status (boolean)
    :param dat:
    :return:
    """
    return np.sign(dat.logFC.astype(float)) != np.sign(dat.me_mediandelta.astype(float))


def joint_de_dmr_table2(joint_de_dmr_result, associated_result=None):
    """
    Second version. This includes a column marking the concordancy of DE and DMR
    :param joint_de_dmr_result: Nested dictionary, as follows: PID -> DMR class -> concordancy_status
    e.g. ['054']['island']['concordant']
    This is achieved by extracting a Venn set before calling the function (e.g. 'pair_only').
    In practice, we can skip over concordancy status and use the column de_dmr_concordant to summarise
    :return: Dict of tables keyed by PID, all classes and concordancies summarised in columns
    """

    de_cols = [
        ('Gene Symbol', 'DE gene symbol'),
        ('logFC', 'DE logFC'),
        ('logCPM', 'DE logCPM'),
        ('PValue', 'DE pvalue'),
        ('FDR', 'DE FDR'),
    ]

    de_dtypes = [
        'object',
        'float',
        'float',
        'float',
        'float',
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

    dmr_dtypes = [
        'uint16',
        'uint8',
        'uint8',
        'uint8',
        'object',
        'float',
        'float',
    ]

    dmr_col_map = [
        ('DMR clusterID', 'me_cid'),
        ('DMR chr', 'chr'),
        ('DMR median delta M', 'me_mediandelta'),
        ('DMR padj', 'me_fdr'),
    ]

    de_dmr_cols = [
        'DE direction',
        'DMR direction',
        'DE DMR concordant',
    ]

    de_dmr_dtypes = [
        'object',
        'object',
        'object',
    ]

    dmr_class_series = pd.Series(0, index=['DMR class %s' % t for t in dmr.CLASSES])
    direction_series = pd.Series('-', index=['DE direction, DMR direction'])

    all_cols = [t[1] for t in de_cols] + dmr_cols + de_dmr_cols
    all_dtypes = de_dtypes + dmr_dtypes + de_dmr_dtypes

    out = {}

    for pid in joint_de_dmr_result:
        df = pd.DataFrame(columns=all_cols)
        df = df.astype(dict(zip(all_cols, all_dtypes)))
        row_lookup = {}
        i = 0
        for cls in dmr.CLASSES:

            if associated_result is not None:
                assoc_tbl = associated_result[pid][cls]
                assoc_hsh = de_dmr_hash(assoc_tbl)
                assoc_status = [row.loc['de_dmr_concordant'] for _, row in assoc_tbl.iterrows()]
                assoc_conc = dict(zip(assoc_hsh, assoc_status))

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

                    # direction
                    new_row.loc['DE direction'] = 'U' if float(row.loc['logFC']) >= 0 else 'D'
                    new_row.loc['DMR direction'] = 'U' if float(row.loc['me_mediandelta']) >= 0 else 'D'

                    # concordance
                    c = ['T' if row.loc['de_dmr_concordant'] else 'F']
                    if associated_result is None:
                        # assume no matching result (e.g. pair_only)
                        pass
                    else:
                        # get the associated concordance (e.g. pair_and_ref)
                        c += ['T' if assoc_conc[hsh] else 'F']

                    new_row.loc['DE DMR concordant'] = ','.join(c)

                    df.loc[i] = new_row

                    row_lookup[hsh] = i
                    i += 1

        # sort by abs DE logFC
        idx = np.abs(df.loc[:, 'DE logFC']).sort_values(ascending=False).index
        out[pid] = df.loc[idx]

    return out



def results_to_excel(blocks, fn, write_index=True):
    """

    :param blocks: Dictionary containing the different comparisons to save. Values are pandas dataframes.
    :param fn: Output file
    :param write_index: If True (default) then the index is written as the first column. May wish to disable this if it
    has no meaning.
    :return:
    """
    xl_writer = pd.ExcelWriter(fn)
    # sort the keys for a more sensible order
    keys = sorted(blocks.keys())
    for k in keys:
        bl = blocks[k]
        bl.to_excel(xl_writer, k, index=write_index)
    xl_writer.save()


def results_to_ipa_format(
        blocks,
        outdir,
        incl_cols=('logFC', 'FDR'),
        identifier='Ensembl',
):
    incl_cols = list(incl_cols)

    for k, bl in blocks.iteritems():
        fn = os.path.join(outdir, "%s.txt" % k)
        header = [
            ['Key', 'Value'],
            ['observation_name', k],
            ['date_created', datetime.datetime.now().isoformat()]
        ]
        if identifier is not None:
            header += [['identifier_types', identifier]]
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
    # RTK II
    rtkii_pids = ['017', '050', '054']
    # all n=2 samples and RTK II samples
    pids = ['018', '044', '049', '050', '052', '054', '061']

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
    # We load pre-computed results if a JSON file with the correct filename is found
    # Otherwise this is written after computing the results
    loaded = False
    if DMR_LOAD_DIR is not None:
        fn_in = os.path.join(DMR_LOAD_DIR, 'dmr_results.pkl')
        if os.path.isfile(fn_in):
            dmr_res = dmr.DmrResultCollection.from_pickle(fn_in, anno=anno)
            loaded = True

    if not loaded:
        dmr_res = compute_dmr2(me_data, me_meta, anno, pids, dmr_params)
        # Save DMR results to disk
        fout = os.path.join(outdir, "dmr_results.pkl")
        dmr_res.to_pickle(fout, include_annotation=False)
        print "Saved DMR results to %s" % fout

    test_results_relevant = dmr_res.results_relevant_by_class()
    test_results_significant = dmr_res.results_significant_by_class()

    # Venn showing the distribution of probe classes amongst all regions (incl. not significantly DM)
    me_plots.dmr_cluster_count_by_class(dmr_res[pids[0]]['matched'].clusters_by_class(), outdir=outdir)
    # Array of Venns showing the same distribution but of *significant* clusters in pair_all and ref_all
    me_plots.dmr_cluster_count_array(test_results_significant, comparisons=('matched', 'ref'),
                                     comparison_labels=('Matched', 'Ref'), outdir=outdir)

    # summarise DMR results in a table
    table_cluster_numbers = tabulate_dmr_results(dmr_res, outdir=outdir)

    # plot the DMR counts, classified into hyper and hypomethylation
    me_plots.venn_dmr_counts(dmr_res.results_significant, comparisons=('matched', 'ref'), outdir=outdir, figname="dmr_venn_all")
    for cls in dmr.CLASSES:
        me_plots.venn_dmr_counts(
            dmr_res.results_significant_by_class(cls),
            comparisons=('matched', 'ref'),
            outdir=outdir,
            figname="dmr_venn_%s" % cls
        )

    # these will fail if we have too many individuals (a limit of the Venn diagram)
    try:
        me_plots.dmr_overlap(
            dmr_res.results_significant,
            comparisons=('matched', 'ref'),
            outdir=outdir,
            figname="dmr_overlap_all"
        )
        for cls in dmr.CLASSES:
            me_plots.dmr_overlap(dmr_res.results_significant_by_class(cls), outdir=outdir, figname="dmr_overlap_%s" % cls)
    except Exception as exc:
        print "Unable to produce DMR overlap plot: %s" % repr(exc)

    # side story: check DMR overlap between individuals and group
    # skip this for now but useful to know!

    # convert DMR results to tables
    dmr_tables = dmr_res.apply(lambda x: x.to_table(include='significant'))

    # integrate the two results
    # We first look for overlap, then generate Venn sets
    dmr_matched = dictionary.dict_by_sublevel(dmr_res, 2, 'matched')
    dmr_ref = dictionary.dict_by_sublevel(dmr_res, 2, 'ref')

    # TODO: update from here

    # TODO: move to a function
    # these are the two inputs and the output

    de_dmr_matched = compute_joint_de_dmr(dmr_matched, de_matched)
    de_dmr_ref = compute_joint_de_dmr(dmr_ref, de_ref)
    joint_de_dmr_sets = {}

    overlap_types = [
        'pair_all',
        'pair_only',
        'pair_and_ref',
        'ref_all',
        'ref_only',
        'ref_and_pair'
    ]

    # now compute the intersections
    # we do this separately for concordant and discordant DE / DMR matches

    for pid in pids:
        for cls in ['all'] + list(dmr.CLASSES):
            for ot in overlap_types:
                joint_de_dmr_sets.setdefault(ot, {}).setdefault(pid, {}).setdefault(cls, {})
            x_pair = de_dmr_matched[pid][cls]
            x_ref = de_dmr_ref[pid][cls]

            # split into discordant and concordant (and add a column to the original data)
            conc_pair = np.sign(x_pair.logFC) != np.sign(x_pair.me_mediandelta)
            conc_ref = np.sign(x_ref.logFC) != np.sign(x_ref.me_mediandelta)
            x_pair.loc[:, 'de_dmr_concordant'] = conc_pair
            x_ref.loc[:, 'de_dmr_concordant'] = conc_ref

            # define the minimum unique identifier
            de_dmr_id_pair = x_pair.loc[:, ['Gene Symbol', 'chr', 'me_cid']].apply(tuple, axis=1)
            de_dmr_id_ref = x_ref.loc[:, ['Gene Symbol', 'chr', 'me_cid']].apply(tuple, axis=1)

            v_all, ct_all = setops.venn_from_arrays(
                de_dmr_id_pair,
                de_dmr_id_ref,
            )

            # sanity checks
            if de_dmr_id_pair.isin(v_all['01']).any():
                raise BasicLogicException("de_dmr_id_pair.isin(v_all['01']).any()")

            if de_dmr_id_ref.isin(v_all['10']).any():
                raise BasicLogicException("de_dmr_id_ref.isin(v_all['10']).any()")

            if de_dmr_id_pair.isin(v_all['10']).sum() + de_dmr_id_pair.isin(v_all['11']).sum() != de_dmr_id_pair.size:
                raise BasicLogicException("de_dmr_id_pair.isin(v_all['10']).sum() + de_dmr_id_pair.isin(v_all['11']).sum() != de_dmr_id_pair.size")

            if de_dmr_id_ref.isin(v_all['01']).sum() + de_dmr_id_ref.isin(v_all['11']).sum() != de_dmr_id_ref.size:
                raise BasicLogicException("de_dmr_id_ref.isin(v_all['01']).sum() + de_dmr_id_ref.isin(v_all['11']).sum() != de_dmr_id_ref.size")

            # no need to reduce to unique results here because the venn operations use sets
            # concordant
            v_conc, ct_conc = setops.venn_from_arrays(
                de_dmr_id_pair[conc_pair],
                de_dmr_id_ref[conc_ref],
            )

            # discordant
            disc_pair = ~conc_pair
            disc_ref = ~conc_ref

            v_disc, ct_disc = setops.venn_from_arrays(
                de_dmr_id_pair[disc_pair],
                de_dmr_id_ref[disc_ref],
            )

            # finally, mixed (two ways to achieve this)
            v_mix = set(de_dmr_id_pair[conc_pair]).intersection(de_dmr_id_ref[disc_ref])
            v_mix = v_mix.union(
                set(de_dmr_id_pair[disc_pair]).intersection(de_dmr_id_ref[conc_ref])
            )

            # stitch tables back together
            # include all rows for main output
            # for IPA we will reduce to DE logFC only (because no other values needed)

            joint_de_dmr_sets['pair_all'][pid][cls]['all'] = x_pair
            joint_de_dmr_sets['pair_only'][pid][cls]['all'] = x_pair.loc[de_dmr_id_pair.isin(v_all['10'])]
            joint_de_dmr_sets['pair_and_ref'][pid][cls]['all'] = x_pair.loc[de_dmr_id_pair.isin(v_all['11'])]
            joint_de_dmr_sets['ref_all'][pid][cls]['all'] = x_ref
            joint_de_dmr_sets['ref_only'][pid][cls]['all'] = x_ref.loc[de_dmr_id_ref.isin(v_all['01'])]
            joint_de_dmr_sets['ref_and_pair'][pid][cls]['all'] = x_ref.loc[de_dmr_id_ref.isin(v_all['11'])]

            de_dmr_id_pair_all_conc = v_conc['10'] + v_conc['11']
            joint_de_dmr_sets['pair_all'][pid][cls]['concordant'] = x_pair.loc[de_dmr_id_pair.isin(v_conc['10'] + v_conc['11'])]
            joint_de_dmr_sets['pair_only'][pid][cls]['concordant'] = x_pair.loc[de_dmr_id_pair.isin(v_conc['10'])]
            joint_de_dmr_sets['pair_and_ref'][pid][cls]['concordant'] = x_pair.loc[de_dmr_id_pair.isin(v_conc['11'])]
            joint_de_dmr_sets['ref_all'][pid][cls]['concordant'] = x_ref.loc[de_dmr_id_ref.isin(v_conc['01'] + v_conc['11'])]
            joint_de_dmr_sets['ref_only'][pid][cls]['concordant'] = x_ref.loc[de_dmr_id_ref.isin(v_conc['01'])]
            joint_de_dmr_sets['ref_and_pair'][pid][cls]['concordant'] = x_ref.loc[de_dmr_id_ref.isin(v_conc['11'])]

            joint_de_dmr_sets['pair_all'][pid][cls]['discordant'] = x_pair.loc[de_dmr_id_pair.isin(v_disc['10'] + v_disc['11'])]
            joint_de_dmr_sets['pair_only'][pid][cls]['discordant'] = x_pair.loc[de_dmr_id_pair.isin(v_disc['10'])]
            joint_de_dmr_sets['pair_and_ref'][pid][cls]['discordant'] = x_pair.loc[de_dmr_id_pair.isin(v_disc['11'])]
            joint_de_dmr_sets['ref_all'][pid][cls]['discordant'] = x_ref.loc[de_dmr_id_ref.isin(v_disc['01'] + v_disc['11'])]
            joint_de_dmr_sets['ref_only'][pid][cls]['discordant'] = x_ref.loc[de_dmr_id_ref.isin(v_disc['01'])]
            joint_de_dmr_sets['ref_and_pair'][pid][cls]['discordant'] = x_ref.loc[de_dmr_id_ref.isin(v_disc['11'])]

            joint_de_dmr_sets['pair_and_ref'][pid][cls]['mixed'] = x_pair.loc[de_dmr_id_pair.isin(v_mix)]
            joint_de_dmr_sets['ref_and_pair'][pid][cls]['mixed'] = x_ref.loc[de_dmr_id_ref.isin(v_mix)]

            # sanity checks:
            for k, v in joint_de_dmr_sets.items():
                if not v[pid][cls]['concordant'].de_dmr_concordant.all():
                    raise BasicLogicException("concordant %s: not v.de_dmr_concordant.all()" % k)
                if v[pid][cls]['discordant'].de_dmr_concordant.any():
                    raise BasicLogicException("discordant %s: v.de_dmr_concordant.any()" % k)

            aa = set(joint_de_dmr_sets['ref_and_pair'][pid][cls]['mixed'].loc[:, 'Gene Symbol'])
            bb = set(joint_de_dmr_sets['pair_and_ref'][pid][cls]['mixed'].loc[:, 'Gene Symbol'])
            if aa != bb:
                raise BasicLogicException("IDs in mixed ref / pair do not match")

    # write DE / DMR results to disk
    # 1. Excel (human-readable)
    # For this, we ignore the concordance result and go through all results
    blocks = {}
    for typ in ['pair_all', 'pair_only', 'ref_all', 'ref_only']:
        t = joint_de_dmr_table2(dictionary.dict_by_sublevel(joint_de_dmr_sets[typ], 3, 'all'))
        for pid in t:
            blocks['GBM%s_%s' % (pid, typ)] = t[pid]

    t = joint_de_dmr_table2(
        dictionary.dict_by_sublevel(joint_de_dmr_sets['pair_and_ref'], 3, 'all'),
        associated_result=dictionary.dict_by_sublevel(joint_de_dmr_sets['ref_and_pair'], 3, 'all')
    )
    for pid in t:
        blocks['GBM%s_pair_and_ref' % pid] = t[pid]

    t = joint_de_dmr_table2(
        dictionary.dict_by_sublevel(joint_de_dmr_sets['ref_and_pair'], 3, 'all'),
        associated_result=dictionary.dict_by_sublevel(joint_de_dmr_sets['pair_and_ref'], 3, 'all')
    )
    for pid in t:
        blocks['GBM%s_ref_and_pair' % pid] = t[pid]

    outfile = os.path.join(outdir, 'individual_gene_lists_de_dmr.xlsx')
    results_to_excel(blocks, outfile, write_index=False)

    # 2. IPA format
    # Here, we reduce the output to genes.
    # A gene is concordant if >=1 corresponding DE/DMR rows are concordant
    ipa_outdir = os.path.join(outdir, 'ipa_dmr')
    if not os.path.exists(ipa_outdir):
        os.makedirs(ipa_outdir)

    # We want several different versions
    def prepare_de_dmr_concordant_blocks(dat, key='all'):
        incl_cols = ['logFC', 'FDR']
        blocks = {}
        for k, v in dat.items():
            for pid in v:
                bl = v[pid][key]
                # reduce to concordant results
                dmr_sign = np.sign(bl.loc[:, 'me_mediandelta'].astype(float))
                de_sign = np.sign(bl.loc[:, 'logFC'].astype(float))
                bl = bl.loc[dmr_sign != de_sign]
                # only keep one row per gene
                idx = ~pd.Index(bl.loc[:, 'Gene Symbol']).duplicated()
                bl = bl.loc[idx].set_index('Gene Symbol')
                bl = bl.loc[:, incl_cols]
                blocks['GBM%s_%s' % (pid, k)] = bl
        return blocks

    # a) All probe classes, providing they are concordant
    ipa_subdir = os.path.join(ipa_outdir, "all_concordant")
    if not os.path.exists(ipa_subdir):
        os.makedirs(ipa_subdir)
    blocks = prepare_de_dmr_concordant_blocks(dictionary.dict_by_sublevel(joint_de_dmr_sets, 4, 'concordant'), key='all')
    results_to_ipa_format(blocks, ipa_subdir, identifier=None)

    # # b) TSS, providing they are concordant
    # ipa_subdir = os.path.join(ipa_outdir, "tss_concordant")
    # if not os.path.exists(ipa_subdir):
    #     os.makedirs(ipa_subdir)
    # blocks = prepare_de_dmr_concordant_blocks(joint_de_dmr, key='tss')
    # results_to_ipa_format(blocks, ipa_subdir, identifier=None)
    #
    # # c) CpG island, providing they are concordant
    # ipa_subdir = os.path.join(ipa_outdir, "island_concordant")
    # if not os.path.exists(ipa_subdir):
    #     os.makedirs(ipa_subdir)
    # blocks = prepare_de_dmr_concordant_blocks(joint_de_dmr, key='island')
    # results_to_ipa_format(blocks, ipa_subdir, identifier=None)


    # integrate the two results
    dmr_matched = dictionary.dict_by_sublevel(test_results_significant, 2, 'matched')
    dmr_matched_only = dictionary.dict_by_sublevel(test_results_exclusive, 2, 'matched')
    dmr_ref = dictionary.dict_by_sublevel(test_results_significant, 2, 'gibco')
    dmr_ref_only = dictionary.dict_by_sublevel(test_results_exclusive, 2, 'gibco')
    dmr_intersection = dictionary.dict_by_sublevel(test_results_inclusive, 2, 'matched')

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

    # write DE / DMR results to disk
    # 1. Excel (human-readable)
    blocks = {}

    for k, v in joint_de_dmr.items():
        this_tbl = joint_de_dmr_table(v)
        for pid in this_tbl:
            blocks['GBM%s_%s' % (pid, k)] = this_tbl[pid]

    outfile = os.path.join(outdir, 'individual_gene_lists_de_dmr.xlsx')
    results_to_excel(blocks, outfile, write_index=False)

    # 2. IPA format
    incl_cols = ['logFC', 'FDR']

    ipa_outdir = os.path.join(outdir, 'ipa_dmr')
    if not os.path.exists(ipa_outdir):
        os.makedirs(ipa_outdir)

    # We want several different versions
    def prepare_de_dmr_concordant_blocks(dat, key='all'):
        blocks = {}
        for k, v in dat.items():
            for pid in v:
                bl = v[pid][key]
                # reduce to concordant results
                dmr_sign = np.sign(bl.loc[:, 'me_mediandelta'].astype(float))
                de_sign = np.sign(bl.loc[:, 'logFC'])
                bl = bl.loc[dmr_sign != de_sign]
                # only keep one row per gene
                idx = ~pd.Index(bl.loc[:, 'Gene Symbol']).duplicated()
                bl = bl.loc[idx].set_index('Gene Symbol')
                bl = bl.loc[:, incl_cols]
                blocks['GBM%s_%s' % (pid, k)] = bl
        return blocks

    # a) All probe classes, providing they are concordant
    ipa_subdir = os.path.join(ipa_outdir, "all_concordant")
    if not os.path.exists(ipa_subdir):
        os.makedirs(ipa_subdir)
    blocks = prepare_de_dmr_concordant_blocks(joint_de_dmr, key='all')
    results_to_ipa_format(blocks, ipa_subdir, identifier=None)

    # b) TSS, providing they are concordant
    ipa_subdir = os.path.join(ipa_outdir, "tss_concordant")
    if not os.path.exists(ipa_subdir):
        os.makedirs(ipa_subdir)
    blocks = prepare_de_dmr_concordant_blocks(joint_de_dmr, key='tss')
    results_to_ipa_format(blocks, ipa_subdir, identifier=None)

    # c) CpG island, providing they are concordant
    ipa_subdir = os.path.join(ipa_outdir, "island_concordant")
    if not os.path.exists(ipa_subdir):
        os.makedirs(ipa_subdir)
    blocks = prepare_de_dmr_concordant_blocks(joint_de_dmr, key='island')
    results_to_ipa_format(blocks, ipa_subdir, identifier=None)


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
    #         tmp = dictionary.dict_by_sublevel(meth_de, 2, cls)
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
