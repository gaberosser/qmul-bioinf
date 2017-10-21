from load_data import rnaseq_data, methylation_array
from rnaseq.differential_expression import edger
from rnaseq.filter import filter_by_cpm
from methylation import process, dmr
from integrator.rnaseq_methylationarray import compute_joint_de_dmr
from integrator import plots
import pandas as pd
import numpy as np
from scipy import stats
import references
import os
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


def compute_de(rnaseq_obj, pids):
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

            de_matched[pid] = edger(the_data, the_groups, the_contrast)
            add_gene_symbols(de_matched[pid])

            # repeat with gibco reference
            # use the same genes, rather than filtering again
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.type == 'GBM')) | rnaseq_obj.meta.index.str.contains('GIBCO')
            the_data = rnaseq_obj.data.loc[the_genes, the_idx]
            the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
            the_contrast = "GBM - NSC"

            de_gibco[pid] = edger(the_data, the_groups, the_contrast)
            add_gene_symbols(de_gibco[pid])

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


if __name__ == "__main__":
    # if this is specified, we load the DMR results from a JSON rather than recomputing them to save time
    DMR_LOAD_DIR = None

    outdir = output.unique_output_dir("paired_rnaseq")
    ref_name = 'GIBCONSC_P4'
    pids = ['017', '050', '054', '061']

    dmr_params = {
        'core_min_sample_overlap': 3,  # 3 / 4 samples must match
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'fdr': 0.05,
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

    de_res = compute_de(rnaseq_obj, pids)
    de = de_res['de']
    de_matched = de_res['de_matched']
    de_gibco = de_res['de_gibco']

    # split results by Venn group
    de_matched_only = dict([
        (pid, de_matched[pid].loc[de[pid]['10']]) for pid in pids
    ])
    de_gibco_only = dict([
        (pid, de_gibco[pid].loc[de[pid]['01']]) for pid in pids
    ])
    # DE intersection: use the matched comparison values
    de_intersection = dict([
        (pid, de_matched[pid].loc[de[pid]['11']]) for pid in pids
    ])


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
    ## TODO: load from JSON if specified
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

    # side story: check overlap between individuals and group
    # skip this for now but useful to know!

    # integrate the two results
    dmr_matched = dmr.dict_by_sublevel(test_results_significant, 2, 'matched')
    dmr_matched_only = dmr.dict_by_sublevel(test_results_exclusive, 2, 'matched')
    dmr_gibco = dmr.dict_by_sublevel(test_results_significant, 2, 'gibco')
    dmr_gibco_only = dmr.dict_by_sublevel(test_results_exclusive, 2, 'gibco')
    dmr_intersection = dmr.dict_by_sublevel(test_results_inclusive, 2, 'matched')

    ## FIXME: failed to add data error

    joint_de_dmr = {
        ('matched', 'matched'): compute_joint_de_dmr(dmr_matched, de_matched),
        ('matched_only', 'matched_only'): compute_joint_de_dmr(dmr_matched_only, de_matched_only),
        ('insc_ref', 'insc_ref'): compute_joint_de_dmr(dmr_intersection, de_intersection),
    }
