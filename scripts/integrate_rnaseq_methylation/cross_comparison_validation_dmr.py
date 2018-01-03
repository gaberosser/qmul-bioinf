from methylation import dmr, process
import os
import re
import collections
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import references
from settings import LOCAL_DATA_DIR
from utils import output, setops
from load_data import methylation_array


def compute_cross_dmr(me_data, me_meta, anno, pids, dmr_params, external_references=(('GIBCO', 'NSC'),)):

    obj = dmr.DmrResults(anno=anno)
    obj.identify_clusters(**dmr_params)
    res = {}

    # loop over GBM groups
    for pid1 in pids:
        res.setdefault(pid1, {})
        the_idx1 = me_meta.index.str.contains(pid1) & (me_meta.loc[:, 'type'] == 'GBM')
        # loop over iNSC groups
        for pid2 in pids:
            the_idx2 = me_meta.index.str.contains(pid2) & (me_meta.loc[:, 'type'] == 'iNSC')
            the_idx = the_idx1 | the_idx2
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
            res[pid1][pid2] = the_obj

        # loop over external reference NSC groups
        for er, er_type in external_references:
            the_idx2 = me_meta.index.str.contains(er) & (me_meta.loc[:, 'type'] == er_type)
            the_idx = the_idx1 | the_idx2
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
            res[pid1][er] = the_obj

    return dmr.DmrResultCollection(**res)


def plot_methylation_heatmap(
    data,
        cluster_ids,
        dmr_res,
        ref_key='GIBCONSC_P4',
        cmap='RdYlGn_r',
        yticklabels=False,
        plot_cluster_dividers=True,
):
    probe_ids = reduce(lambda x, y: x + y, (dmr_res.clusters[t].pids for t in cluster_ids), [])
    # add manual breaks
    # this needs to be reversed due to the method of plotting a heatmap
    break_idx = np.cumsum([len(dmr_res.clusters[t].pids) for t in cluster_ids][::-1])[:-1]

    po_dat = data.loc[probe_ids]

    # rearrange columns
    cols = (
        sorted(po_dat.columns[po_dat.columns.str.contains('GBM')].tolist()) +
        [ref_key] +
        sorted(po_dat.columns[po_dat.columns.str.contains('DURA')].tolist())
    )
    po_dat = po_dat.loc[:, cols]
    # insert spacing columns
    idx = np.where(po_dat.columns == ref_key)[0][0]
    po_dat.insert(idx, '', np.nan)
    po_dat.insert(idx + 2, ' ', np.nan)

    fig = plt.figure(figsize=(7, 10))
    ax = fig.add_subplot(111)
    ax = sns.heatmap(po_dat, cmap=cmap, ax=ax, yticklabels=yticklabels)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    fig.tight_layout()
    if plot_cluster_dividers:
        [ax.plot([0, len(cols) + 2], [t, t], c='k', lw=1., ls='--', alpha=0.6) for t in break_idx]
    return ax

if __name__ == "__main__":
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'cross_validate_dmr')
    outdir = output.unique_output_dir("cross_validate_dmr", reuse_empty=True)
    ref_name = 'GIBCONSC_P4'
    pids = ['017', '019', '030', '031', '050', '054']
    subgroups = {
        'RTK I': ['019', '030', '031'],
        'RTK II': ['017', '050', '054'],
    }
    cmap = 'RdYlGn_r'

    dmr_params = {
        'core_min_sample_overlap': 3,  # 3 / 4 samples must match
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'fdr': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': 4,
    }

    intersecter = lambda x, y: set(x).intersection(y)
    unioner = lambda x, y: set(x).union(y)

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
    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    hash_elements = tuple(sorted(pids)) + (
        dmr_params['d_max'],
        dmr_params['n_min'],
        dmr_params['delta_m_min'],
        dmr_params['dmr_test_method'],
    ) + tuple(dmr_params['test_kwargs'].items())
    filename = 'dmr_results.%d.pkl' % hash(hash_elements)
    fn = os.path.join(DMR_LOAD_DIR, filename)

    loaded = False
    if DMR_LOAD_DIR is not None:
        if os.path.isfile(fn):
            dmr_res = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
            loaded = True

    if not loaded:
        dmr_res = compute_cross_dmr(me_data, me_meta, anno, pids, dmr_params)
        # Save DMR results to disk
        dmr_res.to_pickle(fn, include_annotation=False)
        print "Saved DMR results to %s" % fn

    # table of sig. cluster IDs
    dmr_sign = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    for pid in pids:
        for pid2 in pids + ['GIBCO']:
            dmr_sign.loc[pid, pid2] = sorted(dmr_res[pid][pid2].clusters_significant.keys())

    dmr_counts = dmr_sign.applymap(len)

    # pair only
    pair_only = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    ref_only = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    pair_and_ref_concordant = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    pair_and_ref_discordant = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    for pid in pids:
        for pid2 in pids + ['GIBCO']:
            p = dmr_sign.loc[pid, pid]
            r = dmr_sign.loc[pid, pid2]
            pres = dmr_res[pid][pid].results_significant
            rres = dmr_res[pid][pid2].results_significant
            x, _ = setops.venn_from_arrays(p, r)
            pair_only.loc[pid, pid2] = x['10']
            ref_only.loc[pid, pid2] = x['01']
            # ref and pair IDs
            pr_id = x['11']
            # signs
            pmed_change_sign = np.array([np.sign(pres[t]['median_change']) for t in pr_id])
            rmed_change_sign = np.array([np.sign(rres[t]['median_change']) for t in pr_id])

            pair_and_ref_concordant.loc[pid, pid2] = list(
                np.array(x['11'])[pmed_change_sign == rmed_change_sign]
            )

            pair_and_ref_discordant.loc[pid, pid2] = list(
                np.array(x['11'])[pmed_change_sign != rmed_change_sign]
            )

    po_counts = pair_only.applymap(len)
    ro_counts = ref_only.applymap(len)

    # identify probes that are present in every ref comparison

    po_each = [
        sorted(
            reduce(intersecter, pair_only.loc[pid, ~pair_only.columns.str.contains(pid)])
        ) for pid in pids
    ]
    po_each = pd.Series(po_each, index=pids)

    # now relax this requirement: which probes would be included if we require their inclusion in N of the cells
    # (rather than all)?
    possible_counts = range(1, pair_only.shape[1])
    po_each_threshold = pd.DataFrame(index=pids, columns=possible_counts)
    for pid in pids:
        this_counter = collections.Counter()
        # iterate over each column
        # we can include the empty diagonal cell, since it will not affect the counting
        for col in pair_only.columns:
            for e in pair_only.loc[pid, col]:
                this_counter[e] += 1
        # progressively filter the gene list based on counts
        the_genes = this_counter.keys()
        for i in possible_counts:
            the_genes = [k for k in this_counter if this_counter[k] >= i]
            po_each_threshold.loc[pid, i] = the_genes

    # ...how many of these are shared between patients?
    # consider all, K -1 and K-2
    K = len(pids)
    for i in possible_counts:
        _, cts = setops.venn_from_arrays(*po_each_threshold.loc[:, i].values)
        this_tally = []

        print "N = %d" % i
        for j in [K, K - 1, K - 2, K - 3]:
            this_ct = sum([cts[k] for k in setops.binary_combinations_sum_gte(K, j)])
            print "%d DMRs shared by >=%d patients" % (this_ct, j)
        # also look at the overlap within the subgroups
        for grp_name, grp_members in subgroups.items():
            # get the group member results
            this_po_each_threshold = po_each_threshold.loc[grp_members]
            _, cts = setops.venn_from_arrays(*this_po_each_threshold.loc[:, i].values)
            the_idx = ''.join(['1'] * len(grp_members))
            print "%d DMRs shared by all patients in subgroup %s" % (cts[the_idx], grp_name)

    # for reference: what do these numbers look like in the Gibco comparison (only)?
    po_gibco_common_counts = pd.Series(index=possible_counts)
    _, cts = setops.venn_from_arrays(*pair_only.loc[:, 'GIBCO'].values)
    for j in possible_counts:
        po_gibco_common_counts.loc[j] = sum([cts[k] for k in setops.binary_combinations_sum_gte(K, j)])
        print "%d DMRs shared by >=%d patients in the pair-only Gibco comparison" % (
            int(po_gibco_common_counts.loc[j]),
            j
        )

    # what is present in X vs Y_i that isn't in X vs any other Y?
    po_diff = pd.DataFrame(index=pair_only.index, columns=pair_only.columns)
    for pid in pids:
        for pid2 in pair_only.columns:
            the_ref = pair_only.loc[pid, pid2]
            all_else = pair_only.loc[pid, pair_only.columns != pid2]
            union_all_else = reduce(set.union, all_else, set())
            po_diff.loc[pid, pid2] = sorted(set(the_ref).difference(union_all_else))

    ro_diff = pd.DataFrame(index=ref_only.index, columns=ref_only.columns)
    for pid in pids:
        for pid2 in ref_only.columns:
            the_ref = ref_only.loc[pid, pid2]
            all_else = ref_only.loc[pid, ref_only.columns != pid2]
            union_all_else = reduce(set.union, all_else, set())
            ro_diff.loc[pid, pid2] = sorted(set(the_ref).difference(union_all_else))

    po_intersection_insc = pd.Series(index=pids)
    for pid in pids:
        # first, find the genes that are always PO when an iNSC reference is used
        tmp = reduce(intersecter, pair_only.loc[pid, pair_only.index[pair_only.index != pid]])
        # now exclude any that are also DE when the Gibco reference is used
        po_intersection_insc.loc[pid] = tmp.difference(pair_only.loc[pid, 'GIBCO'])

    po_specific_to_reference = [
        sorted(
            reduce(lambda x, y: set(x).intersection(y), po_diff.loc[~po_diff.index.str.contains(pid), pid])
        ) for pid in pids + ['GIBCO']
        ]
    po_specific_to_reference = pd.Series(po_specific_to_reference, index=pids + ['GIBCO'])

    ro_intersection_insc = pd.Series(index=pids)
    for pid in pids:
        # first, find the genes that are always PO when an iNSC reference is used
        tmp = reduce(intersecter, ref_only.loc[pid, ref_only.index[ref_only.index != pid]])
        # now exclude any that are also DE when the Gibco reference is used
        ro_intersection_insc.loc[pid] = tmp.difference(ref_only.loc[pid, 'GIBCO'])

    ro_specific_to_reference = [
        sorted(
            reduce(lambda x, y: set(x).intersection(y), ro_diff.loc[~ro_diff.index.str.contains(pid), pid])
        ) for pid in pids + ['GIBCO']
        ]
    ro_specific_to_reference = pd.Series(ro_specific_to_reference, index=pids + ['GIBCO'])


    # get the clusters that consistently differ in the pair comparison only and NOT in Gibco (across all patients)
    # these will have a methylation pattern in Gibco similar to GBM, so that they do NOT appear
    po_gibco_diff = po_specific_to_reference.loc['GIBCO']
    ro_gibco_diff = ro_specific_to_reference.loc['GIBCO']

    # can also look separately at different probe classes
    # po_gibco_diff_tss = [t for t in po_gibco_diff if 'tss' in dmr_res.clusters[t].cls]
    # po_gibco_diff_island = [t for t in po_gibco_diff if 'island' in dmr_res.clusters[t].cls]

    ax = plot_methylation_heatmap(me_data, po_gibco_diff, dmr_res)
    ax.figure.savefig(os.path.join(outdir, "consistently_in_pair_only.png"), dpi=200)

    ax = plot_methylation_heatmap(me_data, ro_gibco_diff, dmr_res, plot_cluster_dividers=False)
    ax.figure.savefig(os.path.join(outdir, "consistently_in_ref_only.png"), dpi=200)