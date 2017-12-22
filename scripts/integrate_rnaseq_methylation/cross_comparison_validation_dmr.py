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


if __name__ == "__main__":
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'integrate_rnaseq_methylation')
    outdir = output.unique_output_dir("cross_validate_dmr", reuse_empty=True)
    ref_name = 'GIBCONSC_P4'
    # all n=2 samples and RTK II samples
    pids = ['017', '019', '030', '031', '050', '054']
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
    po_dmr = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    ro_dmr = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    pr_concordant_dmr = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    pr_discordant_dmr = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    for pid in pids:
        for pid2 in pids + ['GIBCO']:
            p = dmr_sign.loc[pid, pid]
            r = dmr_sign.loc[pid, pid2]
            pres = dmr_res[pid][pid].results_significant
            rres = dmr_res[pid][pid2].results_significant
            x, _ = setops.venn_from_arrays(p, r)
            po_dmr.loc[pid, pid2] = x['10']
            ro_dmr.loc[pid, pid2] = x['01']
            # ref and pair IDs
            pr_id = x['11']
            # signs
            pmed_change_sign = np.array([np.sign(pres[t]['median_change']) for t in pr_id])
            rmed_change_sign = np.array([np.sign(rres[t]['median_change']) for t in pr_id])

            pr_concordant_dmr.loc[pid, pid2] = list(
                np.array(x['11'])[pmed_change_sign == rmed_change_sign]
            )

            pr_discordant_dmr.loc[pid, pid2] = list(
                np.array(x['11'])[pmed_change_sign != rmed_change_sign]
            )
