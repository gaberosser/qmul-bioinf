from methylation import loader, dmr, process
from plotting import clustering, common
from stats import transformations
import pandas as pd
import numpy as np
import copy
import os
from utils import output
import references
from scipy.stats import zscore
import multiprocessing as mp


def load_methylation(pids, ref_names=None, norm_method='swan', ref_name_filter=None):
    """
    Load and prepare the Illumina methylation data
    """
    # patient data
    obj = loader.load_by_patient(pids, norm_method=norm_method)
    anno = loader.load_illumina_methylationepic_annotation()

    # reference data
    if ref_names is not None:
        ref_obj = loader.load_reference(ref_names, norm_method=norm_method)
        if ref_name_filter is not None:
            ref_obj.filter_by_sample_name(ref_name_filter, exact=True)
        obj = loader.loader.MultipleBatchLoader([obj, ref_obj])

    me_data = obj.data.dropna()
    me_data = process.m_from_beta(me_data)

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(me_data.index)

    anno = anno.loc[common_probes]
    # dmr.add_merged_probe_classes(anno)
    me_data = me_data.loc[common_probes]
    obj.data = me_data

    return obj, anno


def pair_dmr(me_meta, me_data, dmr_clusters, pids, type1='iPSC', type2='FB', **dmr_params):
    res = {}

    for pid in pids:
        this = dmr_clusters.copy()
        the_idx1 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == type1)
        the_idx2 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == type2)
        the_idx = the_idx1 | the_idx2
        the_groups = me_meta.loc[the_idx, 'type'].values
        the_samples = me_meta.index[the_idx].groupby(the_groups)
        the_samples = [the_samples[type1], the_samples[type2]]

        this.test_clusters(
            me_data,
            samples=the_samples,
            n_jobs=dmr_params['n_jobs'],
            min_median_change=dmr_params['delta_m_min'],
            method=dmr_params['dmr_test_method'],
            alpha=dmr_params['alpha'],
            **dmr_params['test_kwargs']
        )
        res[pid] = this
    return dmr.DmrResultCollection(**res)

if __name__ == "__main__":
    pids = ['019', '030', '031', '050', '054']

    norm_method = 'pbc'
    dmr_params = {
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'alpha': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': mp.cpu_count(),
    }

    # our data
    me_obj, anno = load_methylation(pids, norm_method=norm_method)

    # ref data
    refs = [
        ('Kim et al.', loader.gse38216(norm_method=norm_method, samples=['H9 ESC 1', 'H9 ESC 2', 'H9 NPC 1', 'H9 NPC 2'])),
        ('Zimmerlin et al.', loader.gse65214(norm_method=norm_method)),
        ('Encode EPIC', loader.encode_epic(norm_method=norm_method)),
        ('Encode 450k', loader.encode_450k(norm_method=norm_method)),
    ]
    for bid, r in refs:
        r.batch_id = bid
    ref_obj = loader.loader.MultipleBatchLoader([t[1] for t in refs])

    # drop unneeded ref data
    ref_obj.meta = ref_obj.meta.loc[ref_obj.meta.index.str.contains('ESC')]

    clusters = []
    cid = 0
    for cc in anno.CHR.unique():
        coords = anno.loc[anno.CHR == cc, 'MAPINFO'].sort_values()
        this_clust = dmr.identify_cluster(coords, dmr_params['n_min'], dmr_params['d_max'])
        for cl in this_clust.values():
            clusters.append(
                dmr.ProbeCluster(cl, anno, cluster_id=cid, chr=cc)
            )
            cid += 1
    dmr_clusters = dmr.DmrResults(clusters=clusters, anno=anno)

    dmr_res = pair_dmr(me_obj.meta, me_obj.data, dmr_clusters, pids, **dmr_params)

