from plotting import bar, common, pie
from methylation import loader, dmr, process
import pandas as pd
from statsmodels.sandbox.stats import multicomp
from utils import output, setops, genomics, log
import multiprocessing as mp
import os
import collections
import pickle
import numpy as np
from scipy import stats, cluster
import matplotlib
from matplotlib import pyplot as plt, patches
from matplotlib.colors import Normalize
from matplotlib import cm
from sklearn.neighbors import KernelDensity

import seaborn as sns
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
logger = log.get_console_logger()

"""
Here we analyse the direction and genomic locations of DMRs in a 450K validation cohort (GSE92462).
We re-analyse our own (EPIC) samples, including only the 450K probes.
We test whether we see similar patterns in the validation cohort.
"""

def run_dmr_analyses(data, comparisons, anno, dmr_params, verbose=True):
    """
    Compute DMRs for paired GBM-iNSC comparisons (defined in that order) for all patients
    :param me_data: Pandas dataframe containing M values, columns are samples and rows are probes.
    :param comparisons: Dictionary, each key is a title, each value is a 2-element iterable containing lists of
    sample names. The comparison is run as group 1 - group 2 in each case.
    :param anno:
    :param pids:
    :param dmr_params:
    :return:
    """
    dmr_res_obj = dmr.DmrResults(anno=anno)
    dmr_res_obj.identify_clusters(**dmr_params)
    dmr_res = {}

    for the_ttl, the_samples in comparisons.items():
        logger.info(
            "Comparison %s. Group 1: %s. Group 2: %s.",
            the_ttl,
            ','.join(the_samples[0]),
            ','.join(the_samples[1]),
        )
        the_obj = dmr_res_obj.copy()
        the_obj.test_clusters(data,
                              samples=the_samples,
                              n_jobs=dmr_params['n_jobs'],
                              min_median_change=dmr_params['delta_m_min'],
                              method=dmr_params['dmr_test_method'],
                              alpha=dmr_params['alpha'],
                              **dmr_params['test_kwargs']
                              )
        dmr_res[the_ttl] = the_obj

    return dmr.DmrResultCollection(**dmr_res)


def bar_plot(res, keys=None):
    from scripts.hgic_final import analyse_dmrs_s1_direction_distribution as addd

    if keys is None:
        keys = sorted(res.keys())

    # bar plot showing balance of DMR direction
    fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(5.5, 5.5))

    addd.direction_of_dm_bar_plot(
        res,
        pids=keys,
        as_pct=True,
        ax=axs[0],
        legend=False
    )
    axs[0].set_ylabel('% DMRs')
    addd.direction_of_dm_bar_plot(
        res,
        pids=keys,
        as_pct=False,
        ax=axs[1],
        legend=False
    )
    axs[1].set_ylabel('Number DMRs')
    plt.setp(axs[1].xaxis.get_ticklabels(), rotation=90)

    return {
        'axs': axs,
        'fig': fig
    }


if __name__ == "__main__":
    pids = consts.PIDS
    norm_method = 'swan'
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    # load our data
    our_obj = loader.load_by_patient(pids, norm_method=norm_method, samples=consts.S1_METHYL_SAMPLES)
    anno = loader.load_illumina_methylationepic_annotation()
    our_obj.meta.insert(0, 'patient_id', our_obj.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    # load validation data
    val_obj = loader.load_reference('GSE92462_450k', norm_method=norm_method)
    # filter
    val_obj.filter_samples(val_obj.meta.type.isin(['GBM (GSC)', 'NSC']))

    # TODO: upload to the classifier and run (toggle this so it's only run once)

    # combine and reduce probes
    obj = loader.loader.MultipleBatchLoader([our_obj, val_obj])
    dat = process.m_from_beta(obj.data)
    meta = obj.meta
    common_probes = anno.index.intersection(dat.index)
    dat = dat.reindex(common_probes)
    anno = anno.reindex(common_probes)

    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method

    the_hash = tsgd.dmr_results_hash(meta.sort_index().index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_450k_validation.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        dmr_res = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        logger.info("Unable to locate pre-existing results. Computing from scratch (this can take a while).")
        comparisons = {}
        # consistency (1): our GIC vs our iNSC (to check that the reduced probes still show the phenotypes)
        for pid in pids:
            our_inscs = meta.index[(meta.type == 'iNSC') & (meta.patient_id == pid)]
            our_gics = meta.index[(meta.type == 'GBM') & (meta.patient_id == pid)]
            comparisons['syngeneic_%s' % pid] = [our_gics, our_inscs]

        # consistency (2): our GIC vs validation references
        # (to check that we still see the phenotype with a different comparator)
        val_refs = meta.index[meta.type == 'NSC']
        for pid in pids:
            our_gics = meta.index[(meta.type == 'GBM') & (meta.patient_id == pid)]
            comparisons['consistency_%s' % pid] = [our_gics, val_refs]

        # validation: validation GIC vs validation references
        val_gics = meta.index[meta.type == 'GBM (GSC)']
        for g in val_gics:
            comparisons['validation_%s' % g] = [pd.Index([g]), val_refs]

        dmr_res = run_dmr_analyses(dat, comparisons, anno, dmr_params)
        # Save DMR results to disk
        dmr_res.to_pickle(fn, include_annotation=False)
        logger.info("Saved DMR results to %s", fn)

    dmr_res_all = dmr_res.results_significant

    # 1. check the phenomenon is still observed in our syngeneic comparisons
    # full list
    for_plot = dict([
        (pid, dmr_res_all['syngeneic_%s' % pid]) for pid in pids
    ])
    plt_dict = bar_plot(for_plot, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "syngeneic_full_list_directions.png"), dpi=200)

    # specific list
    spec_ix = setops.specific_features(*[for_plot[pid].keys() for pid in pids])
    for_plot = dict([
        (pid, dict([(k, for_plot[pid][k]) for k in s])) for pid, s in zip(pids, spec_ix)
    ])
    plt_dict = bar_plot(for_plot, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "syngeneic_specific_list_directions.png"), dpi=200)

    # 2. check that we see the same phenomenon when we switch to the validation cohort comparator
    # full list
    for_plot = dict([
        (pid, dmr_res_all['consistency_%s' % pid]) for pid in pids
    ])
    plt_dict = bar_plot(for_plot, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "consistency_full_list_directions.png"), dpi=200)

    # specific list
    spec_ix = setops.specific_features(*[for_plot[pid].keys() for pid in pids])
    for_plot = dict([
        (pid, dict([(k, for_plot[pid][k]) for k in s])) for pid, s in zip(pids, spec_ix)
    ])
    plt_dict = bar_plot(for_plot, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "consistency_specific_list_directions.png"), dpi=200)

    # 3. how does the same plot look in the validation cohort?
    # full list
    for_plot = dict([
        (k.replace('validation_', ''), dmr_res_all[k]) for k in dmr_res_all if 'validation' in k
    ])
    plt_dict = bar_plot(for_plot)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "validation_full_list_directions.png"), dpi=200)

    # specific list
    keys = sorted(for_plot.keys())
    spec_ix = setops.specific_features(*[for_plot[k].keys() for k in keys])
    for_plot = dict([
        (j, dict([(k, for_plot[j][k]) for k in s])) for j, s in zip(keys, spec_ix)
    ])
    plt_dict = bar_plot(for_plot)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "validation_specific_list_directions.png"), dpi=200)




