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
from scripts.hgic_final import analyse_dmrs_s1_direction_distribution as addd
from plotting import common

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
logger = log.get_console_logger()


def get_hashed_filename(
    samples,
    norm_method='swan',
    dmr_params=consts.DMR_PARAMS,
    load_dir=os.path.join(output.OUTPUT_DIR, 'dmr')
):
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method

    the_hash = tsgd.dmr_results_hash(samples, dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash

    if load_dir is not None:
        return os.path.join(load_dir, filename)
    else:
        return filename


def load_dmr_results(anno, samples, norm_method='swan', dmr_params=consts.DMR_PARAMS):
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method

    the_hash = tsgd.dmr_results_hash(samples, dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        return dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise IOError("Unable to find pre-computed results file %s." % fn)


def load_methylation_data(samples, anno, norm_method='swan'):
    me_obj = loader.load_by_sample_name(
        samples,
        norm_method=norm_method
    )
    common_probes = anno.index.intersection(me_obj.data.index)

    this_anno = anno.loc[common_probes]
    me_obj.data = me_obj.data.loc[common_probes]

    return me_obj, this_anno


if __name__ == "__main__":
    """
    Here we carry out similar analysis to that in analyse_dmr_direction_and_distribution, but we focus on using 
    different comparators. 
    Where syngeneic GIC-iNSC comparisons were used in the original script, here we use different syngeneic cell types. 
    We also run (non-syngeneic) cross-comparisons.
    """

    pids = consts.PIDS
    norm_method = 'swan'
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    # set this to True if output bed files are required (this is quite slow due to the large number of combinations)
    write_bed_files = False

    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    # we'll need the same array annotation data for all our results
    anno = loader.load_illumina_methylationepic_annotation()

    # this will help us retain all our loaded objects
    data_loaded = {}
    pids_included = {'GIC-iNSC syn': pids}
    all_results = {}

    # load previous syngeneic results
    dmr_res_insc_syn = load_dmr_results(anno, consts.S1_METHYL_SAMPLES_GIC + consts.S1_METHYL_SAMPLES_INSC)
    all_results['GIC-iNSC syn'] = dmr_res_insc_syn

    # GIC-FB (syngeneic)
    k = 'GIC-FB syn'
    samples = consts.S1_METHYL_SAMPLES_GIC + consts.S1_METHYL_SAMPLES_FB
    this_pids = ['018', '019', '030', '031', '017', '050', '054', '026', '052']
    pids_included[k] = this_pids

    try:
        dmr_res_fb_syn = load_dmr_results(anno, samples)
    except IOError:
        logger.info("GIC - FB (syngeneic). Computing results.")
        fn = get_hashed_filename(samples, norm_method=norm_method, dmr_params=dmr_params)
        me_obj, this_anno = load_methylation_data(samples, anno, norm_method=norm_method)
        me_data = process.m_from_beta(me_obj.data)

        data_loaded['GIC-FB syn'] = me_obj
        dmr_res_fb_syn = tsgd.paired_dmr(
            me_data,
            me_obj.meta,
            this_anno,
            this_pids,
            dmr_params,
            type1='GBM',
            type2='FB'
        )
        dmr_res_fb_syn.to_pickle(fn, include_annotation=False)
        logger.info("Saved DMR results to %s", fn)

    all_results[k] = dmr_res_fb_syn

    # GIC-iAPC (syngeneic)
    k = 'GIC-iAPC syn'
    samples = consts.S1_METHYL_SAMPLES_GIC + consts.S1_METHYL_SAMPLES_IAPC
    this_pids = ['019', '031', '050', '052']
    pids_included[k] = this_pids

    try:
        dmr_res_iapc_syn = load_dmr_results(anno, samples)
    except IOError:
        logger.info("GIC - iOPC (syngeneic). Computing results.")
        fn = get_hashed_filename(samples, norm_method=norm_method, dmr_params=dmr_params)
        me_obj, this_anno = load_methylation_data(samples, anno, norm_method=norm_method)
        me_data = process.m_from_beta(me_obj.data)

        data_loaded['GIC-iAPC syn'] = me_obj
        dmr_res_iapc_syn = tsgd.paired_dmr(
            me_data,
            me_obj.meta,
            this_anno,
            this_pids,
            dmr_params,
            type1='GBM',
            type2='iAPC'
        )
        dmr_res_iapc_syn.to_pickle(fn, include_annotation=False)
        logger.info("Saved DMR results to %s", fn)

    all_results[k] = dmr_res_iapc_syn

    # GIC-iOPC (syngeneic)
    k = 'GIC-iOPC syn'
    samples = consts.S1_METHYL_SAMPLES_GIC + consts.S1_METHYL_SAMPLES_IOPC
    this_pids = ['019', '031', '050', '052']
    pids_included[k] = this_pids

    try:
        dmr_res_iopc_syn = load_dmr_results(anno, samples)
    except IOError:
        logger.info("GIC - iOPC (syngeneic). Computing results.")
        fn = get_hashed_filename(samples, norm_method=norm_method, dmr_params=dmr_params)
        me_obj, this_anno = load_methylation_data(samples, anno, norm_method=norm_method)
        me_data = process.m_from_beta(me_obj.data)

        data_loaded['GIC-iOPC syn'] = me_obj
        dmr_res_iopc_syn = tsgd.paired_dmr(
            me_data,
            me_obj.meta,
            this_anno,
            this_pids,
            dmr_params,
            type1='GBM',
            type2='iOPC'
        )
        dmr_res_iopc_syn.to_pickle(fn, include_annotation=False)
        logger.info("Saved DMR results to %s", fn)

    all_results[k] = dmr_res_iopc_syn

    # bar plot showing balance of DMR direction
    for k in pids_included:
        fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(5.5, 5.5))

        addd.direction_of_dm_bar_plot(
            all_results[k].results_significant,
            pids=pids_included[k],
            as_pct=True,
            ax=axs[0],
            legend=False
        )
        axs[0].set_ylabel('% DMRs')
        addd.direction_of_dm_bar_plot(
            all_results[k].results_significant,
            pids=pids_included[k],
            as_pct=False,
            ax=axs[1],
            legend=False
        )
        axs[0].set_title(k)
        axs[1].set_ylabel('Number DMRs')
        fig.savefig(os.path.join(outdir, '%s_direction.png' % k.replace(' ', '_')), dpi=200)