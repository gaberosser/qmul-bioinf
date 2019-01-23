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


def load_methylation_data(samples, anno=None, norm_method='swan'):
    me_obj = loader.load_by_sample_name(
        samples,
        norm_method=norm_method
    )

    if anno is None:
        return me_obj
    else:
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
        logger.info("%s. Computing results.", k)
        fn = get_hashed_filename(samples, norm_method=norm_method, dmr_params=dmr_params)
        me_obj, this_anno = load_methylation_data(samples, anno, norm_method=norm_method)
        me_data = process.m_from_beta(me_obj.data)

        data_loaded[k] = me_obj
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
        logger.info("%s. Computing results.", k)
        fn = get_hashed_filename(samples, norm_method=norm_method, dmr_params=dmr_params)
        me_obj, this_anno = load_methylation_data(samples, anno, norm_method=norm_method)
        me_data = process.m_from_beta(me_obj.data)

        data_loaded[k] = me_obj
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
        logger.info("%s. Computing results.", k)
        fn = get_hashed_filename(samples, norm_method=norm_method, dmr_params=dmr_params)
        me_obj, this_anno = load_methylation_data(samples, anno, norm_method=norm_method)
        me_data = process.m_from_beta(me_obj.data)

        data_loaded[k] = me_obj
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

    # limit to the PIDs that are present in all comparisons
    # compute specific DMRs and plot direction
    common_pids = sorted(setops.reduce_intersection(*pids_included.values()))
    dmrs_specific = {}
    for k in pids_included:
        dmrs_specific[k] = {}
        this = [all_results[k].results_significant[p] for p in common_pids]
        this_specific = setops.specific_features(*this)
        for p, cids, spec_dict in zip(common_pids, this_specific, this):
            if p not in dmrs_specific[k]:
                dmrs_specific[k][p] = {}
            for cid in cids:
                dmrs_specific[k][p][cid] = spec_dict[cid]

    for k in pids_included:
        fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(5.5, 4.6))

        addd.direction_of_dm_bar_plot(
            all_results[k].results_significant,
            pids=common_pids,
            as_pct=True,
            ax=axs[0, 0],
            legend=False
        )
        axs[0, 0].set_ylabel('% DMRs')
        addd.direction_of_dm_bar_plot(
            all_results[k].results_significant,
            pids=common_pids,
            as_pct=False,
            ax=axs[1, 0],
            legend=False
        )
        axs[1, 0].set_ylabel('Number DMRs')
        axs[0, 0].set_title('Full DMR list')

        addd.direction_of_dm_bar_plot(
            dmrs_specific[k],
            pids=common_pids,
            as_pct=True,
            ax=axs[0, 1],
            legend=False
        )
        axs[0, 1].set_title('Specific DMR list')
        addd.direction_of_dm_bar_plot(
            dmrs_specific[k],
            pids=common_pids,
            as_pct=False,
            ax=axs[1, 1],
            legend=False
        )
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, '%s_direction.png' % k.replace(' ', '_')), dpi=200)

    # Same KS plot we've used on the full syngeneic dataset
    from scripts.methylation.analyse_dmr_direction_and_distribution import ks_test_dmr_locations

    # we'll need chromosome lengths
    chroms = [str(t) for t in range(1, 23)]
    fa_fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human/ensembl/GRCh38.release87/fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    )
    cl = genomics.feature_lengths_from_fasta(fa_fn, features=chroms)

    ks_res = collections.defaultdict(dict)
    ks_res_specific = collections.defaultdict(dict)

    for k in pids_included:
        # BG: same for all patients, so just pick the first
        tmp = addd.get_dmr_locations(
            {pids_included[k][0]: all_results[k].results[pids_included[k][0]]},
            all_results[k].clusters,
            chrom_lengths=cl,
            split_by_direction=False
        )
        this_bg_loci = tmp['dmr_loci'][pids_included[k][0]]

        # full DMR list
        this_res = all_results[k].results_significant
        tmp = addd.get_dmr_locations(
            all_results[k].results_significant,
            all_results[k].clusters,
            chrom_lengths=cl,
            split_by_direction=True
        )
        this_hyper = tmp['dmr_loci_hyper']
        this_hypo = tmp['dmr_loci_hypo']

        hyper_fdr, hyper_stat = ks_test_dmr_locations(
            tmp['dmr_loci_hyper'],
            this_bg_loci
        )
        ks_res[k]['hyper'] = pd.DataFrame(hyper_fdr).loc[chroms, pids_included[k]]

        hypo_fdr, hypo_stat = ks_test_dmr_locations(
            tmp['dmr_loci_hypo'],
            this_bg_loci
        )
        ks_res[k]['hypo'] = pd.DataFrame(hypo_fdr).loc[chroms, pids_included[k]]

        # specific DMR list
        tmp = addd.get_dmr_locations(
            dmrs_specific[k],
            all_results[k].clusters,
            chrom_lengths=cl,
            split_by_direction=True
        )
        this_hyper_specific = tmp['dmr_loci_hyper']
        this_hypo_specific = tmp['dmr_loci_hypo']

        hyper_specific_fdr, _ = ks_test_dmr_locations(
            tmp['dmr_loci_hyper'],
            this_bg_loci
        )
        ks_res_specific[k]['hyper'] = pd.DataFrame(hyper_fdr).loc[chroms, common_pids]

        hypo_fdr, hypo_stat = ks_test_dmr_locations(
            tmp['dmr_loci_hypo'],
            this_bg_loci
        )
        ks_res_specific[k]['hypo'] = pd.DataFrame(hypo_fdr).loc[chroms, common_pids]

    # plot width scaling
    width_slope = 0.75
    width_offset = 2.
    for k in pids_included:
        the_key = k.lower().replace(' ', '_')
        vmin = .5
        vmax = 3.
        alpha = 0.1
        gs = plt.GridSpec(1, 3, width_ratios=[9, 9, 1])
        # fig = plt.figure(figsize=(width_offset + width_slope * len(pids_included[k]), 5.5))
        fig = plt.figure(figsize=(width_offset + width_slope * len(common_pids), 5.5))
        ax_hyper = fig.add_subplot(gs[0])
        ax_hypo = fig.add_subplot(gs[1], sharex=ax_hyper, sharey=ax_hyper)
        ax_cbar = fig.add_subplot(gs[2])

        # hyper_to_plot = -np.log10(ks_res[k]['hyper'])
        hyper_to_plot = -np.log10(ks_res[k]['hyper'].loc[:, common_pids])
        hyper_to_plot = hyper_to_plot[hyper_to_plot > -np.log10(alpha)]
        # hypo_to_plot = -np.log10(ks_res[k]['hypo'])
        hypo_to_plot = -np.log10(ks_res[k]['hypo'].loc[:, common_pids])
        hypo_to_plot = hypo_to_plot[hypo_to_plot > -np.log10(alpha)]

        sns.heatmap(
            hyper_to_plot,
            cmap='Reds',
            vmin=vmin,
            vmax=vmax,
            ax=ax_hyper,
            cbar=False,
            linewidths=1.,
            linecolor='w'
        )
        sns.heatmap(
            hypo_to_plot,
            cmap='Reds',
            vmin=vmin,
            vmax=vmax,
            ax=ax_hypo,
            linewidths=1.,
            linecolor='w',
            cbar_ax=ax_cbar
        )
        plt.setp(ax_hypo.yaxis.get_ticklabels(), visible=False)
        ax_hyper.set_ylabel('Chromosome')
        ax_hyper.set_xlabel('Patient')
        ax_hyper.set_title('Hypermethylated DMRs')
        ax_hypo.set_xlabel('Patient')
        ax_hypo.set_title('Hypomethylated DMRs')
        ax_cbar.set_ylabel('$-\log_{10}(p)$')
        plt.setp(ax_hyper.yaxis.get_ticklabels(), rotation=0)

        gs.update(left=0.1, right=0.88, top=0.95, wspace=0.1)
        fig.savefig(os.path.join(outdir, "ks_fdr_all_%s.png" % the_key), dpi=200)
        fig.savefig(os.path.join(outdir, "ks_fdr_all_%s.tiff" % the_key), dpi=200)

        fig = plt.figure(figsize=(width_offset + width_slope * len(common_pids), 5.5))
        ax_hyper = fig.add_subplot(gs[0])
        ax_hypo = fig.add_subplot(gs[1], sharex=ax_hyper, sharey=ax_hyper)
        ax_cbar = fig.add_subplot(gs[2])

        hyper_to_plot = -np.log10(ks_res_specific[k]['hyper'])
        hyper_to_plot = hyper_to_plot[hyper_to_plot > -np.log10(alpha)]
        hypo_to_plot = -np.log10(ks_res_specific[k]['hypo'])
        hypo_to_plot = hypo_to_plot[hypo_to_plot > -np.log10(alpha)]

        sns.heatmap(
            hyper_to_plot,
            cmap='Reds',
            vmin=vmin,
            vmax=vmax,
            ax=ax_hyper,
            cbar=False,
            linewidths=1.,
            linecolor='w'
        )
        sns.heatmap(
            hypo_to_plot,
            cmap='Reds',
            vmin=vmin,
            vmax=vmax,
            ax=ax_hypo,
            linewidths=1.,
            linecolor='w',
            cbar_ax=ax_cbar
        )
        plt.setp(ax_hypo.yaxis.get_ticklabels(), visible=False)
        ax_hyper.set_ylabel('Chromosome')
        ax_hyper.set_xlabel('Patient')
        ax_hyper.set_title('Hypermethylated DMRs')
        ax_hypo.set_xlabel('Patient')
        ax_hypo.set_title('Hypomethylated DMRs')
        ax_cbar.set_ylabel('$-\log_{10}(p)$')
        plt.setp(ax_hyper.yaxis.get_ticklabels(), rotation=0)

        fig.savefig(os.path.join(outdir, "ks_fdr_specific_%s.png" % the_key), dpi=200)
        fig.savefig(os.path.join(outdir, "ks_fdr_specific_%s.tiff" % the_key), dpi=200)
