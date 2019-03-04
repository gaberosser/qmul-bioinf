from plotting import bar, common, pie, polar
from methylation import loader, dmr, process
import pandas as pd
from stats import nht
import hgic_consts
from utils import output, setops, genomics, log, dictionary
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
from scripts.dmr_direction_bias_story import analyse_dmr_direction_and_distribution as addd

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
logger = log.get_console_logger()

"""
Here we check for the same biased differential methylation phenotype in the bulk FFPE samples, using the syngeneic
cell lines as a comparator. 
"""


if __name__ == "__main__":
    pids = consts.PIDS
    norm_method = 'swan'
    alpha = 0.05
    pk_alpha = -np.log10(alpha)

    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    # load our data
    cc_obj = loader.load_by_patient(pids, norm_method=norm_method, samples=consts.S1_METHYL_SAMPLES)
    ffpe_obj = loader.load_by_patient(pids, norm_method=norm_method, type='ffpe')

    anno = loader.load_illumina_methylationepic_annotation()
    # add patient ID column to metadata
    cc_obj.meta.insert(0, 'patient_id', cc_obj.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))
    ffpe_obj.meta.insert(0, 'patient_id', [hgic_consts.NH_ID_TO_PATIENT_ID_MAP[t] for t in ffpe_obj.meta.index])
    ffpe_obj.meta.insert(1, 'type', 'ffpe')

    dat_cc = process.m_from_beta(cc_obj.data).sort_index()
    # replace CC data with those normed differently (in R)
    # dat_cc = pd.read_csv('cell_culture_swan_one_norm.csv', header=0, index_col=0).sort_index()

    dat_ffpe = process.m_from_beta(ffpe_obj.data).sort_index()
    dat = pd.concat((dat_cc, dat_ffpe), axis=1, join='inner')
    meta = pd.concat(
        (cc_obj.meta, ffpe_obj.meta), axis=0, join='outer', sort=True
    )
    meta.loc[meta.batch.isnull(), 'batch'] = meta.loc[meta.batch.isnull(), 'batch_1']

    common_probes = anno.index.intersection(dat.index)
    dat = dat.reindex(common_probes)
    anno = anno.reindex(common_probes)

    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method

    the_hash = tsgd.dmr_results_hash(dat.sort_index(axis=1).columns.tolist(), dmr_hash_dict)
    filename = 'dmr_results_ffpe_vs_cell_culture.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        dmr_res = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        logger.info("Unable to locate pre-existing results. Computing from scratch (this can take a while).")
        comparisons = {}
        for pid in pids:  # GIC
            gic = meta.index[(meta.type == 'GBM') & (meta.patient_id == pid)]
            insc = meta.index[(meta.type == 'iNSC') & (meta.patient_id == pid)]
            gbm = meta.index[(meta.type == 'ffpe') & (meta.patient_id == pid)]
            comparisons["%s_ffpe-iNSC" % pid] = [gbm, insc]
            comparisons["%s_GIC-iNSC" % pid] = [gic, insc]
        dmr_res = addd.run_dmr_analyses(dat, comparisons, anno, dmr_params)
        # Save DMR results to disk
        dmr_res.to_pickle(fn, include_annotation=False)
        logger.info("Saved DMR results to %s", fn)

    dmr_res_all = dmr_res.results_significant

    # look at distn of M values (before any differences)
    # two figures, one for each replicate
    cols = common.get_best_cmap(len(pids))
    fig, axs = plt.subplots(ncols=2, sharex=True, sharey=True)

    for i, pid in enumerate(pids):
        c = cols[i]
        this_ix = meta.index[(meta.patient_id == pid) & (meta.type == 'GBM')]
        for j, ix in enumerate(this_ix):
            this_dat = dat.loc[:, ix]
            this_ax = axs[j]
            sns.kdeplot(this_dat, color=c, label=pid, ax=this_ax)

    [ax.set_xlim([-10, 10]) for ax in axs]
    [ax.set_xlabel("M value") for ax in axs]
    axs[0].set_ylabel("Density (a.u.)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "m_value_density_gic.png"), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for i, pid in enumerate(pids):
        c = cols[i]
        ix = meta.index[(meta.patient_id == pid) & (meta.type == 'ffpe')]
        this_dat = dat.loc[:, ix].squeeze()
        sns.kdeplot(this_dat, color=c, label=pid, ax=ax)
    ax.set_xlim([-10, 10])
    ax.set_xlabel("M value")
    ax.set_ylabel("Density (a.u.)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "m_value_density_ffpe.png"), dpi=200)

    # uh oh - batch effects? replot with colour indicating batch
    cols = dict(zip(meta.batch.unique(), common.get_best_cmap(len(meta.batch.unique()))))
    fig, axs = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(10, 5.5))
    batches_seen = collections.defaultdict(set)

    for i, pid in enumerate(pids):
        this_ix = meta.index[(meta.patient_id == pid) & (meta.type == 'GBM')]
        for j, ix in enumerate(this_ix):
            this_batch = meta.batch[ix]
            lbl = None
            if this_batch not in batches_seen[j]:
                lbl = this_batch
                batches_seen[j].add(this_batch)
            this_dat = dat.loc[:, ix].values
            this_ax = axs[j]
            sns.kdeplot(this_dat, color=cols[this_batch], ax=this_ax, label=lbl)

    [ax.set_xlim([-10, 10]) for ax in axs]
    [ax.set_xlabel("M value") for ax in axs]
    axs[0].set_ylabel("Density (a.u.)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "m_value_density_gic_by_batch.png"), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    batches_seen = set()

    for i, pid in enumerate(pids):
        ix = meta.index[(meta.patient_id == pid) & (meta.type == 'ffpe')]
        this_dat = dat.loc[:, ix].squeeze().values
        this_batch = meta.batch[ix].squeeze()
        lbl = None
        if this_batch not in batches_seen:
            lbl = this_batch
            batches_seen.add(this_batch)
        sns.kdeplot(this_dat, color=cols[this_batch], label=lbl, ax=ax)
    ax.set_xlim([-10, 10])
    ax.set_xlabel("M value")
    ax.set_ylabel("Density (a.u.)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "m_value_density_ffpe_by_batch.png"), dpi=200)


    # bar chart, syngeneic full
    for_plot = dict([
        (pid, dmr_res_all['%s_GIC-iNSC' % pid]) for pid in pids
    ])
    plt_dict = addd.bar_plot(for_plot, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "syngeneic_full_list_directions.png"), dpi=200)

    # bar chart, syngeneic specific
    spec_ix = setops.specific_features(*[for_plot[pid].keys() for pid in pids])
    for_plot = dict([
        (pid, dict([(k, for_plot[pid][k]) for k in s])) for pid, s in zip(pids, spec_ix)
    ])
    plt_dict = addd.bar_plot(for_plot, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "syngeneic_specific_list_directions.png"), dpi=200)

    # bar chart, FFPE full
    for_plot = dict([
        (pid, dmr_res_all['%s_ffpe-iNSC' % pid]) for pid in pids
    ])
    plt_dict = addd.bar_plot(for_plot, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "ffpe-iNSC_full_list_directions.png"), dpi=200)

    # bar chart, FFPE specific
    spec_ix = setops.specific_features(*[for_plot[pid].keys() for pid in pids])
    for_plot = dict([
        (pid, dict([(k, for_plot[pid][k]) for k in s])) for pid, s in zip(pids, spec_ix)
    ])
    plt_dict = addd.bar_plot(for_plot, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "ffpe-iNSC_specific_list_directions.png"), dpi=200)
