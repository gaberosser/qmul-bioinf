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
from scripts.methylation.analyse_dmr_direction_validation_cohort import run_dmr_analyses, bar_plot

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
logger = log.get_console_logger()

"""
Here we analyse the direction (and genomic locations?) of DMRs in a cross-comparison.
We compare our GIC lines with (non-matching) iNSC lines.
"""

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

    dat = process.m_from_beta(our_obj.data)
    meta = our_obj.meta
    common_probes = anno.index.intersection(dat.index)
    dat = dat.reindex(common_probes)
    anno = anno.reindex(common_probes)

    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method

    the_hash = tsgd.dmr_results_hash(meta.sort_index().index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_cross_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        dmr_res = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        logger.info("Unable to locate pre-existing results. Computing from scratch (this can take a while).")
        comparisons = {}
        for pid in pids:  # GIC
            for pid2 in pids:  # iNSC
                gics = meta.index[(meta.type == 'GBM') & (meta.patient_id == pid)]
                inscs = meta.index[(meta.type == 'iNSC') & (meta.patient_id == pid2)]
                comparisons['-'.join([pid, pid2])] = [gics, inscs]
        dmr_res = run_dmr_analyses(dat, comparisons, anno, dmr_params)
        # Save DMR results to disk
        dmr_res.to_pickle(fn, include_annotation=False)
        logger.info("Saved DMR results to %s", fn)

    dmr_res_all = dmr_res.results_significant

    # number of DMRs, split by direction

    n_dmr = pd.DataFrame(index=pids, columns=pids)
    n_dmr.index.name = 'GIC'
    n_dmr.columns.name = 'iNSC'
    n_hyper = n_dmr.copy()
    n_hypo = n_dmr.copy()

    for k, v in dmr_res_all.items():
        p1, p2 = k.split('-')
        n_dmr.loc[p1, p2] = len(v)
        n_hyper.loc[p1, p2] = len([t for t in v.values() if t['median_change'] > 0])
        n_hypo.loc[p1, p2] = len([t for t in v.values() if t['median_change'] < 0])

    n_hyper_pct= n_hyper / n_dmr * 100.

    # box whisker plot

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = sns.boxplot(n_hyper_pct.transpose(), whis=5, ax=ax)
    ax.scatter(range(len(pids)), np.diagonal(n_hyper_pct), facecolor='k')
    ax.set_ylim([0, 100])

    # run down the rows or columns and generate an 'overlap spectrum' for each one
    # rows: check the effect of varying the iNSC line (CONSISTENCY)
    # cols: check the effect of varying the GIC line (non-syngeneic DIFFERENCE)
    # also repeat for the columns, which is just the S1 approach (SYNGENEIC)

    row_collapse = pd.DataFrame(
        dict([
            (
                pid,
                setops.quantify_feature_membership(
                    setops.venn_from_arrays(
                        *[dmr_res_all['%s-%s' % (pid, p)].keys() for p in pids]
                    )[1]
                )
            )
            for pid in pids
        ])
    )[pids]

    col_collapse = pd.DataFrame(
        dict([
            (
                pid,
                setops.quantify_feature_membership(
                    setops.venn_from_arrays(
                        *[dmr_res_all['%s-%s' % (p, pid)].keys() for p in pids]
                    )[1]
                )
            )
            for pid in pids
        ])
    )[pids]

    syn_dist = setops.quantify_feature_membership(
        setops.venn_from_arrays(
            *[dmr_res_all['%s-%s' % (p, p)].keys() for p in pids]
        )[1]
    )