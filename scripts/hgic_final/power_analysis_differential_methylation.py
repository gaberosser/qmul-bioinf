import multiprocessing as mp
import os
import pandas as pd
import numpy as np
import pickle
from utils import output, setops, excel, log
from methylation import dmr, process, loader as methylation_loader, annotation_gene_to_ensembl
from rnaseq import loader as rnaseq_loader, differential_expression, general, filter
from integrator import rnaseq_methylationarray
from analysis import cross_comparison
from load_data import loader
from plotting import venn
from matplotlib import pyplot as plt, text, patches
import seaborn as sns
import statsmodels.stats.power as smp
from scripts.hgic_final import consts, two_strategies_grouped_dispersion as tsgd

logger = log.get_console_logger()


if __name__ == "__main__":
    outdir = output.unique_output_dir()

    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    norm_method_s1 = 'swan'

    pids = consts.PIDS
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    me_obj, anno = tsgd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
    me_data = me_obj.data
    me_meta = me_obj.meta

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    the_hash = tsgd.dmr_results_hash(me_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise NotImplementedError("Require pre-computed DMR results (S1)")

    # joint distribution of cluster sizes and std among them (one per sample)
    cluster_sizes = pd.Series(dict([
        (k, len(v.pids)) for k, v in dmr_res_s1.clusters.items()
    ]))
    cluster_size_dist = cluster_sizes.value_counts()
    n_sample = me_data.columns.size
    n_cluster = len(dmr_res_s1.clusters)
    n_std = 20  # number of stdev bins to use
    stdev_bins = np.linspace(0, 10, n_std)

    res = {}

    for j, cid in enumerate(dmr_res_s1.clusters.keys()):
        this_pids = dmr_res_s1.clusters[cid].pids
        n = len(this_pids)
        res.setdefault(n, {})
        this_dat = me_data.loc[this_pids]
        this_std = this_dat.std(axis=0)
        for i, sname in enumerate(me_data.columns):
            res[n].setdefault(sname, []).append(this_std[sname])

    aa = pd.DataFrame([[np.mean(res[i][k]) for k in me_data.columns] for i in sorted(res.keys())], columns=me_data.columns, index=sorted(res.keys()))
    mean_std = aa.mean(axis=1)

    # we know:
    # the minimum effect size (1.4)
    # the desired FDR value (0.01)
    # the distribution of the number of probes
    # the distribution of stdev within those probes

    # weighted estimate: power(cluster size 6) * n(cluster size 6) + ...
    # simplify by assuming a t test for the power calculation

    fdr = dmr_params['alpha']
    delta_m = dmr_params['delta_m_min']
    power_calc = {}
    for clen in aa.index:
        eff_size = delta_m / mean_std[clen]
        power_calc[clen] = smp.tt_ind_solve_power(effect_size=eff_size, nobs1=clen * 2, ratio=1., alpha=fdr)

    weighted_mean = (pd.Series(power_calc) * cluster_size_dist).sum() / float(n_cluster)

