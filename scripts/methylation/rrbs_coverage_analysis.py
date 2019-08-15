import json
import os
import gzip
import multiprocessing as mp
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from utils import log, output
from settings import INTERMEDIATE_DIR
from StringIO import StringIO
clogger = log.get_console_logger("rrbs_coverage_analysis")


def parse_one_result(cpg, perms):
    if len(cpg) == 0:
        return

    tab = pd.read_csv(StringIO(cpg), sep='\t', header=0, dtype=int)
    cpg_chr = tab.iloc[0, 0]
    cpg_start = tab.iloc[0, 1]
    cpg_arr = tab.iloc[:, -1].values

    # permutations
    n_bp = len(cpg_arr)
    n_perm = len(perms)
    perm_arr = np.zeros((n_bp, n_perm))
    perm_loc = []
    for j in range(n_perm):
        if len(perms[j]) == 0:
            # no coverage at all - skip
            continue
        tab_perm = pd.read_csv(
            StringIO(perms[j]), sep='\t', header=0, dtype=int
        )
        the_chr = tab_perm.iloc[0, 0]
        the_start = tab_perm.iloc[0, 1]
        the_arr = tab_perm.iloc[:, -1].values
        perm_loc.append((the_chr, the_start))
        if len(the_arr) == n_bp:
            perm_arr[:, j] = the_arr
        else:
            perm_arr[:len(the_arr), j] = the_arr
            perm_arr[len(the_arr):, j] = np.nan

    return (cpg_chr, cpg_start), cpg_arr, perm_loc, perm_arr


if __name__ == "__main__":
    outdir = output.unique_output_dir()
    indir = os.path.join(INTERMEDIATE_DIR, 'rrbs_coverage_sampling', 'GC-CV-7163')
    flist = ['GC-CV-7163-{n}_S{n}.cpg_coverage.json.gz'.format(n=i) for i in range(1, 7)]
    titles = [
        'eNSC3',
        'eNSC5',
        'eNSC6',
        'mDura iNSC3',
        'mDura iNSC5',
        'mDura iNSC6',
    ]

    cpg_island_cov = {}
    cpg_island_loc = {}

    cpg_island_perm = {}
    cpg_island_perm_loc = {}

    for fname in flist:
        fn = os.path.join(indir, fname)
        clogger.info("Loading JSON file %s", fn)
        try:
            with gzip.open(fn, 'rb') as f:
                res = json.load(f)
            clogger.info("File %s has %d results.", fname, len(res['cpg_islands']))
        except Exception:
            clogger.exception("Failed to load file %s", fn)
            continue
        cpg_island_cov[fname] = {}
        cpg_island_loc[fname] = {}
        cpg_island_perm[fname] = {}
        cpg_island_perm_loc[fname] = {}

        pool = mp.Pool()
        jobs = {}

        this_cpg = cpg_island_cov[fname]
        this_loc = cpg_island_loc[fname]
        this_perm = cpg_island_perm[fname]
        this_perm_loc = cpg_island_perm_loc[fname]

        n_skipped = 0

        for i in range(len(res['cpg_islands'])):
            jobs[i] = pool.apply_async(
                parse_one_result,
                args=(res['cpg_islands'][i], res['permutations'][i])
            )
        pool.close()
        for i, j in jobs.items():
            this_res = j.get(1e2)
            if this_res is None:
                n_skipped += 1
            else:
                cpg_loc, cpg_arr, perm_loc, perm_arr = this_res
                this_cpg[i] = cpg_arr
                this_loc[i] = cpg_loc
                this_perm_loc[i] = perm_loc
                this_perm[i] = perm_arr

    # create an ensemble trace for each file
    n_pt = 200
    x_pt = np.linspace(0, 1, n_pt)

    fig, axs = plt.subplots(2, 3, figsize=(16, 8.5))
    mean_cov_trace = {}
    mean_perm_trace = {}

    for i, (fname, ttl) in enumerate(zip(flist, titles)):
        ax = axs[0, i] if i < 3 else axs[1, i - 3]

        n = len(cpg_island_cov[fname])
        cov_trace = np.zeros(n_pt)
        cov_sq_trace = np.zeros(n_pt)
        cov_ct = 0
        perm_trace = np.zeros(n_pt)
        perm_sq_trace = np.zeros(n_pt)
        perm_ct = 0

        for yp in cpg_island_cov[fname].values():
            xi = (yp.size - 1) * x_pt
            xp = np.arange(yp.size)
            yi = np.interp(xi, xp, yp)
            cov_trace += yi
            cov_sq_trace += yi ** 2
            cov_ct += 1

        for parr in cpg_island_perm[fname].values():
            xi = (parr.shape[0] - 1) * x_pt
            xp = np.arange(parr.shape[0])
            for i in range(parr.shape[1]):
                yp = parr[:, i]
                yi = np.interp(xi, xp, yp)
                perm_trace += yi
                perm_sq_trace += yi ** 2
                perm_ct += 1

        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        ax.plot(x_pt, cov_trace / float(cov_ct), label='CpG island')
        ax.plot(x_pt, perm_trace / float(perm_ct), label='Random region')
        ax.legend(loc='upper right')
        ax.set_xlabel('Normalised coordinate')
        ax.set_ylabel('Mean coverage (read / bp)')
        ax.set_title(ttl)

        mean_cov_trace[ttl] = cov_trace / float(cov_ct)
        mean_perm_trace[ttl] = perm_trace / float(perm_ct)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "mean_coverage_traces.png"), dpi=200)


    # repeat but with mean relative enrichment

    fig, axs = plt.subplots(2, 3, figsize=(16, 8.5))

    for i, (fname, ttl) in enumerate(zip(flist, titles)):
        ax = axs[0, i] if i < 3 else axs[1, i - 3]
        m = mean_cov_trace[ttl] / mean_perm_trace[ttl]
        ax.plot(x_pt, m)
        ax.set_xlabel('Normalised coordinate')
        ax.set_ylabel('Relative coverage enrichment in CpG islands')
        ax.set_title(ttl)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "relative_coverage_enrichment_traces.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "relative_coverage_enrichment_traces.png"), dpi=200)

