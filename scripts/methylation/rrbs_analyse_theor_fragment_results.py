import os
import pandas as pd
import pysam
import numpy as np
from settings import OUTPUT_DIR
import pickle
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("dark")


def compl_ecdf(cpg_cov_all, values):
    cov = []
    ecdf = []
    for val in values:
        cov.append(val)
        ecdf.append((cpg_cov_all >= val).sum())
    return np.array(cov), np.array(ecdf)



if __name__ == "__main__":
    indir = os.path.join(OUTPUT_DIR, 'rrbs_theor_fragments.0')
    # NB: plots will go in the same directory!
    outdir = indir
    fnames = ['GC-CV-7163-{0}_S{0}.coverage.pkl'.format(i) for i in range(1, 7)]
    sample_names = [
        'eNSC3',
        'eNSC5',
        'eNSC6',
        'mDura3Human',
        'mDura5Human',
        'mDura6Human',
    ]
    res = {}
    for fname, sname in zip(fnames, sample_names):
        fn = os.path.join(indir, fname)
        with open(fn, 'rb') as f:
            res[sname] = pickle.load(f)


    for sname in res:
        x = res[sname]
        n_cpg = float(sum([len(v) for v in x['cpg_cov'].values()]))

        cpg_cov_all_nz = np.array(reduce(lambda x, y: x + y, [[t for t in z if t > 0] for z in x['cpg_cov'].values()]))

        # distribution of CpG coverage: low coverage region
        cov, ecdf = compl_ecdf(cpg_cov_all_nz, range(1, 31))

        # inverse CDF of low coverage region
        fig = plt.figure(figsize=(8.5, 5), num="%s_low" % sname)
        ax1 = fig.add_subplot(111)
        ax1.bar(cov, ecdf / n_cpg * 100)
        ax1.set_xticks(cov)
        ax1.set_xlabel('Minimum coverage')
        ax1.set_ylabel('% CpG sites')
        ax2 = ax1.twinx()
        h = ax2.plot(cov, ecdf / 1e6, 'x')
        ax2.set_ylim(np.array(ax1.get_ylim()) / 100 * n_cpg / 1e6)
        ax2.set_ylabel("Number of CpG sites (millions)")
        h[0].set_visible(False)
        fig.tight_layout()

        fig.savefig(os.path.join(outdir, "%s_cpg_coverage_low.png" % sname), dpi=200)

        # distribution of CpG coverage: high coverage region
        cov, ecdf = compl_ecdf(cpg_cov_all_nz, [20, 30, 40, 50, 60, 70, 80, 90, 100, 200])

        fig = plt.figure(figsize=(8.5, 5), num="%s_high" % sname)
        ax1 = fig.add_subplot(111)
        ax1.bar(range(len(cov)), ecdf / n_cpg * 100)
        ax1.set_xticks(range(len(cov)))
        ax1.set_xticklabels(cov)
        ax1.set_xlabel('Minimum coverage')
        ax1.set_ylabel('% CpG sites')
        ax2 = ax1.twinx()
        h = ax2.plot(range(len(cov)), ecdf / 1e3, 'x')
        ax2.set_ylim(np.array(ax1.get_ylim()) / 100 * n_cpg / 1e3)
        ax2.set_ylabel("Number of CpG sites (thousands)")
        h[0].set_visible(False)
        fig.tight_layout()

        fig.savefig(os.path.join(outdir, "%s_cpg_coverage_high.png" % sname), dpi=200)

        # coverage vs fragment size for predicted fragments
        chroms = sorted(x['fragment_coords'], key=int)
        fsize = np.concatenate(
            [np.diff(x['fragment_coords'][t]) for t in chroms]
        )
        fcov = np.concatenate([np.array(x['fragment_coverage'][t])[:, 1] for t in chroms])

