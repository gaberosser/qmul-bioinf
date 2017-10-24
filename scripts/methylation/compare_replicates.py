from load_data import methylation_array
from methylation import process
from stats import transformations
import pandas as pd
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
from utils import output
import os


def scatter_plots(dat, pids, mad=None, top_n_by_mad=None):
    fig, axs = plt.subplots(nrows=2, ncols=len(pids), sharex=True, sharey=True, figsize=(11, 7))

    for i, pid in enumerate(pids):
        for j, lbl in enumerate(['GBM', 'DURA']):
            cols = dat.columns[dat.columns.str.contains(lbl + pid)]
            if len(cols) != 2:
                axs[j, i].set_title("n = 1")
                continue
            if top_n_by_mad is None:
                axs[j, i].scatter(dat.loc[:, cols[0]], dat.loc[:, cols[1]])
                axs[j, i].set_title("%s%s r=%.3f" % (
                    lbl, pid, stats.linregress(
                        dat.loc[:, cols[0]], dat.loc[:, cols[1]]
                    ).rvalue
                ))
            else:
                if mad is None:
                    raise AttributeError("Must supply `mad` if `top_n_by_mad` is not None.")
                idx = mad.sort_values(ascending=False).index[:top_n_by_mad]
                axs[j, i].scatter(dat.loc[idx, cols[0]], dat.loc[idx, cols[1]])
                axs[j, i].set_title("%s%s r=%.3f" % (
                    lbl, pid, stats.linregress(
                        dat.loc[idx, cols[0]], dat.loc[idx, cols[1]]
                    ).rvalue
                ))
    fig.tight_layout()
    return fig, axs


if __name__ == '__main__':
    outdir = output.unique_output_dir("methylation_replicates", reuse_empty=True)
    pids = ['017', '050', '054', '061']
    b, me_meta = methylation_array.load_by_patient(pids)
    m = process.m_from_beta(b)
    mad = transformations.median_absolute_deviation(m).sort_values(ascending=False)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(mad.values)
    ax.set_xlabel("Probe rank by MAD")
    ax.set_ylabel("MAD value")
    ax.axvline(50000, ls='--', c='r')
    fig.savefig(os.path.join(outdir, "MAD_sorted.png"), dpi=200)

    fig1, axs1 = scatter_plots(m, pids)
    fig1.savefig(os.path.join(outdir, "correlation_all_probes.png"), dpi=200)
    fig2, axs2 = scatter_plots(m, pids, mad=mad, top_n_by_mad=50000)
    fig2.savefig(os.path.join(outdir, "correlation_top_50000.png"), dpi=200)
    fig3, axs3 = scatter_plots(m, pids, mad=mad, top_n_by_mad=100000)