import os
import pandas as pd
import pysam
import numpy as np
from settings import OUTPUT_DIR, GIT_LFS_DATA_DIR
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


def multi_grouped_boxplot(data, sharex=True, sharey=True, subplot_titles=None, w_pad=0.2, figsize=None, **kwargs):
    """
    Generate a grouped boxplot with multiple observations per group separated
    :param mat: Data to be plotted.
    Each element of args is a list of lists representing multiple observations for one group. Every list should have the
    same length.
    :param kwargs:
    :return:
    """
    if subplot_titles is not None:
        if len(subplot_titles) != len(data):
            raise AttributeError("len of subplot_titles must match len of data")

    fig, axs = plt.subplots(ncols=len(data), sharey=sharey, sharex=sharex, figsize=figsize)
    big_ax = fig.add_subplot(111, frameon=False)
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.grid(False)

    for i in range(len(data)):
        ax = axs[i]
        ax.boxplot(data[i], **kwargs)
        if subplot_titles is not None:
            ax.set_title(subplot_titles[i])

    # set shared maximum
    ylims = np.array([ax.get_ylim() for ax in axs])
    for ax in axs:
        ax.set_ylim([ylims[:, 0].min(), ylims[:, 1].max()])

    fig.tight_layout(w_pad=w_pad)

    return fig, axs, big_ax


if __name__ == "__main__":
    indir = os.path.join(OUTPUT_DIR, 'rrbs_theor_fragments')
    # NB: plots will go in the same directory!
    outdir = indir
    fnames = ['GC-CV-7163-{0}_S{0}'.format(i) for i in range(1, 7)]
    sample_names = [
        'eNSC3',
        'eNSC5',
        'eNSC6',
        'mDura3Human',
        'mDura5Human',
        'mDura6Human',
    ]
    fragment_length_bins = [0, 100, 500, 1000, 2000, 1e9]
    process_cpg = True
    process_ccgg = True

    fcov_binned = fcov_normed_binned = flen_binned = None

    # fcov_binned = [pd.DataFrame(columns=sample_names) for i in range(len(fragment_length_bins) - 1)]
    # fcov_norm_binned = [pd.DataFrame(columns=sample_names) for i in range(len(fragment_length_bins) - 1)]
    # flen_binned = [pd.DataFrame(columns=sample_names) for i in range(len(fragment_length_bins) - 1)]

    for fname, sname in zip(fnames, sample_names):
        coverage_fn = os.path.join(indir, "%s.cpg_coverage.gz" % fname)
        counts_fn = os.path.join(indir, "%s.mspi_fragments.counts" % fname)

        ccgg_coverage = pd.read_csv(counts_fn, sep='\t', comment='#', header=0, index_col=0)
        cpg_depth = pd.read_csv(coverage_fn, sep='\t', header=None, index_col=None)
        cpg_depth.columns = ['chr', 'coord', 'depth']

        if process_cpg:

            n_cpg = float(cpg_depth.shape[0])
            cpg_cov_all_nz = cpg_depth.loc[cpg_depth.depth > 0]

            # distribution of CpG coverage: low coverage region
            cov, ecdf = compl_ecdf(cpg_cov_all_nz.depth, range(1, 31))

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
            cov, ecdf = compl_ecdf(cpg_cov_all_nz.depth, [20, 30, 40, 50, 60, 70, 80, 90, 100, 200])

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

        if process_ccgg:

            flen = ccgg_coverage.Length
            fcov = ccgg_coverage.iloc[:, -1]  # column name is full path to file, so unwieldy

            # bin the data
            flen_binned_idx = [
                ((flen >= fragment_length_bins[i]) & (flen < fragment_length_bins[i+1]))
                for i in range(len(fragment_length_bins) - 1)
            ]

            if fcov_binned is None:
                # set all now
                fcov_binned = [pd.DataFrame(columns=sample_names, index=flen.index[t]) for t in flen_binned_idx]
                fcov_normed_binned = [pd.DataFrame(columns=sample_names, index=flen.index[t]) for t in flen_binned_idx]
                flen_binned = [flen[t] for t in flen_binned_idx]

            for i, t in enumerate(flen_binned_idx):
                fcov_binned[i].loc[fcov[t].index, sname] = fcov[t]
                fcov_normed_binned[i].loc[fcov[t].index, sname] = fcov[t] / flen[t].astype(float)

    if process_ccgg:
        dat = [[t.loc[:, sname].values for t in fcov_binned] for sname in sample_names]
        fig, axs, big_ax = multi_grouped_boxplot(dat, subplot_titles=sample_names, showfliers=False, figsize=(10.2, 4.4))
        x_lbls = [
            "%d - %d" % (
                fragment_length_bins[i],
                fragment_length_bins[i+1] - 1
            )
            if fragment_length_bins[i+1] < 1e6
            else ">= %d" % fragment_length_bins[i]
            for i in range(len(fragment_length_bins) - 1)
        ]
        for ax in axs:
            ax.set_xticklabels(x_lbls, rotation=90)
        big_ax.set_ylabel('Coverage')
        big_ax.set_xlabel('Theoretical fragment size', labelpad=50)
        fig.subplots_adjust(bottom=0.3)
        fig.savefig(os.path.join(outdir, "mspi_theoretical_fragment_coverage_binned.png"), dpi=200)
