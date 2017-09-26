from plotting import common, heatmap
from utils.output import unique_output_dir
import pandas as pd
import os
from settings import DATA_DIR
from matplotlib import pyplot as plt


def fixed_plot_fun(data, cmap='RdBu', vmin=0., vmax=20.):
    ax, cax = heatmap.single_heatmap(
        data,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        fig_kwargs={'figsize': [6, 5.5]}
    )
    cax.ax.set_title('$\Delta C_t$', fontsize=16, position=[1., 1.03])
    plt.setp(ax.xaxis.get_ticklabels(), fontsize=14)
    plt.setp(ax.yaxis.get_ticklabels(), fontsize=14)
    plt.setp(cax.ax.yaxis.get_ticklabels(), fontsize=14)
    plt.tight_layout()

    return ax, cax


if __name__ == "__main__":
    SAVEFIG = True

    if SAVEFIG:
        outdir = unique_output_dir('sb_qpcr_heatmap')

    infile = os.path.join(DATA_DIR, 'sara_b_qpcr', 'tumour_xenograft_vs_healthy2.csv')
    data = pd.read_csv(infile, header=0, index_col=0).transpose()

    ax, cax = fixed_plot_fun(data, cmap='RdBu')
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(outdir, "qpcr_heatmap_all.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "qpcr_heatmap_all.pdf"), dpi=200)

    ax, cax = fixed_plot_fun(data, cmap='RdBu_r')
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(outdir, "qpcr_heatmap_all_reverse.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "qpcr_heatmap_all_reverse.pdf"), dpi=200)

    ax, cax = fixed_plot_fun(data.iloc[:, 3:], cmap='RdBu')
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(outdir, "qpcr_heatmap_restricted.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "qpcr_heatmap_restricted.pdf"), dpi=200)

    ax, cax = fixed_plot_fun(data.iloc[:, 3:], cmap='RdBu_r')
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(outdir, "qpcr_heatmap_restricted_reverse.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "qpcr_heatmap_restricted_reverse.pdf"), dpi=200)

    # same but plotting 20 - \Delta C_t so that "high is high"
    a = 20. - data

    ax, cax = fixed_plot_fun(a, cmap='RdBu_r')
    cax.ax.set_title('Expression', fontsize=14, position=[1., 1.03])
    cax.set_ticks([0, 20.])
    cax.set_ticklabels(['Low', 'High'])
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(outdir, "qpcr_expr_heatmap_all.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "qpcr_expr_heatmap_all.pdf"), dpi=200)

    ax, cax = fixed_plot_fun(a.iloc[:, 3:], cmap='RdBu_r')
    cax.ax.set_title('Expression', fontsize=14, position=[1., 1.03])
    cax.set_ticks([0, 20.])
    cax.set_ticklabels(['Low', 'High'])
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(outdir, "qpcr_expr_heatmap_restricted.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "qpcr_expr_heatmap_restricted.pdf"), dpi=200)

    # standardise per gene
    m = a.mean(axis=1)
    s = a.std(axis=1)
    za = a.subtract(m, axis=0).divide(s, axis=0)

    ax, cax = fixed_plot_fun(za, cmap='RdBu_r', vmin=-1.5, vmax=1.5)
    cax.ax.set_title('Rel. expression', fontsize=14, position=[.5, 1.03])
    cax.set_ticks([-1.5, 1.5])
    cax.set_ticklabels(['Low', 'High'])
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(outdir, "qpcr_norm_expr_heatmap_all.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "qpcr_norm_expr_heatmap_all.pdf"), dpi=200)

    b = a.iloc[:, 3:]
    m = b.mean(axis=1)
    s = b.std(axis=1)
    zb = b.subtract(m, axis=0).divide(s, axis=0)

    ax, cax = fixed_plot_fun(zb, cmap='RdBu_r', vmin=-1.5, vmax=1.5)
    cax.ax.set_title('Rel. expression', fontsize=14, position=[.5, 1.03])
    cax.set_ticks([-1.5, 1.5])
    cax.set_ticklabels(['Low', 'High'])
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(outdir, "qpcr_norm_expr_heatmap_restricted.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "qpcr_norm_expr_heatmap_restricted.pdf"), dpi=200)