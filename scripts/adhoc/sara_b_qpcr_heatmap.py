from plotting import utils, heatmap
import pandas as pd
import os
from settings import DATA_DIR
from matplotlib import pyplot as plt

SAVEFIG = True

def fixed_plot_fun(data, cmap='RdBu'):
    ax, cax = heatmap.single_heatmap(
        data,
        vmin=0,
        vmax=20.,
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
    if SAVEFIG:
        OUTDIR = 'temp_results.0'
        i = 1
        while os.path.exists(OUTDIR):
            OUTDIR = 'temp_results.%d' % i
            i += 1
        print "Creating temp output dir %s" % OUTDIR
        os.makedirs(OUTDIR)

    infile = os.path.join(DATA_DIR, 'sara_b_qpcr', 'tumour_xenograft_vs_healthy.csv')
    data = pd.read_csv(infile, header=0, index_col=0).transpose()

    ax, cax = fixed_plot_fun(data, cmap='RdBu')
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_all.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_all.pdf"), dpi=200)

    ax, cax = fixed_plot_fun(data, cmap='RdBu_r')
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_all_reverse.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_all_reverse.pdf"), dpi=200)

    ax, cax = fixed_plot_fun(data.iloc[:, 3:], cmap='RdBu')
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_restricted.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_restricted.pdf"), dpi=200)

    ax, cax = fixed_plot_fun(data.iloc[:, 3:], cmap='RdBu_r')
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_restricted_reverse.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_restricted_reverse.pdf"), dpi=200)

    # same but plotting 20 - \Delta C_t so that "high is high"
    ax, cax = fixed_plot_fun(20. - data, cmap='RdBu_r')
    cax.ax.set_title('Expression', fontsize=14, position=[1., 1.03])
    cax.set_ticks([0, 20.])
    cax.set_ticklabels(['Low', 'High'])
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(OUTDIR, "qpcr_expr_heatmap_all.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "qpcr_expr_heatmap_all.pdf"), dpi=200)

    ax, cax = fixed_plot_fun(20. - data.iloc[:, 3:], cmap='RdBu_r')
    cax.ax.set_title('Expression', fontsize=14, position=[1., 1.03])
    cax.set_ticks([0, 20.])
    cax.set_ticklabels(['Low', 'High'])
    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(OUTDIR, "qpcr_expr_heatmap_restricted.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "qpcr_expr_heatmap_restricted.pdf"), dpi=200)

