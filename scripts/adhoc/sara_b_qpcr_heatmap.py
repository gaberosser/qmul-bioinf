from plotting import utils, heatmap
import pandas as pd
import os
from settings import DATA_DIR
from matplotlib import pyplot as plt

SAVEFIG = True


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
    ax, cax = heatmap.single_heatmap(
        data,
        vmin=0,
        vmax=20.,
        cmap='RdBu',
        fig_kwargs={'figsize': [6, 5.5]}
    )
    cax.ax.set_title('$\Delta C_t$', fontsize=16, position=[1., 1.03])
    plt.setp(ax.xaxis.get_ticklabels(), fontsize=14)
    plt.setp(ax.yaxis.get_ticklabels(), fontsize=14)
    plt.setp(cax.ax.yaxis.get_ticklabels(), fontsize=14)
    plt.tight_layout()

    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_all.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_all.pdf"), dpi=200)

    ax, cax = heatmap.single_heatmap(
        data.iloc[:, 3:],
        vmin=0,
        vmax=20.,
        cmap='RdBu_r',
        fig_kwargs={'figsize': [6, 5.5]}
    )
    cax.ax.set_title('$\Delta C_t$', fontsize=16, position=[1., 1.03])
    plt.setp(ax.xaxis.get_ticklabels(), fontsize=14)
    plt.setp(ax.yaxis.get_ticklabels(), fontsize=14)
    plt.setp(cax.ax.yaxis.get_ticklabels(), fontsize=14)
    plt.tight_layout()

    if SAVEFIG:
        fig = ax.figure
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_restricted.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "qpcr_heatmap_restricted.pdf"), dpi=200)