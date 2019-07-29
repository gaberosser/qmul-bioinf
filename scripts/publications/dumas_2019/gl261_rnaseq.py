from utils import output, excel
from plotting import common, scatter, clustering
from sklearn.decomposition import PCA
from stats import transformations
from rnaseq import loader, differential_expression, filter, general

import os
import pandas as pd
import collections
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from scripts.hgic_final.cluster_lines_methylation import plot_pca


if __name__ == '__main__':

    # parameters
    alpha = 0.05
    min_logfc = 0
    eps = 0.01
    n_by_mad = 3000
    min_cpm = 0.01  # for filtering purposes
    fontsize = 12

    treatment_colour = {
        'WT': '0.2',
        'Rheb KO': '0.8'
    }

    outdir = output.unique_output_dir()

    # load our data
    obj_star = loader.load_references('wtchg_p190202', alignment_subdir='mouse', tax_id=10090)
    obj_salmon = loader.load_references('wtchg_p190202', alignment_subdir='mouse', source='salmon', tax_id=10090)

    # dump to file for sharing
    dat = obj_salmon.data.copy()
    general.add_gene_symbols_to_ensembl_data(dat, tax_id=10090)
    dat.to_excel(os.path.join(outdir, 'salmon_tpm_all_data.xlsx'))
    dat = obj_star.data.copy()
    general.add_gene_symbols_to_ensembl_data(dat, tax_id=10090)
    dat.to_excel(os.path.join(outdir, 'star_counts_all_data.xlsx'))

    # In our original analysis (see anaelle/analyse_mouse_gl261_rnaseq.py) we noted that sample WT76 is possibly an
    # outlier, so we decided to exclude it. Hmm, not sure that's very good science...

    ix = ~obj_star.data.columns.isin(['WT_76'])
    dat = filter.filter_by_cpm(obj_star.data.loc[:, ix], min_cpm=min_cpm, min_n_samples=2)

    # run DE analysis (3 vs 2)

    the_groups = obj_star.meta.treatment.str.replace('Rheb KO', 'Rheb_KO')[ix]  # group names must be valid in R
    the_comparison = ('Rheb_KO', 'WT')
    de_res = differential_expression.run_one_de(
        dat,
        the_groups,
        the_comparison,
        tax_id=10090,
        method='GLM',
        fdr=alpha,
        lfc=min_logfc,
        return_full=True
    )
    de_res.to_csv(os.path.join(outdir, "de_res_minus_wt76_full.csv"))

    # PCA
    meta = obj_salmon.meta.loc[obj_salmon.meta.index != 'WT_76']
    dat = obj_salmon.data.loc[:, meta.index]
    log_dat = np.log10(dat + eps)
    colour_subgroups = meta.treatment

    # change text for plot
    colour_subgroups = colour_subgroups.replace('WT', 'fl/fl')
    colour_subgroups = colour_subgroups.replace('Rheb KO', r'$\Delta$/$\Delta$')

    cmap = collections.OrderedDict([
        ('fl/fl', 'k'),
        (r'$\Delta$/$\Delta$', 'w')
    ])

    marker_subgroups = colour_subgroups
    mmap = collections.OrderedDict([
        ('fl/fl', 'o'),
        (r'$\Delta$/$\Delta$', 's')
    ])

    p = PCA()
    pc_dat = p.fit_transform(log_dat.transpose())
    variance_explained = p.explained_variance_ratio_ * 100.

    ax = scatter.scatter_with_colour_and_markers(
        pc_dat[:, [0, 1]],
        colour_subgroups=None,
        colour_map=None,
        marker_subgroups=marker_subgroups,
        marker_map=mmap,
        default_colour='k',
        ms=50
    )
    ax.set_xlabel("PCA component %s (%.1f%%)" % (1, variance_explained[0]), fontsize=fontsize)
    ax.set_ylabel("PCA component %s (%.1f%%)" % (2, variance_explained[1]), fontsize=fontsize)

    leg = ax.legend()
    fr = leg.get_frame()
    fr.set_edgecolor('k')
    plt.setp(leg.get_texts(), fontsize=fontsize)
    plt.setp(ax.xaxis.get_ticklabels(), fontsize=fontsize)
    plt.setp(ax.yaxis.get_ticklabels(), fontsize=fontsize)

    ax.figure.set_size_inches(5., 4.8)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "pca_gl261_rnaseq.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "pca_gl261_rnaseq.tiff"), dpi=200)

    # clustermap
    colour_bar = pd.DataFrame(treatment_colour['Rheb KO'], index=dat.columns, columns=[''])
    colour_bar.loc[meta.treatment == 'WT'] = treatment_colour['WT']

    this_mad = transformations.median_absolute_deviation(log_dat).sort_values(ascending=False)
    this_log_dat = log_dat.loc[this_mad.index[:n_by_mad]]

    # version 1: cluster rows
    cm = clustering.plot_clustermap(
        this_log_dat,
        cmap='RdYlBu_r',
        metric='correlation',
        col_colors=colour_bar,
        vmin=-2,
        vmax=2,
    )

    leg_dict = {
        'fl/fl': {
            'class': 'patch',
            'edgecolor': 'none',
            'facecolor': treatment_colour['WT'],
        },
        r'$\Delta$/$\Delta$': {
            'class': 'patch',
            'edgecolor': 'none',
            'facecolor': treatment_colour['Rheb KO'],
        },
    }
    common.add_custom_legend(cm.ax_heatmap, leg_dict, loc_outside=True, fontsize=fontsize)
    leg = cm.ax_heatmap.get_legend()
    leg.set_frame_on(False)
    cm.ax_heatmap.xaxis.set_ticklabels([])

    cm.fig.set_size_inches((5, 6.))
    cm.gs.update(bottom=0.02, right=0.78)
    cm.savefig(os.path.join(outdir, "clustermap_gl261_rnaseq.png"), dpi=200)
    cm.savefig(os.path.join(outdir, "clustermap_gl261_rnaseq.tiff"), dpi=200)


    # manually change the order of the yaxis (genes) so that it reflects mean across KO (descending)
    ix = this_log_dat.loc[:, meta.treatment == 'Rheb KO'].sum(axis=1).sort_values(ascending=False).index

    cm = clustering.plot_clustermap(
        this_log_dat.loc[ix],
        cmap='RdYlBu_r',
        metric='correlation',
        col_colors=colour_bar,
        vmin=-2,
        vmax=2,
        row_cluster=False
    )

    leg_dict = {
        'fl/fl': {
            'class': 'patch',
            'edgecolor': 'none',
            'facecolor': treatment_colour['WT'],
        },
        r'$\Delta$/$\Delta$': {
            'class': 'patch',
            'edgecolor': 'none',
            'facecolor': treatment_colour['Rheb KO'],
        },
    }
    common.add_custom_legend(cm.ax_heatmap, leg_dict, loc_outside=True, fontsize=fontsize)
    leg = cm.ax_heatmap.get_legend()
    leg.set_frame_on(False)
    cm.ax_heatmap.xaxis.set_ticklabels([])

    cm.fig.set_size_inches((5, 6.))
    cm.gs.update(bottom=0.02, right=0.78)
    cm.savefig(os.path.join(outdir, "clustermap_gl261_rnaseq_ordered.png"), dpi=200)
    cm.savefig(os.path.join(outdir, "clustermap_gl261_rnaseq_ordered.tiff"), dpi=200)

