from utils import output, excel
from plotting import common, rnaseq, clustering
from sklearn.decomposition import PCA
from stats import transformations
from rnaseq import loader, differential_expression, filter

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

    treatment_colour = {
        'WT': '0.8',
        'Rheb KO': '0.2'
    }

    outdir = output.unique_output_dir()

    # load our data
    obj_star = loader.load_references('wtchg_p190202', alignment_subdir='mouse', tax_id=10090)
    obj_salmon = loader.load_references('wtchg_p190202', alignment_subdir='mouse', source='salmon', tax_id=10090)

    # load Bowman data

    # plots with only our samples
    dat = filter.filter_by_cpm(obj_salmon.data, min_cpm=min_cpm, min_n_samples=2)
    log_dat = np.log10(obj_salmon.data + eps)

    # ECDF
    ax = rnaseq.log_cpm_ecdf_plot(dat, units='tpm', min_cpm=min_cpm)
    ax.figure.set_size_inches(6, 4)
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "cdf_our_samples.png"), dpi=200)

    # PCA
    colour_subgroups = obj_salmon.meta.treatment

    cmap = collections.OrderedDict(zip(
        colour_subgroups,
        common.get_best_cmap(len(colour_subgroups)),
    ))

    p = PCA()
    pc_dat = p.fit_transform(log_dat.transpose())

    p, ax = plot_pca(
        log_dat,
        colour_subgroups,
        colour_map=cmap,
        p=p
    )

    for i, col in enumerate(dat.columns):
        ax.text(pc_dat[i, 0], pc_dat[i, 1], col)

    ax.figure.set_size_inches(5.9, 4.8)
    ax.figure.subplots_adjust(right=0.8, left=0.12, bottom=0.1, top=0.98)
    ax.figure.savefig(os.path.join(outdir, "pca_our_samples.png"), dpi=200)

    # clustermap: just our samples
    colour_bar = pd.DataFrame(treatment_colour['Rheb KO'], index=dat.columns, columns=[''])
    colour_bar.loc[colour_subgroups == 'WT'] = treatment_colour['WT']

    this_mad = transformations.median_absolute_deviation(log_dat).sort_values(ascending=False)
    this_log_dat = log_dat.loc[this_mad.index[:n_by_mad]]
    cm = clustering.plot_clustermap(
        this_log_dat,
        cmap='RdYlBu_r',
        metric='correlation',
        col_colors=colour_bar,
        vmin=-2,
        vmax=2
    )
    cm.fig.set_size_inches(5, 8.4)
    cm.gs.update(bottom=0.15, right=0.98)
    cm.savefig(os.path.join(outdir, "clustermap_our_samples.png"), dpi=200)

    # DE 3 vs 3
    dat = filter.filter_by_cpm(obj_star.data, min_cpm=min_cpm, min_n_samples=2)

    the_groups = obj_star.meta.treatment.str.replace('Rheb KO', 'Rheb_KO')  # group names must be valid in R
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

    de_res.to_csv(os.path.join(outdir, "de_res_full.csv"))

    # try removing one Rheb case
    # no real reason to justify this at present
    if False:
        ix = obj_star.data.columns != 'KO_Rheb_74'
        dat = filter.filter_by_cpm(obj_star.data.loc[:, ix], min_cpm=min_cpm, min_n_samples=2)

        the_groups = obj_star.meta.treatment.str.replace('Rheb KO', 'Rheb_KO')[ix]  # group names must be valid in R
        the_comparison = ('Rheb_KO', 'WT')
        de_res_2 = differential_expression.run_one_de(
            dat,
            the_groups,
            the_comparison,
            tax_id=10090,
            method='GLM',
            fdr=alpha,
            lfc=min_logfc,
            return_full=True
        )
