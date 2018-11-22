import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import pandas as pd
import itertools
import collections
import os

from rnaseq import loader, filter
from scripts.hgic_final import consts
from plotting import common
from utils import output
import references

if __name__ == "__main__":
    """
    Idea here: recreate the analysis Sven carried out, generating biplots for the RNA-Seq data.
    We can then extend this idea to methylation (?)

    Absolutely amazing resource here:
    https://www.fbbva.es/microsite/multivariate-statistics/biplots.html

    Code snippet inspiration here:
    https://stackoverflow.com/questions/39216897/plot-pca-loadings-and-loading-in-biplot-in-sklearn-like-rs-autoplot
    """
    outdir = output.unique_output_dir()

    eps = 1.  # offset applied during log transform

    # load data for iNSC and GBM (Salmon TPM)
    obj = loader.load_by_patient(consts.PIDS, source='salmon')
    obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES_GIC + consts.S1_RNASEQ_SAMPLES_INSC)

    scatter_colours = {
        'GBM': 'red',
        'iNSC': 'blue'
    }

    cmap = common.get_best_cmap(len(consts.PIDS))
    scatter_colours = dict(zip(consts.PIDS, cmap))

    scatter_markers = {
        'GBM': 's',
        'iNSC': 'o'
    }

    dat = filter.filter_by_cpm(obj.data, min_n_samples=2)
    # TODO: include VST or similar here
    dat = np.log(dat + eps)

    n = float(dat.shape[1])

    # standardise: mean centred data required for sensible decomposition
    # standardisation occurs along the FEATURES axis, which is dim 1
    scaler = StandardScaler(with_std=False)

    # features on the COLS, mean centre by PATIENT
    # scaler.fit(dat)
    # X = scaler.transform(dat).transpose()

    # features on the ROWS, mean centre by PATIENT
    # scaler.fit(dat)
    # X = scaler.transform(dat)

    # features on the ROWS, mean centre by gene
    scaler = scaler.fit(dat.transpose())
    X = scaler.transform(dat.transpose()).transpose()

    # features on the COLS, mean centre by gene
    # scaler.fit(dat.transpose())
    # X = scaler.transform(dat.transpose())

    # SVD
    u, s, vh = np.linalg.svd(X, full_matrices=False)

    # checked this against the sklearn PCA code
    explained_variance = (s ** 2) / n
    explained_variance_ratio = explained_variance / explained_variance.sum()

    preserve_dist_ax = 'samples'
    dims = (0, 1)  # dimensions / components to include in plot
    scale = .05

    # preserve inter-sample distances

    # project gene data into PCA
    # this matches the output of pca.transform() (except for possible sign switch)
    us = u.dot(np.diag(s))

    feat_x = scale * us[:, dims[0]]
    feat_y = scale * us[:, dims[1]]

    sample_x = vh[dims[0]]
    sample_y = vh[dims[1]]

    ## NB if we wish to preserve inter-feature distances, we would need to multiply vh by s instead

    # track legend entries to avoid double labelling
    legend_added = set()

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(feat_x, feat_y, c='gray', s=10, edgecolor='none', alpha=0.7)

    arr_ix, typ = obj.meta.type.factorize()
    for i, t in enumerate(typ):
        for j, pid in enumerate(consts.PIDS):
            if pid in legend_added:
                lbl=None
            else:
                lbl = pid
                legend_added.add(pid)
            marker = scatter_markers[t]
            col = scatter_colours[pid]
            the_ix = (arr_ix == i) & (obj.meta.index.str.contains(pid))
            ax.scatter(
                sample_x[the_ix],
                sample_y[the_ix],
                c=col,
                marker=marker,
                edgecolor='k',
                linewidth=1,
                label=lbl,
                zorder=10
            )

    # connect patients
    for pid in consts.PIDS:
        # for t0, t1 in itertools.combinations(typ):
        ix = obj.meta.index[obj.meta.index.str.contains(pid)]
        for t0, t1 in itertools.combinations(typ, 2):
            # draw all possible connections between these two cell types (in one direction only)
            ix0 = obj.meta.index[obj.meta.index.str.contains(pid) & (obj.meta.type == t0)]
            ix1 = obj.meta.index[obj.meta.index.str.contains(pid) & (obj.meta.type == t1)]

            for a, b in itertools.product(ix0, ix1):
                ax.plot(
                    [sample_x[obj.meta.index == a][0], sample_x[obj.meta.index == b][0]],
                    [sample_y[obj.meta.index == a][0], sample_y[obj.meta.index == b][0]],
                    lw=1.5,
                    color=scatter_colours[pid],
                    zorder=9
                )

    # custom legend outside of plot
    line_kwargs = {
        'class': 'line',
        'markerfacecolor': 'none',
        'markeredgecolor': 'k',
        'markeredgewidth': 1.0,
        'linestyle': 'none'
    }
    patch_kwargs = {
        'class': 'patch',
        'edgecolor': 'k',
        'linewidth': 1.
    }
    legend_dict = {'Patient': collections.OrderedDict(), 'Cell type': collections.OrderedDict()}
    for pid in consts.PIDS:
        ll = dict(patch_kwargs)
        ll['facecolor'] = scatter_colours[pid]
        legend_dict['Patient'][pid] = ll
    for t in typ:
        pp = dict(line_kwargs)
        pp['marker'] = scatter_markers[t]
        legend_dict['Cell type'][t] = pp

    common.add_custom_legend(ax, legend_dict, loc_outside=True)
    ax.set_xlabel('PC%d (%.2f %%)' % (dims[0] + 1, explained_variance_ratio[dims[0]] * 100.))
    ax.set_ylabel('PC%d (%.2f %%)' % (dims[1] + 1, explained_variance_ratio[dims[1]] * 100.))
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d.png" % dims), dpi=200)

    # annotate most influential genes
    selection_radius = 0.6
    selected = (feat_x ** 2 + feat_y ** 2) ** .5 > selection_radius
    genes_selected = dat.index[selected]
    symbols_selected = references.ensembl_to_gene_symbol(genes_selected)


    pass