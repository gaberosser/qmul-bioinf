# from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist
import pandas as pd
from load_data import microarray_data, allen_human_brain_atlas, rnaseq_data
from microarray import process
from matplotlib import pyplot as plt, gridspec
from matplotlib.colors import ListedColormap
import seaborn as sns
from settings import DATA_DIR
import numpy as np
import os
from scripts.comparison_rnaseq_microarray import consts


def dendro_heatmap(Z, labels, show_xlabels=True, cmap=None, dendro_kwargs=None, heatmap_kwargs=None):
    """
    :param labels: pd.Series, indexed by sample names, with values giving group names
    """
    if cmap is None:
        cmap = "Greys"
    if dendro_kwargs is None:
        dendro_kwargs = {}
    if heatmap_kwargs is None:
        heatmap_kwargs = {}

    gs_kw = dict(
        left=0.1,
        right=0.9,
        top=0.98,
        bottom=0.1,
        wspace=0.05,
        hspace=0.05,
        height_ratios=[24, 1]
    )
    gs = gridspec.GridSpec(nrows=2, ncols=1, **gs_kw)
    fig = plt.figure()
    ax = fig.add_subplot(gs[0])
    subax = fig.add_subplot(gs[1])

    den = dendrogram(Z, leaf_rotation=90, ax=ax, no_labels=True, labels=labels.index, **dendro_kwargs)

    # reorder the labels according to the dendrogram sorting
    labels = labels.loc[den['ivl']]
    label_ind, label_str = labels.factorize(sort=True)  # sorting ensures the same label order each time

    xticklabels = []
    if show_xlabels:
        xticklabels = labels.index.values

    sns.heatmap([label_ind], ax=subax, cbar=False, cmap=cmap, xticklabels=xticklabels, **heatmap_kwargs)

    for tick in subax.get_xticklabels():
        tick.set_rotation(90)

    return fig, ax, subax, gs


if __name__ == '__main__':
    # TODO: fix data loading functions

    ncott, ncott_meta = microarray_data.load_annotated_microarray_gse37382(
        aggr_field='SYMBOL',
        aggr_method='max'
    )
    sort_idx = ncott_meta.subgroup.sort_values().index
    ncott_meta = ncott_meta.loc[sort_idx]
    ncott = ncott.loc[:, sort_idx]

    # X = ncott.copy()
    # m = ncott_meta.copy()


    # load Kool dataset
    kool, kool_meta = microarray_data.load_annotated_microarray_gse10327(
        aggr_field='SYMBOL',
        aggr_method='max',
    )
    sort_idx = kool_meta.subgroup.sort_values().index
    kool_meta = kool_meta.loc[sort_idx]
    kool = kool.loc[:, sort_idx]
    kool_meta.loc[:, 'subgroup'] = (
        kool_meta.loc[:, 'subgroup'].str
            .replace('A', 'WNT')
            .replace('B', 'SHH')
            .replace('E', 'Group 3')
            .replace('C', 'Group 4')
            .replace('D', 'Group 4')
    )

    # X = kool.copy()
    # m = kool_meta.copy()

    # load Robinson dataset
    robi, robi_meta = microarray_data.load_annotated_microarray_gse37418(aggr_field='SYMBOL', aggr_method='max')
    robi_meta = robi_meta.loc[~robi_meta.subgroup.isin(['U', 'SHH OUTLIER'])]
    sort_idx = robi_meta.subgroup.sort_values().index
    robi_meta = robi_meta.loc[sort_idx]
    robi = robi.loc[:, sort_idx]
    robi_meta.loc[:, 'subgroup'] = robi_meta.subgroup.str.replace('G3', 'Group 3').replace('G4', 'Group 4')

    X = robi.copy()
    m = robi_meta.copy()



    # YuGene
    X = process.yugene_transform(X)

    # the samples must be represented in ROWS
    X = X.transpose()

    Z = linkage(X, method='average', metric='correlation')
    # fig, ax, subax, gs = dendro_heatmap(Z, m.subgroup)
    # check the cophenetic distance: the sorrelation between ACTUAL pairwise distances between samples and the
    # distances according to the hierarchical clustering
    # c, coph_dists = cophenet(Z, pdist(X))
    # print "Correlation coeff between actual pdist and hierarchical pdist is %.2f" % c

    # pick top high stdev genes
    s = X.std(axis=0).sort_values(ascending=False)
    idx = s.index[:1500]

    Xs = X.loc[:, idx]
    Zs = linkage(Xs, method='average', metric='correlation')
    fig, ax, subax, gs = dendro_heatmap(Zs, m.subgroup, show_xlabels=False)
    ax.set_ylabel("Correlation distance")

    c_s, coph_dists_s = cophenet(Zs, pdist(Xs))
    print "Correlation coeff between actual pdist and hierarchical pdist is %.2f" % c_s

    nano_genes = []
    for grp, arr in consts.NANOSTRING_GENES:
        # skip WNT
        if grp == 'WNT':
            continue
        nano_genes.extend(arr)
    nano_genes.remove('EGFL11')
    nano_genes.append('EYS')

    # Xn = X.loc[:, nano_genes]
    # Zn = linkage(Xn, method='average', metric='correlation')

    # fig, ax, subax, gs = dendro_heatmap(Zn, m.subgroup, show_xlabels=False)

    # c_s, coph_dists_s = cophenet(Zn, pdist(Xn))
    # print "Correlation coeff between actual pdist and hierarchical pdist is %.2f" % c_s
