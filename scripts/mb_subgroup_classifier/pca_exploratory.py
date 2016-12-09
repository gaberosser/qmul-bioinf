from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
from load_data import microarray_data, allen_human_brain_atlas, rnaseq_data
from microarray import process
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from scripts.comparison_rnaseq_microarray import consts
from scipy.stats import chi2
import os
import operator
from settings import DATA_DIR


def combine_expr(*args):
    common_genes = args[0].index
    sample_names = args[0].columns
    for t in args[1:]:
        common_genes = common_genes.intersection(t.index)
        sample_names = sample_names.append(t.columns)
    res = pd.DataFrame(index=common_genes, columns=sample_names)
    for t in args:
        res.loc[:, t.columns] = t.loc[common_genes]
    return res


def cluster_ellipsoid_2d(data, p=0.99):
    """
    Compute the parameters required to plot a cluster ellipsoid, based on a multivariate normal assumption.
    :param data: Rows are observations, columns are PCA components
    :param p: The probability cutoff defining the centroid boundary.
    :return: (loc, width, height, angle)
    """
    chi2_cutoff = chi2.ppf(p, 2)  # 2nd arg is DOF

    loc = data.mean(axis=0)
    covar = np.cov(data.transpose())
    eigval, eigvec = np.linalg.eig(covar)
    ang = np.arctan2(eigvec[0, 0], eigvec[1, 0])
    # define radii
    sx = np.sqrt(eigval[0] * chi2_cutoff)
    sy = np.sqrt(eigval[1] * chi2_cutoff)

    return loc, 2 * sx, 2 * sy, 180. * ang / np.pi


if __name__ == '__main__':

    ELL_P = 0.95

    # it's useful to maintain a list of known upregulated genes
    nano_genes = []
    for grp, arr in consts.NANOSTRING_GENES:
        if grp != 'WNT':
            nano_genes.extend(arr)
    nano_genes.remove('EGFL11')
    nano_genes.append('EYS')

    # load Ncott data (285 non-WNT MB samples)
    ncott, ncott_meta = microarray_data.load_annotated_microarray_gse37382(
        aggr_field='SYMBOL',
        aggr_method='max'
    )
    sort_idx = ncott_meta.subgroup.sort_values().index
    ncott_meta = ncott_meta.loc[sort_idx]
    ncott = ncott.loc[:, sort_idx]

    # load Allen (healthy cerebellum)

    he, he_meta = allen_human_brain_atlas.cerebellum_microarray_reference_data(agg_field='gene_symbol', agg_method='max')
    he_meta.loc[:, 'subgroup'] = 'control'

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

    # load Robinson dataset
    robi, robi_meta = microarray_data.load_annotated_microarray_gse37418(aggr_field='SYMBOL', aggr_method='max')
    robi_meta = robi_meta.loc[~robi_meta.subgroup.isin(['U', 'SHH OUTLIER'])]
    sort_idx = robi_meta.subgroup.sort_values().index
    robi_meta = robi_meta.loc[sort_idx]
    robi = robi.loc[:, sort_idx]
    robi_meta.loc[:, 'subgroup'] = robi_meta.subgroup.str.replace('G3', 'Group 3').replace('G4', 'Group 4')

    combs = []
    grps = []
    # combs.append(he); grps.extend(he_meta.subgroup.values)
    # combs.append(ncott); grps.extend(ncott_meta.subgroup.values)
    # combs.append(kool); grps.extend(kool_meta.subgroup.values)
    combs.append(robi); grps.extend(robi_meta.subgroup.values)

    X = combine_expr(*combs)
    m = pd.DataFrame(data=grps, index=X.columns, columns=['subgroup'])

    # PCA fitting requires sample to be represented in the ROWS
    X = process.yugene_transform(X).transpose()
    pca = PCA(n_components=10)
    pca.fit(X)
    y = pca.transform(X)

    # load RNA-Seq data
    xz_data = np.log2(rnaseq_data.gse83696(index_by='Approved Symbol') + 1e-8)
    xz_data[xz_data < 0] = 0.

    X_rna = pd.DataFrame(data=xz_data, columns=xz_data.columns, index=X.columns)
    X_rna.fillna(0, inplace=True)
    X_rna = process.yugene_transform(X_rna).transpose()

    y_rna = pca.transform(X_rna)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(range(1, 11), np.cumsum(pca.explained_variance_), '-o')
    plt.axis([0.8, 10.2, 0, 100.])
    ax.set_xticks(range(1, 11))
    plt.xlabel('Principal component')
    plt.ylabel('Cumulative % variance explained')
    fig.savefig("pca_variance_explained.png", dpi=200)

    # get labels and plot by subgroup
    idx, labels = m.subgroup.factorize()
    colour_map = {
        'WNT': 'c',
        'SHH': 'k',
        'Group 3': 'b',
        'Group 4': 'r',
        'control': 'gray'
    }

    # 2D ellipsoids
    ell = []
    for i, l in enumerate(labels):
        if l not in colour_map:
            continue
        j = idx == i
        ell.append(cluster_ellipsoid_2d(y[j, :2], p=ELL_P))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, l in enumerate(labels):
        if l not in colour_map:
            continue
        c = colour_map[l]
        j = idx == i
        ax.scatter(y[j, 0], y[j, 1], c=c, label=l)
    ax.scatter(y_rna[:, 0], y_rna[:, 1], c='g', label='XZ')
    plt.legend(frameon=False, loc='upper right')
    plt.xlabel("PCA component 1")
    plt.ylabel("PCA component 2")
    fig.savefig("pca_1_2.png", dpi=200)

    for i, l in enumerate(labels):
        if l not in colour_map:
            continue
        c = colour_map[l]
        e = Ellipse(*ell[i], alpha=0.25, facecolor=c, edgecolor=c, lw=2)
        ax.add_artist(e)
    fig.savefig("pca_1_2_ellipses.png", dpi=200)

    # 2D ellipsoids
    ell = []
    for i, l in enumerate(labels):
        if l not in colour_map:
            continue
        j = idx == i
        ell.append(cluster_ellipsoid_2d(y[:, [0, 2]][j], p=ELL_P))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, l in enumerate(labels):
        if l not in colour_map:
            continue
        c = colour_map[l]
        j = idx == i
        ax.scatter(y[j, 0], y[j, 2], c=c, label=l)
    ax.scatter(y_rna[:, 0], y_rna[:, 2], c='g', label='XZ')
    plt.legend(frameon=False, loc='upper right')
    plt.xlabel("PCA component 1")
    plt.ylabel("PCA component 3")
    fig.savefig("pca_1_3.png", dpi=200)

    for i, l in enumerate(labels):
        if l not in colour_map:
            continue
        c = colour_map[l]
        e = Ellipse(*ell[i], alpha=0.25, facecolor=c, edgecolor=c, lw=2)
        ax.add_artist(e)
    fig.savefig("pca_1_3_ellipses.png", dpi=200)

    fig = plt.figure()
    ax_3d = fig.add_subplot(111, projection='3d')
    for i, l in enumerate(labels):
        if l not in colour_map:
            continue
        c = colour_map[l]
        j = idx == i
        ax_3d.scatter(y[j, 0], y[j, 1], y[j, 2], c=c, label=l)
    ax_3d.scatter(y_rna[:, 0], y_rna[:, 1], y_rna[:, 2], c='g', label='XZ')
    plt.legend(frameon=False, loc='upper right')

    ax_3d.set_xlabel("PCA component 1")
    ax_3d.set_ylabel("PCA component 2")
    ax_3d.set_zlabel("PCA component 3")

    fig.savefig("pca_3d.png", dpi=200)
