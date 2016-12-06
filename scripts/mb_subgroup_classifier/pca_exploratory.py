from sklearn.decomposition import PCA
import pandas as pd
from load_data import microarray_data, allen_human_brain_atlas
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from scripts.comparison_rnaseq_microarray import consts
import os
from settings import DATA_DIR


if __name__ == '__main__':

    # it's useful to maintain a list of known upregulated genes
    nano_genes = []
    for grp, arr in consts.NANOSTRING_GENES:
        if grp != 'WNT':
            nano_genes.extend(arr)
    nano_genes.remove('EGFL11')
    nano_genes.append('EYS')

    meta_fn = os.path.join(DATA_DIR, 'microarray_GSE37382', 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=0, sep=',')
    infile = os.path.join(DATA_DIR, 'microarray_GSE37382', 'data.ann.txt.gz')
    res = pd.read_csv(infile, sep='\t', header=0, index_col=0)
    # fix sample names
    res.columns = res.columns.str.replace('.', '_')
    # aggregate by gene symbol, using median to resolve redundancy
    res = res.drop(['ENTREZID', 'GENENAME', 'ACCNUM', 'ENSEMBL'], axis=1)
    res = res.loc[~res.SYMBOL.isnull()]
    # res = res.groupby('SYMBOL', axis=0).median()
    res = res.groupby('SYMBOL', axis=0).max()

    sort_idx = meta.subgroup.sort_values().index
    meta = meta.loc[sort_idx]
    labels = meta.loc[sort_idx, 'subgroup']  # if not sorting, still need to align here
    res = res.loc[:, sort_idx]

    X = res.copy()
    m = meta.copy()

    # PCA fitting requires sample to be represented in the ROWS
    X = res.transpose()
    pca = PCA(n_components=5)
    pca.fit(X)
    y = pca.transform(X)

    # get labels and plot by subgroup
    fig = plt.figure()
    ax = fig.add_subplot(111)
    idx, labels = meta.subgroup.factorize()
    colour_map = {
        'SHH': 'k',
        'Group 3': 'b',
        'Group 4': 'r'
    }
    for i, l in enumerate(labels):
        c = colour_map[l]
        j = idx == i
        ax.scatter(y[j, 0], y[j, 1], c=c, label=l)
    plt.legend(frameon=False, loc='upper right')
    plt.xlabel("PCA component 1")
    plt.ylabel("PCA component 2")
    fig.savefig("pca_all_by_gene.png", dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i, l in enumerate(labels):
        c = colour_map[l]
        j = idx == i
        ax.scatter(y[j, 0], y[j, 1], y[j, 2], c=c, label=l)
    plt.legend(frameon=False, loc='upper right')

    ax.set_xlabel("PCA component 1")
    ax.set_ylabel("PCA component 2")
    ax.set_zlabel("PCA component 3")


    # s = X.std(axis=0).sort_values(ascending=False)
    # k = s.index[:1500]
    # Xs = X.loc[:, k]
    #
    # pca_s = PCA(n_components=5)
    # pca_s.fit(Xs)
    # ys = pca_s.transform(Xs)
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # idx, labels = meta.subgroup.factorize()
    # colour_map = {
    #     'SHH': 'k',
    #     'Group 3': 'b',
    #     'Group 4': 'r'
    # }
    # for i, l in enumerate(labels):
    #     c = colour_map[l]
    #     j = idx == i
    #     ax.scatter(ys[j, 0], ys[j, 1], c=c, label=l)
    # plt.legend(frameon=False, loc='upper right')
    # plt.xlabel("PCA component 1")
    # plt.ylabel("PCA component 2")
    # fig.savefig("pca_top1500_by_gene.png", dpi=200)

