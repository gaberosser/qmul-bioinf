from sklearn.decomposition import PCA
import pandas as pd
from load_data import microarray_data, allen_human_brain_atlas
from matplotlib import pyplot as plt
import seaborn as sns


if __name__ == '__main__':
    res, meta = microarray_data.load_annotated_microarray_gse37382(index_field='gene_symbol')
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

    s = X.std(axis=0).sort_values(ascending=False)
    k = s.index[:1500]
    Xs = X.loc[:, k]

    pca_s = PCA(n_components=5)
    pca_s.fit(Xs)
    ys = pca_s.transform(Xs)

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
        ax.scatter(ys[j, 0], ys[j, 1], c=c, label=l)
    plt.legend(frameon=False, loc='upper right')
    plt.xlabel("PCA component 1")
    plt.ylabel("PCA component 2")
    fig.savefig("pca_top1500_by_gene.png", dpi=200)

