# from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist
import pandas as pd
from load_data import microarray_data, allen_human_brain_atlas
from matplotlib import pyplot as plt
import seaborn as sns



if __name__ == '__main__':
    res, meta = microarray_data.load_annotated_microarray_gse37382(index_field='gene_symbol')
    # clustering requires sample to be represented in the ROWS
    X = res.transpose()

    # pick top high stdev genes
    s = X.std(axis=0).sort_values(ascending=False)
    idx = s.index[:1500]


    Z = linkage(X, method='average', metric='correlation')

    # check the cophenetic distance: the sorrelation between ACTUAL pairwise distances between samples and the
    # distances according to the hierarchical clustering
    c, coph_dists = cophenet(Z, pdist(X))
    print "Correlation coeff between actual pdist and hierarchical pdist is %.2f" % c
    plt.figure()
    dendrogram(Z, leaf_rotation=90)

    Xs = X.loc[:, idx]
    Zs = linkage(Xs, method='average', metric='correlation')
    plt.figure()
    dendrogram(Zs, leaf_rotation=90)

    c_s, coph_dists_s = cophenet(Zs, pdist(Xs))
    print "Correlation coeff between actual pdist and hierarchical pdist is %.2f" % c_s