import pandas as pd
from load_data import microarray_data, rnaseq_data
from microarray import process
from scipy.cluster import hierarchy


if __name__ == '__main__':
    n_genes = 1500
    # load Thompson unannotated dataset
    # Where redundant probes are present, use the one with the highest stdev
    data, _ = microarray_data.load_annotated_thompson2006(aggr_field='SYMBOL', aggr_method='max_std')
    mad = process.median_absolute_deviation(data, axis=1)
    top_genes = mad.sort_values(ascending=False).index[:n_genes]
    X = data.loc[top_genes]

    # unsupervised hierarchical clustering
    z = hierarchy.linkage(X.transpose(), method='average', metric='correlation')

    den = hierarchy.dendrogram(z)