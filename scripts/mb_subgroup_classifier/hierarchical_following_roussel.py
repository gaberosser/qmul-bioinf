import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from load_data import microarray_data, rnaseq_data
from microarray import process
from scipy.cluster import hierarchy
from scripts.comparison_rnaseq_microarray import consts


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

    # plot: key gene expression plus clustering
    g = [
        'WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2',
        'PDLIM3', 'EYA1', 'ATOH1', 'SFRP1',
        'IMPG2', 'GABRA5', 'NRL', 'MAB21L2', 'NPR3',
        'KCNA1', 'KHDRBS2', 'OAS1'
    ]
    expr = data.loc[g]
    row_colors = pd.Series(['#2D438E'] * 5 + ['#E5161C'] * 4 + ['#F2EA00'] * 5 + ['#2A8C43'] * 3, index=g, name='')

    cg = sns.clustermap(
        expr,
        cmap='RdBu_r',
        row_colors=row_colors,
        row_cluster=None,
        col_linkage=z,
        z_score=0,
        xticklabels=False,
    )
    cg.ax_heatmap.yaxis.label.set_visible(False)
    plt.setp(
        cg.ax_heatmap.yaxis.get_ticklabels(),
        rotation=0,
        fontsize=14
    )

    # now we assume that 4 clusters have been found and use those for colours
    lbl = hierarchy.fcluster(z, 4, criterion='maxclust')
    col_colors = pd.Series(index=expr.columns, name='')
    col_colors.loc[lbl == 1] = '#000000'
    col_colors.loc[lbl == 2] = '#a5a5a5'
    col_colors.loc[lbl == 3] = '#6d6d6d'
    col_colors.loc[lbl == 4] = '#dbdbdb'

    cg = sns.clustermap(
        expr,
        cmap='RdBu_r',
        row_colors=row_colors,
        col_colors=col_colors,
        row_cluster=None,
        col_linkage=z,
        z_score=0,
        xticklabels=False,
    )
    cg.ax_heatmap.yaxis.label.set_visible(False)
    plt.setp(
        cg.ax_heatmap.yaxis.get_ticklabels(),
        rotation=0,
        fontsize=14
    )