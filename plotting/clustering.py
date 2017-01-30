from load_data import rnaseq_data
from scripts.comparison_rnaseq_microarray import consts
from microarray import process
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy as hc


def plot_correlation_clustermap(data, row_colors=None, n_gene=None, method='average'):
    """
    :param n_gene: If supplied, this is the number of genes to use, ordered by descending MAD
    """
    if n_gene is not None:
        # reduce data to the specified number using MAD
        mad = process.median_absolute_deviation(data).sort_values(ascending=False)
        genes = mad.index[:n_gene]
        data = data.loc[genes]
    corr = 1. - data.corr()
    z = hc.linkage(corr, method=method)
    cg = sns.clustermap(corr, cmap='RdBu_r', row_colors=row_colors, col_colors=row_colors, row_linkage=z, col_linkage=z)
    plt.setp(cg.ax_heatmap.get_xticklabels(), rotation=90, fontsize=14)
    plt.setp(cg.ax_heatmap.get_yticklabels(), rotation=0, fontsize=14)
    # shift the margins a bit to fit axis tick labels
    cg.gs.update(bottom=0.2, right=0.8, top=0.99, left=0.01)
    return cg
