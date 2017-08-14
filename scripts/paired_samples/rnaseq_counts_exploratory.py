from load_data import rnaseq_data
from scripts.comparison_rnaseq_microarray import consts
from stats import transformations
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
        mad = transformations.median_absolute_deviation(data).sort_values(ascending=False)
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


if __name__ == '__main__':
    tpm, meta = rnaseq_data.gbm_paired_samples(annotate_by=None, units='tpm')
    # dendrogram / hierarchical clustering with correlation distance
    corr = 1. - tpm.corr()
    z = hc.linkage(corr, method='average')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    dendro = hc.dendrogram(z, labels=corr.columns, color_threshold=.15, leaf_rotation=45, ax=ax)

    # set some colours for easy visual clarification of the groupings (black and grey)
    cols = pd.Series(data=['#000000'] * 5 + ['#BDBDBD'] * 5, index=corr.index, name='')
    cg = plot_correlation_clustermap(tpm, row_colors=cols)
    cg.savefig("clustermap_gbm_rnaseq_pairs.pdf")
    cg.savefig("clustermap_gbm_rnaseq_pairs.png", dpi=300)

    # repeat but only keep top 1500 most variable genes
    n_gene = 1500
    cg = plot_correlation_clustermap(tpm, row_colors=cols, n_gene=n_gene)
    cg.savefig("clustermap_gbm_rnaseq_pairs.top%d.pdf" % n_gene)
    cg.savefig("clustermap_gbm_rnaseq_pairs.top%d.png" % n_gene, dpi=300)

    # repeat but remove GBM24 and DURA24 to get greater sensitivity among the others
    tpm_no24 = tpm.drop(['GBM024', 'DURA024N28_NSC'], axis=1)
    cols = pd.Series(index=tpm_no24.columns, name='')
    cols.loc[tpm_no24.columns.str.contains('GBM')] = '#000000'
    cols.loc[tpm_no24.columns.str.contains('DURA')] = '#BDBDBD'
    cg = plot_correlation_clustermap(tpm_no24, row_colors=cols, n_gene=n_gene)
    cg.savefig("clustermap_gbm_rnaseq_pairs.remove24.top%d.pdf" % n_gene)
    cg.savefig("clustermap_gbm_rnaseq_pairs.remove24.top%d.png" % n_gene, dpi=300)