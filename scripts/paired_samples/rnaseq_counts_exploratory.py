from load_data import rnaseq_data
from scripts.comparison_rnaseq_microarray import consts
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy as hc


if __name__ == '__main__':
    tpm = rnaseq_data.gbm_paired_samples(annotate_by=None, units='tpm')
    # dendrogram / hierarchical clustering with correlation distance
    corr = 1. - tpm.corr()
    z = hc.linkage(corr, method='average')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    dendro = hc.dendrogram(z, labels=corr.columns, color_threshold=.15, leaf_rotation=45, ax=ax)

    # set some colours for easy visual clarification of the groupings (black and grey)
    cols = pd.Series(data=['#000000'] * 5 + ['#BDBDBD'] * 5, index=corr.index, name='')

    cg = sns.clustermap(corr, cmap='RdBu_r', row_colors=cols, col_colors=cols, row_linkage=z, col_linkage=z)
    plt.setp(cg.ax_heatmap.get_xticklabels(), rotation=90, fontsize=14)
    plt.setp(cg.ax_heatmap.get_yticklabels(), rotation=0, fontsize=14)
    # shift the margins a bit to fit axis tick labels
    cg.gs.update(bottom=0.2, right=0.8, top=0.99, left=0.01)
    cg.savefig("clustermap_gbm_rnaseq_pairs.pdf")
    cg.savefig("clustermap_gbm_rnaseq_pairs.png", dpi=300)