from microarray import process
import numpy as np
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



def plot_clustermap(
        dat,
        show_gene_labels=False,
        **kwargs
):
    """
    Plot a heatmap of (e.g.) gene expression data, showing clustering of samples and optionally genes.
    :param dat: pandas DataFrame of the gene expression data, with samples on the columns
    :param show_gene_labels: If True, the genes are labelled. This won't fit if the number shown is too large.
    :param kwargs: Passed directly to Seaborn's `clustermap` function.  Examples of some that I've used:
        cmap: string or ColorMap instance
        col_colors=pd.DataFrame, indexes should match the columns in `dat`, columns represent different attributes
        that we want to show. The elements are colour strings (hex). Each column results in a new bar at the top of the
        plot.
        col_linkage: Specify our own linkage array, which is computed using `scipy.cluster.hierarchy.linkage`
        (or similar). This controls the dendrogram at the top of the plot.
        z_score: 0 means compute Z score by row, 1 by column, default is neither.
        xticklabels: boolean
    :return: Instance of Clustermap. This is useful for other functions (below).
    """
    cg = sns.clustermap(
        dat,
        **kwargs
    )
    # check whether x ticks were requested - if so, rotate and raise the bottom of the axes
    if kwargs.get('xticklabels', False):
        plt.setp(cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90)
        bottom = 0.1

    else:
        bottom = 0.02
    # remove useless row dendrogram
    cg.ax_row_dendrogram.set_visible(False)
    # and remove the space created for it
    wr = cg.gs.get_width_ratios()
    wr[0] = 0.035
    wr[1] = 0.02
    cg.gs.set_width_ratios(wr)
    # reduce whitespace
    cg.gs.update(bottom=bottom, top=0.98, left=0.02)

    cg.ax_heatmap.yaxis.label.set_visible(False)
    cg.ax_heatmap.xaxis.label.set_visible(False)
    if show_gene_labels:
        plt.setp(
            cg.ax_heatmap.yaxis.get_ticklabels(),
            rotation=0,
            fontsize=14
        )
    else:
        cg.ax_heatmap.yaxis.set_ticklabels([])

    return cg


def show_dendrogram_cut_by_distance(clustermap, dist, axis=1):
    """
    Add a dashed line to the specified dendrogram showing the cutoff point.
    The line is specified by distance
    :param clustermap: ClusterMap instance from the seaborn `clustermap` function
    :param dist: Distance, in the units used when plotting the clustermap
    :param axis: 0 (row) or 1 (column)
    :return: None
    """
    if axis not in (0, 1):
        raise ValueError("Valid choices for `axis` are (0, 1).")

    if axis == 0:
        dend = clustermap.dendrogram_row
        ax = clustermap.ax_row_dendrogram
    else:
        dend = clustermap.dendrogram_col
        ax = clustermap.ax_col_dendrogram

    xmin = np.array(dend.independent_coord).min()
    xmax = np.array(dend.independent_coord).max()
    rng = xmax - xmin

    plot_kwargs = dict(
        linestyle='--',
        color='k',
        alpha=0.3
    )

    if axis == 0:
        ax.plot(
            [dist, dist],
            [xmin - 0.1 * rng, xmax + 0.1 * rng],
            **plot_kwargs
        )
    else:
        ax.plot(
            [xmin - 0.1 * rng, xmax + 0.1 * rng],
            [dist, dist],
            **plot_kwargs
        )


def show_dendrogram_cut_by_nclust(clustermap, nclust, axis=1, sample_subset_idx=None):
    """
    Add a dashed line to the specified dendrogram showing the cutoff point.
    The line is specified by the desired final number of clusters
    :param clustermap: ClusterMap instance from the seaborn `clustermap` function
    :param nclust: The final number of clusters required. The distance is computed using this.
    :param axis: 0 (row) or 1 (column)
    :param sample_subset_idx: If supplied, this is a list of indices pointing to the samples that we should use. In this
    case, the distance is computed at the point where _only those leaf nodes_ have formed nclust clusters.
    :return: None
    """
    if axis == 0:
        dend = clustermap.dendrogram_row
    else:
        dend = clustermap.dendrogram_col
    z = dend.linkage

    if sample_subset_idx is None:
        # take the mean distance between the nclust-1 and nclust levels
        dist = z[-nclust:(2 - nclust), 2].mean()
    else:
        # call the distance based on a fixed number of clusters DEFINED FOR A SUBSET OF NODES
        node_ids = set(sample_subset_idx)
        n = len(z) + 1
        cutoff_idx0 = None
        cutoff_idx1 = None
        # loop through linkage rows, keeping track of the specified leaf nodes
        for i, (l0, l1) in enumerate(z[:, :2]):
            l0 = int(l0)
            l1 = int(l1)
            if l0 in node_ids:
                node_ids.remove(l0)
                node_ids.add(n + i)
            if l1 in node_ids:
                node_ids.remove(l1)
                # fine to do this twice since we are working with a set
                node_ids.add(n + i)
            if len(node_ids) == nclust:
                cutoff_idx0 = i
            if len(node_ids) == (nclust - 1):
                cutoff_idx1 = i
                break
        if cutoff_idx0 is None or cutoff_idx1 is None:
            raise ValueError("Failed to compute the requested cluster distance")
        dist = z[[cutoff_idx0, cutoff_idx1], 2].mean()

    show_dendrogram_cut_by_distance(clustermap, dist, axis=axis)

    return dist