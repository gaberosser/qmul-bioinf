from microarray import process
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy as hc


def generate_colour_map_dict(meta, colname, contains_arr, label=None, existing=None, sample_names=None, non_matching=False):
    """
    Generates a DataFrame colour map suitable for use as a row_color or col_color input to clustermap.
    Colours are selected automatically based on the number of categories, up to a total of 9.
    If an existing dataframe is supplied, an additional column is added for the new label
    :param meta: The full metadata.
    :param colname: The name of the column to access in metadata. If None, index is used.
    :param contains_arr: An array of strings or compiled regex objects. These will be used to match on `contains`.
    :param label: If supplied, this label is used as the column name in the output. Otherwise, colname is used.
    :param existing: If supplied, this is an existing DataFrame. A new column is added and a copy returned.
    :param sample_names: If supplied, this array contains the samples names, used as the index of the output. Otherwise
    use the columns of metadata.
    :param non_matching: If True, also apportion a colour for cases that don't match anything in contains_arr. If a
    string, use this as the extra colour. If False, ignore non-matching colours.
    :return:
    """
    colour_brewers = {
        2: ['#1f78b4','#b2df8a'],
        3: ['#7fc97f','#beaed4','#fdc086'],
        4: ['#7fc97f','#beaed4','#fdc086','#ffff99'],
        5: ['#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0'],
        6: ['#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f'],
        7: ['#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17'],
        8: ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00'],
        9: ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6'],
        10: ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'],
        11: ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99'],
        12: ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'],
    }
    if len(contains_arr) not in colour_brewers:
        raise AttributeError(
            "Unsupported number of categories: %d" % len(contains_arr)
        )
    if non_matching is True and (len(contains_arr) + 1) not in colour_brewers:
        raise AttributeError(
            "Unsupported number of categories due to non_matching addition: %d + 1" % len(contains_arr)
        )

    if sample_names is None:
        sample_names = meta.index

    if colname is None:
        attr = meta.index
        if label is None:
            label = meta.name
    else:
        attr = meta.loc[:, colname]
        if label is None:
            label = colname

    if existing is not None:
        out = existing.copy()
        if len(out.index.intersection(sample_names)) != len(out.index):
            raise AttributeError("The index of the supplied existing dataframe do not match the expected sample names")
    else:
        out = pd.DataFrame(index=sample_names)

    if non_matching is False:
        cmap = colour_brewers[len(contains_arr)]
        col = pd.Series(index=sample_names, name=label)
    elif non_matching is True:
        cmap = colour_brewers[len(contains_arr) + 1]
        bg = cmap.pop()
        col = pd.Series(data=bg, index=sample_names, name=label)
    else:
        cmap = colour_brewers[len(contains_arr)]
        col = pd.Series(data=non_matching, index=sample_names, name=label)

    for i, r in enumerate(contains_arr):
        # NB must use values here because the indexes may not match
        col.loc[attr.str.contains(r).values] = cmap[i]

    out.loc[:, label] = col

    return out


def dendrogram_with_colours(
    data,
    col_colours,
    method='average',
    metric='correlation'
):
    """
    IMPORTANT: it is assumed that samples are in COLUMNS. If not, this routine might crash due to memory overflow!
    :param data:
    :param method:
    :param metric:
    :param col_colours:
    :return:
    """
    gs = plt.GridSpec(2, 1, height_ratios=(12, 1), hspace=0.01)
    fig = plt.figure()
    axs = [fig.add_subplot(ax) for ax in gs]
    axs[1].yaxis.set_visible(False)
    axs[1].set_facecolor(fig.get_facecolor())

    # plot dendrogram
    z = hc.linkage(data.transpose(), metric=metric, method=method)
    r = hc.dendrogram(z, ax=axs[0])

    # plot colours
    matrix = np.ones((1, data.shape[1]))
    
    sns.matrix.ClusterGrid.color_list_to_matrix_and_cmap(list(col_colours.group.values), range(len(col_colours)))




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