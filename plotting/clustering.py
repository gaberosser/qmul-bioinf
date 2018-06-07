from microarray import process
from stats import transformations
import numpy as np
import pandas as pd
import seaborn as sns
from plotting.common import COLOUR_BREWERS
from matplotlib import pyplot as plt, colors, patches
from scipy.cluster import hierarchy as hc


def generate_colour_map_dict(
        meta,
        colname,
        contains_arr,
        label=None,
        existing=None,
        sample_names=None,
        non_matching=False,
        group_names=None
):
    """
    Generates a DataFrame colour map suitable for use as a row_color or col_color input to clustermap.
    Also generates a dict suitable for dendrogram plot (used to make the legend)
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
    if len(contains_arr) not in COLOUR_BREWERS:
        raise AttributeError(
            "Unsupported number of categories: %d" % len(contains_arr)
        )
    if non_matching is True and (len(contains_arr) + 1) not in COLOUR_BREWERS:
        raise AttributeError(
            "Unsupported number of categories due to non_matching addition: %d + 1" % len(contains_arr)
        )

    if sample_names is None:
        sample_names = meta.index

    if colname is None:
        attr = meta.index
        if label is None:
            label = attr.name
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

    bg = None
    if non_matching is False:
        cmap = COLOUR_BREWERS[len(contains_arr)]
        col = pd.Series(index=sample_names, name=label)
    else:
        if non_matching is True:
            # auto-select a colour for non matching
            cmap = COLOUR_BREWERS[len(contains_arr) + 1]
            bg = cmap.pop()
        else:
            # use supplied value
            cmap = COLOUR_BREWERS[len(contains_arr)]
            bg = non_matching
        col = pd.Series(data=non_matching, index=sample_names, name=label)

    for i, r in enumerate(contains_arr):
        # NB must use values here because the indexes may not match
        # fillna ensures that there are no NaNs in the array used for lookups
        col.loc[attr.str.contains(r).fillna(False).values] = cmap[i]

    out.loc[:, label] = col

    if group_names is not None:
        legend_dict = {label: dict([
            (group_names[i], cmap[i]) for i in range(len(contains_arr))
        ])}
        if bg is not None:
            legend_dict[label]['Other'] = bg
        return out, legend_dict

    return out


def dendrogram_with_colours(
    data,
    col_colours,
    linkage=None,
    method='average',
    metric='correlation',
    distance_threshold=None,
    fig_kws=None,
    legend_labels=None,
    vertical=True,
    show_labels=True,
):
    """
    IMPORTANT: it is assumed that samples are in COLUMNS. If not, this routine might crash due to memory overflow!
    :param data:
    :param linkage: If supplied, this is the output from hc.linkage (or similar). This overrides method and metric.
    :param method: Passed to hc.linkage()
    :param metric: Passed to hc.linkage()
    :param col_colours: DataFrame with one column per classifying feature. One coloured row will be generated per col.
    :param distance_threshold: If supplied, this is the distance below which the dendrogram is broken into different
    colours
    :param fig_kws: Passed to plt.figure() upon figure creation.
    :param legend_labels: If supplied, a legend is generated explaining the interpretation of the col_colours.
    If col_colours contains just one column, this is a dict with keys corresponding to labels and values corresponding
    to the same colours given in col_colours.
    If multiple columns are found in col_colours, this is a dict with keys matching the column names and values
    dictionaries as above for each separate legend block
    :param vertical: If True (default), the root is at the top and descendants travel downwards. Otherwise, the root
    is at the left and descendants travel right.
    :param show_labels: If True (default), include sample labels in the plot. Need to disable this when the number of
    samples is large.
    :return:
    """
    if distance_threshold is None:
        # set distance threshold so everything is above it
        distance_threshold = -1
    if fig_kws is None:
        fig_kws = dict()
    if linkage is not None:
        z = linkage
    else:
        z = hc.linkage(data.transpose(), metric=metric, method=method)

    nsample = data.shape[1]
    ngroup = col_colours.shape[1]

    fig = plt.figure(**fig_kws)
    if vertical:
        orientation = 'top'
        leg_loc = 'upper left'
        gs = plt.GridSpec(2, 1, height_ratios=(12, 1), hspace=0.01)
        dend_ax = fig.add_subplot(gs[0, 0])
        cc_ax = fig.add_subplot(gs[1, 0])
    else:
        orientation = 'left'
        leg_loc = 'upper right'
        gs = plt.GridSpec(1, 2, width_ratios=(12, 1), wspace=0.01)
        dend_ax = fig.add_subplot(gs[0, 0])
        cc_ax = fig.add_subplot(gs[0, 1])

    if legend_labels is not None:
        if ngroup == 1 and len(legend_labels) != 1:
            legend_labels = {col_colours.columns[0]: legend_labels}

    cc_ax.set_axis_bgcolor(fig.get_facecolor())


    # plot dendrogram
    # labels will be added manually afterwards

    r = hc.dendrogram(
        z,
        ax=dend_ax,
        labels=data.columns,
        above_threshold_color='k',
        color_threshold=distance_threshold,
        no_labels=True,
        orientation=orientation
    )

    if not vertical:
        plt.setp(dend_ax.xaxis.get_ticklabels(), rotation=90)

    # get all colours and convert to unique integer
    unique_cols = np.unique(col_colours.values.flatten())

    # generate matrix and cmap for colour bar
    mat = np.ones((nsample, ngroup))
    cmap_colours = []
    for i, c in enumerate(unique_cols):
        cmap_colours.append(colors.colorConverter.to_rgb(c))
        mat[(col_colours == c).values] = i

    # reshape and reorder based on dendrogram
    if vertical:
        mat = mat.transpose().astype(int)
        mat = mat[:, r['leaves']]
    else:
        mat = mat[r['leaves']][::-1]

    # plot coloured bar using heatmap
    sns.heatmap(mat, cmap=colors.ListedColormap(cmap_colours), ax=cc_ax, cbar=False)

    if vertical:
        xax = cc_ax.xaxis
        yax = cc_ax.yaxis
        xrot = 90
        yrot = 0
    else:
        xax = cc_ax.yaxis
        xax.tick_right()
        yax = cc_ax.xaxis
        xrot = 0
        yrot = 90

    if show_labels:
        xax.set_ticks(np.arange(nsample) + 0.5)
        xax.set_ticklabels(r['ivl'], rotation=xrot)
    yax.set_ticks(np.arange(ngroup) + 0.5)
    yax.set_ticklabels(col_colours.columns, rotation=yrot)

    if legend_labels is not None:
        # draw legend outside of the main axis
        handles = []
        for grp, d in legend_labels.items():
            for ttl, c in d.items():
                handles.append(patches.Patch(color=colors.colorConverter.to_rgb(c), label=ttl))
            handles.append(patches.Patch(color='w', alpha=0., label=''))
        handles.pop()

        if vertical:
            dend_ax.legend(handles=handles, loc=leg_loc, borderaxespad=0., bbox_to_anchor=(1., 1.01))
        else:
            dend_ax.legend(handles=handles, loc=leg_loc, borderaxespad=0., bbox_to_anchor=(-0.02, 1.01))

    gs.tight_layout(fig, h_pad=0., w_pad=0.)

    # this is a HACK: we have to force the figure drawing to get the size of the legend
    fig.canvas.draw()

    # now we need to change right coords to keep legend on the screen
    leg = dend_ax.get_legend()
    if leg is not None:
        bb = leg.get_window_extent()  # bbox in window coordinates
        bb_ax = bb.transformed(leg.axes.transAxes.inverted())  # bbox in ax coordinates
        # import ipdb; ipdb.set_trace()
        if vertical:
            # update the rightmost coordinates
            gs.update(right=(1 - bb_ax.width))
        else:
            # update the leftmost coordinates
            new_left = bb.width / fig.get_window_extent().width
            gs.update(left=new_left + 0.01)

    return {
        'fig': fig,
        'gridspec': gs,
        'dendrogram_ax': dend_ax,
        'col_colour_ax': cc_ax,
        'dendrogram': r,
        'linkage': z,
    }


def plot_correlation_clustermap(data,
                                row_colors=None,
                                n_gene=None,
                                method='average',
                                metric='correlation',
                                distance=None,
                                **kwargs):
    """
    :param n_gene: If supplied, this is the number of genes to use, ordered by descending MAD
    :param kwargs: Passed to seaborn's `clustermap`
    """
    if n_gene is not None:
        # reduce data to the specified number using MAD
        mad = transformations.median_absolute_deviation(data).sort_values(ascending=False)
        genes = mad.index[:n_gene]
        data = data.loc[genes]
    rl = None
    if distance is not None:
        rl = hc.linkage(distance)
        sq = hc.distance.squareform(distance)
    else:
        rl = hc.linkage(data.transpose(), method=method, metric=metric)
        sq = hc.distance.squareform(
            hc.distance.pdist(data.transpose(), metric=metric)
        )

    # invert distance so that closer samples have a larger number
    # do this even if distances have been provided directly
    if metric == 'correlation':
        sq = 1 - sq
    else:
        # TODO: add specific versions for other metrics if required
        sq = max(sq.flat) - sq

    # make a dataframe for clustering so that the plot has correct labels
    sq = pd.DataFrame(data=sq, index=data.columns, columns=data.columns)

    cg = sns.clustermap(
        sq,
        cmap='RdBu_r',
        row_colors=row_colors,
        col_colors=row_colors,
        row_linkage=rl,
        col_linkage=rl,
        **kwargs
    )
    plt.setp(cg.ax_heatmap.get_xticklabels(), rotation=90, fontsize=14)
    plt.setp(cg.ax_heatmap.get_yticklabels(), rotation=0, fontsize=14)
    # shift the margins a bit to fit axis tick labels
    cg.gs.update(bottom=0.2, right=0.8, top=0.99, left=0.01)
    return cg



def plot_clustermap(
        dat,
        method='average', metric='correlation',
        show_gene_labels=False,
        rotate_xticklabels=True,
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
        method=method,
        metric=metric,
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

    if rotate_xticklabels:
        plt.setp(cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90)

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


def dendrogram_threshold_by_nclust(linkage, nclust, sample_subset_idx=None):
    """
    Find the cut point that would split the dendrogram into nclust clusters, optionally
    :param linkage:
    :param nclust:
    :param sample_subset_idx: If supplied, this is a list of indices pointing to the samples that we should use. In this
    case, the distance is computed at the point where _only those leaf nodes_ have formed nclust clusters.
    :return:
    """
    if sample_subset_idx is None:
        # take the mean distance between the nclust-1 and nclust levels
        if nclust == 1:
            dist = np.inf
        elif nclust == 2:
            dist = linkage[-2:, 2].mean()
        else:
            dist = linkage[-nclust:(2 - nclust), 2].mean()
    else:
        ## FIXME: nclust == 2 special case?

        # call the distance based on a fixed number of clusters DEFINED FOR A SUBSET OF NODES
        node_ids = set(sample_subset_idx)
        n = len(linkage) + 1
        cutoff_idx0 = None
        cutoff_idx1 = None
        # loop through linkage rows, keeping track of the specified leaf nodes
        for i, (l0, l1) in enumerate(linkage[:, :2]):
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
        dist = linkage[[cutoff_idx0, cutoff_idx1], 2].mean()
    return dist


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
    dist = dendrogram_threshold_by_nclust(z, nclust, sample_subset_idx=sample_subset_idx)
    show_dendrogram_cut_by_distance(clustermap, dist, axis=axis)

    return dist