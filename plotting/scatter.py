from matplotlib import pyplot as plt, patches
import seaborn as sns
import pandas as pd
import collections
from plotting import common
import numpy as np


# these markers tend to look smaller (more vertices?) so need to be boosted a little
FILLED_MARKERS_TO_EXPAND = ('P', '*', 'h', 'H', 'p', 'X', '8')


def scatter_with_colour_and_markers(
    dat,
    colour_subgroups=None,
    colour_map=None,
    marker_subgroups=None,
    marker_map=None,
    ax=None,
    legend='outside',
    default_colour='gray',
    default_marker='o',
    ec='k',
    lw=1.0,
    ms=40
):
    """
    :param dat: Data to be plotted in any array format. Expect two columns (x and y). Can also be a pd.DataFrame.
    :param colour_subgroups:
    :param colour_map:
    :param marker_subgroups:
    :param marker_map:
    :param ax:
    :param legend: Include legend in plot? If True, use the 'best' location (according to matplotlib), if 'outside'
    (default), place outside. If False, do not plot legend.
    :param default_colour:
    :param default_marker:
    :param ec: Edgecolour
    :param lw: Linewidth
    :param ms: Marker size

    :return:
    """

    # cast dat to pd DataFrame, should have two columns
    dat = pd.DataFrame(dat)

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')

    c_has_labels = True
    if colour_subgroups is None:
        c_has_labels = False
        colour_subgroups = pd.Series(default_colour, index=dat.index)

    cidx, clabels = colour_subgroups.factorize()

    m_has_labels = True
    if marker_subgroups is None:
        m_has_labels = False
        marker_subgroups = pd.Series(default_marker, index=dat.index)

    midx, mlabels = marker_subgroups.factorize()

    if colour_map is None:
        cmap = common.get_best_cmap(len(clabels))
        colour_map = dict([
            (k, cmap[i]) for i, k in enumerate(clabels)
        ])

    if marker_map is None:
        # marker_map = dict([(k, 'o') for k in mlabels])
        mmap = common.get_best_marker_map(len(mlabels))
        marker_map = dict([
            (k, mmap[i]) for i, k in enumerate(mlabels)
        ])

    for ic, lc in enumerate(clabels):
        for im, lm in enumerate(mlabels):
            c = colour_map.get(lc, default_colour)
            m = marker_map.get(lm, default_marker)

            if m in FILLED_MARKERS_TO_EXPAND:
                # apply a 10% increase to these markers (only)
                this_ms = 1.1 * ms
            else:
                this_ms = ms

            j = (cidx == ic) & (midx == im)
            if j.sum() != 0:
                if c_has_labels and not m_has_labels:
                    lbl = lc
                elif m_has_labels and not c_has_labels:
                    lbl = lm
                else:
                    lbl = None
                ax.scatter(
                    dat.values[j, 0],
                    dat.values[j, 1],
                    c=c,
                    s=this_ms,
                    label=lbl,
                    marker=m,
                    edgecolor=ec,
                    linewidths=lw
                )

    # set legend manually if it requires two groups
    if c_has_labels and m_has_labels:
        for_legend = []

        # colours: show in patches with no edgecolor
        # for lc in clabels:
        for lc in colour_map.keys():
            if lc in clabels:
                the_patch = patches.Patch(
                    edgecolor='none',
                    facecolor=colour_map.get(lc, default_colour),
                    linewidth=lw,
                    label=lc
                )
                for_legend.append(the_patch)

        # spacer that doesn't show up
        the_spacer = patches.Patch(
                edgecolor='none',
                facecolor='none',
                label=''
        )
        for_legend.append(the_spacer)

        # markers: show with no fill
        # for lm in mlabels:
        for lm in marker_map.keys():
            if lm in mlabels:
                the_line = plt.Line2D(
                    [0],
                    [0],
                    marker=marker_map.get(lm, default_marker),
                    markerfacecolor='none',
                    markeredgewidth=lw,
                    markeredgecolor=ec,
                    # markersize=ms,  # the markersize units are different here, so don't specify
                    linestyle='none',
                    linewidth=0.,
                    label=lm
                )
                for_legend.append(the_line)

        if legend == 'outside':
            common.legend_outside_axes(ax, handles=for_legend)
        elif isinstance(legend, str):
            ax.legend(handles=for_legend, loc=legend)
        elif legend:
            ax.legend(handles=for_legend)

    if (legend is not None) and (legend != False):
        if legend == 'outside':
            common.legend_outside_axes(ax)
        elif isinstance(legend, str):
            ax.legend(loc=legend)
        elif legend:
            ax.legend()

    return ax


def qqplot_two_samples(x, y, ax=None, add_loe=True):
    """
    Produce a QQ plot given two arrays of data
    :param x:
    :param y:
    :return:
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    x = sorted(x)
    y = sorted(y)

    # one quartile for every sample in the smaller of the two arrays
    if len(x) >= len(y):
        n_quantile = len(y)
        x_is_idx = False
    else:
        n_quantile = len(x)
        x_is_idx = True

    quantiles = np.linspace(0, 1, n_quantile)

    qq = []
    for i, q in enumerate(quantiles):
        if x_is_idx:
            this = [x[i], y[int(np.round((len(y) - 1) * q))]]
        else:
            this = [x[int(np.round((len(x) - 1) * q))], y[i]]
        qq.append(this)

    qq = np.array(qq)
    ax.scatter(qq[:, 0], qq[:, 1])
    if add_loe:

        zmin = min(min(x), min(y))
        zmax = max(max(x), max(y))
        ax.plot([zmin, zmax], [zmin, zmax], 'k--')

    return qq, ax


def pie_path_marker(weights, start_angle=0, n_vert=64):
    """
    Generate a custom path marker in the form of a pie chart.
    Each segment is a separate marker and needs to be added with a separate call to plt.scatter().
    Weights equal to zero are 'skipped' by generating an empty segment. This means that a fixed length colour array
    can be used when multiple pie charts are required.
    :param weights: Array of weights, which will be normalised for the purpose of calculating segment sizes.
    :param start_angle: The starting angle in radians.
    :param n_vert: The number of vertices to use in the approximation of a circle. The default looks good even at very
    large sizes (s>50000). It may be possible to halve this for smaller markers if many markers are needed.
    :return:
    """
    weights = np.array(weights).astype(float)

    # normalise weights
    weights_n = weights / weights.sum()
    markers = []

    curr_theta = start_angle
    for w in weights_n:
        if w == 0:
            markers.append([[0, 0]])
        else:
            this_n_vert = max(int(w * n_vert), 4)
            end_theta = curr_theta + 2 * np.pi * w
            this_theta = np.linspace(curr_theta, end_theta, this_n_vert)
            x = [0] + np.cos(this_theta).tolist()
            y = [0] + np.sin(this_theta).tolist()
            markers.append(np.column_stack([x, y]))
            curr_theta = end_theta

    return markers


def plot_one_scatter_pie(xy, weights, colours=None, ax=None, marker_kwargs=None, **scatter_kwargs):
    """
    Add one pie chart path marker to the supplied or current axes.
    :param xy: Array of length two specifying the centre of the pie chart.
    :param weights: Array containing unnormalised weights for the pie segment sizes.
    :param colours: Optionally provide an array of colours. If missing, we choose a the 'best colourmap' option. If
    supplied, the length must match that of weights.
    :param ax: Axis object to plot into. If missing, use gca().
    :param marker_kwargs: If supplied, this dictionary is passed to pie_path_markers.
    :param scatter_kwargs: Passed to the scatter() method of ax.
    :return:
    """
    if len(xy) != 2:
        raise ValueError("xy must have length 2.")
    x, y = xy

    if colours is None:
        colours = common.get_best_cmap(len(weights))

    if colours is not None and len(weights) != len(colours):
        raise ValueError("Length of weights and colours must be equal.")

    if marker_kwargs is None:
        marker_kwargs = {}

    if ax is None:
        ax = plt.gca()

    markers = pie_path_marker(weights, **marker_kwargs)
    handles = []
    for m, c in zip(markers, colours):
        h = ax.scatter(
            [x],
            [y],
            marker=m,
            facecolor=c,
            **scatter_kwargs
        )
        handles.append(h)

    return handles


def scatter_with_pies(xy, weights_arr, colours_arr=None, ax=None, marker_kwargs=None, **scatter_kwargs):
    """
    Generate a scatterplot with pie charts as markers.
    :param xy: Array of length N, where N is the number of markers required. Each entry is a length 2 array giving the
    central coordinate of the pie chart.
    :param weights_arr: Array of length N. Each entry is an array of weights.
    :param colours_arr: Optional.
    Either an array of length N with each entry being an array of the same length as the corresponding entry in
    `weights_arr`
    Or a single array of length M, where M is the length of ALL entries of weights_arr (this is checked).
    :param ax: Axis object to plot into. If missing, use gca().
    :param marker_kwargs: If supplied, this dictionary is passed to pie_path_markers.
    :param scatter_kwargs: Passed to the scatter() method of ax.
    :return: Handles generated by the scatter() calls (array of length N, each entry is an array of PathCollection
    objects).
    """

    if len(xy) != len(weights_arr):
        raise ValueError("Length of xy and weights_arr must be equal.")

    if colours_arr is not None:
        if not hasattr(colours_arr[0], '__iter__'):
            # option 1: all weights have the same length and this is the colours array to use for all
            all_len = np.array([len(t) for t in weights_arr])
            if (all_len == len(colours_arr)).all():
                colours_arr = [colours_arr] * len(weights_arr)
            else:
                raise ValueError("If colours_arr is a single array, all weight_arr entries must have the same length.")

        if len(weights_arr) != len(colours_arr):
            raise ValueError("Length of weights_arr and colours_arr must be equal.")

    if ax is None:
        ax = plt.gca()

    handles = []

    for i in range(len(xy)):
        w = weights_arr[i]
        if colours_arr is None:
            colours = common.get_best_cmap(len(w))
        else:
            colours = colours_arr[i]
        if len(colours) != len(w):
            raise ValueError("Pie number %d: number of weights does not match number of colours." % i)

        handles.append(plot_one_scatter_pie(
            xy[i], w, colours=colours, ax=ax, marker_kwargs=marker_kwargs, **scatter_kwargs
        ))

    return handles
