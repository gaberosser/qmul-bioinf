from matplotlib import pyplot as plt, cm
import seaborn as sns
import numpy as np
import pandas as pd

DEFAULT_COLOURS = ('k', 'b', 'r', 'g', 'y', 'gray', 'c', 'm')


def grouped_bar_chart(series, width=0.8, ax=None, colours=None, ylim=None, mark_clipped=True, labels=None):
    """
    Plot a bar chart with groups.
    :param series: Iterable of data objects. We assume these are ordered correctly.
    :param width: The width of the full group in x-axis space.
    :param ax:
    :param colours: Optional iterable containing the colours for the groups
    :param ylim: Optional iterable of length 2. If supplied, this specifies the values at which y values are clipped.
    :param mark_clipped: If True, any values that are clipped will be filled with a hatched pattern to indicate clipping.
    :return:
    """
    if colours is None:
        colours = DEFAULT_COLOURS
    if ax is None:
        ax = plt.gca()
    ns = len(series)

    dw = width / float(ns)
    max_series_length = 0

    for i, y in enumerate(series):
        # set label if supplied
        lab = None
        if labels is not None:
            lab = labels[i]
        nd = len(y)
        max_series_length = max(max_series_length, nd)
        x = np.arange(nd) + i * dw
        if ylim is None:
            xr = x
            yr = y
        else:
            # clip any extremes
            clip_ind = (y < ylim[0]) | (y > ylim[1])
            # if pandas array, pull out the mask
            if hasattr(clip_ind, 'values'):
                clip_ind = clip_ind.values
            xr = x[~clip_ind]
            yr = y[~clip_ind]
            xc = x[clip_ind]
            yc = y[clip_ind]
            yc[yc < ylim[0]] = ylim[0]
            yc[yc > ylim[1]] = ylim[1]
        col = colours[i % ns]  # cyclic with period n
        # plot unclipped
        ax.bar(xr, yr, dw, facecolor=col, edgecolor='none', label=lab)
        if ylim is not None:
            if mark_clipped:
                hatch = "/"
            else:
                hatch = None
            # plot clipped
            h = ax.bar(xc, yc, dw, facecolor=col, edgecolor='k', hatch=hatch)

    ax.set_xticks(np.arange(max_series_length) + width / 2.)
    # if supplied series are pandas dataframes, use the index to label
    if hasattr(series[0], 'index'):
        ax.set_xticklabels(series[0].index, rotation=90)

    if ylim is not None:
        ax.set_ylim(ylim)


def multi_grouped_bar_chart(dict_of_series_lists, width=0.8, figsize=None, equal_xaxes=True, xlabel_coords=None, **kwargs):
    """
    Produce an array of subplots, all in one row, each of which is a grouped bar chart.
    :param dict_of_series_lists: The keys are used to label the xlabels of the subplots, the values are passed to
    grouped_bar_chart for plotting
    :param width:
    :param figsize: Optional, length 2 iterable with size in inches
    :param equal_xaxes: If True (default) then ensure every subplot xaxis has the same number of xticks
    :param xlabel_coords: If supplied, set the position of the xlabel using these coords. May be necessary to avoid
    clutter beow the axis.
    :param kwargs: Passed to the grouped_bar_chart function
    :return: fig, axs (array of axis handles)
    """
    if figsize is None:
        figsize=(10, 5)

    n_plot = len(dict_of_series_lists)
    max_series_length = max([len(s[0]) for s in dict_of_series_lists.values()])
    fig, axs = plt.subplots(ncols=n_plot, sharey=True, figsize=(10, 5))
    for i, (subplot_lab, arr) in enumerate(dict_of_series_lists.items()):
        ax = axs[i]
        grouped_bar_chart(arr, width=width, ax=ax, **kwargs)
        ax.set_xlabel(subplot_lab)
        if xlabel_coords is not None:
            ax.xaxis.set_label_coords(*xlabel_coords)
        if equal_xaxes:
            ax.set_xlim([width - 1, max_series_length])

    return fig, axs


def stacked_bar_chart(x, y, labels=None, ax=None, colours=None, width=0.9):
    """
    Plot a stacked bar chart. Each group is a row in the matrix y. The number of columns in y is equal to the length
    of x.
    :param x:
    :param y:
    :param labels:
    :param ax:
    :param colours:
    :param width: Relative to the spacing in x (which is assumed constant (TODO?))
    :return:
    """

    n_grp = y.shape[0]
    nx = y.shape[1]
    if len(x) != nx:
        raise ValueError("Length of x must match number of columns in y")
    dx = x[1] - x[0]

    if colours is None:
        colours = [cm.jet(t) for t in np.linspace(0, 1, n_grp)]
    elif len(colours) != n_grp:
        raise ValueError("If supplied, colours must have the same length as the number of rows in y.")

    if labels is not None and len(labels) != n_grp:
        raise ValueError("If supplied, labels must have the same length as the number of rows in y.")

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()

    bottom = np.zeros(nx)

    for i in range(n_grp):
        lbl = labels[i] if labels is not None else None
        ax.bar(x - 0.5 * dx, y[i], bottom=bottom, width=width * dx, color=colours[i], label=lbl, edgecolor='none')
        bottom += y[i]

    if labels is not None:
        ax.legend()

    return fig, ax
