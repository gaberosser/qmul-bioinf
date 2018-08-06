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
        marker_map = dict([(k, 'o') for k in mlabels])

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
        for lc in clabels:
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
        for lm in mlabels:
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

