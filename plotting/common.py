import numpy as np
from stats import basic
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm


COLOUR_BREWERS = {
    2: ['#1f78b4', '#b2df8a'],
    3: ['#7fc97f', '#beaed4', '#fdc086'],
    4: ['#7fc97f', '#beaed4', '#fdc086', '#ffff99'],
    5: ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0'],
    6: ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f'],
    7: ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f', '#bf5b17'],
    8: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00'],
    9: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6'],
    10: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'],
    11: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a',
         '#ffff99'],
    12: ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a',
         '#ffff99', '#b15928'],
}

_FILLED_MARKERS = ['o', 'v', 's', 'd', 'p', 'X', '*', 'D', 'P', '^', '<', '>', 'h', 'H', '8']

FILLED_MARKERS = dict([(i, _FILLED_MARKERS[:i]) for i in range(2, len(_FILLED_MARKERS))])


def legend_outside_axes(ax, horiz_loc='right', vert_loc='centre', **kwargs):
    if horiz_loc == 'right':
        x = 1.
        xa = 'left'
    elif horiz_loc == 'left':
        x = 0.
        xa = 'right'
    elif horiz_loc == 'centre':
        x = 0.5
        xa = 'center'
    else:
        raise NotImplementedError("Unsupported horiz_loc %s" % horiz_loc)

    if vert_loc == 'top':
        y = 1.
        ya = 'upper'
    elif vert_loc == 'bottom':
        y = 0.
        ya = 'lower'
    elif vert_loc == 'centre':
        y = 0.5
        ya = 'center'
    else:
        raise NotImplementedError("Unsupported vert_loc %s" % vert_loc)

    loc_str = ' '.join([ya, xa])
    if loc_str == 'center center':
        loc_str = 'center'

    ax.legend(loc=loc_str, bbox_to_anchor=(x, y), **kwargs)


def get_best_cmap(N, cmap='jet'):
    """
    Get the best colourmap we can create for the given number of categories.
    We preferentially use colour brewer, but if the number of categories exceeds 12 then we revert to a continuous
    cmap
    :param N: Number of categories
    :param cmap: The colourmap to use if we need to work with a continuous scale (not used if N<=12).
    :return: List of colours, either RGB or hex.
    """
    if N in COLOUR_BREWERS:
        return COLOUR_BREWERS[N]
    else:
        func = continuous_cmap(N, cmap=cmap)
        return [func(i) for i in range(N)]


def continuous_cmap(N, cmap='jet'):
    '''
    Returns a function that maps each index in 0, 1, ... N-1 to a distinct
    RGB color, defined by cmap.
    '''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap=cmap)
    def map_index_to_rgb_color(index):
        return colors.to_hex(scalar_map.to_rgba(index))
    return map_index_to_rgb_color


def axis_border(ax, c=None, lw=1):
    """
    Add a border around the supplied axis.
    :param c: The line colour. If not supplied, the default is used.
    :param lw: The line width.
    """
    ax.set_frame_on(True)
    for s in ax.spines.values():
        s.set_visible(True)
        if c is not None:
            s.set_color(c)
        if lw is not None:
            s.set_linewidth(lw)


def align_labels(axs, axis='x'):
    """
    Align the labels on the specified axis.
    This is useful for subplots, where the labels can shift around depending on the size of the ticklabels.
    The algorithm proceeds by finding the most distant label and moving all others to match it.
    :param axs:
    :param axis: 'x' or 'y'
    :return:
    """
    if axis not in ('x', 'y'):
        raise ValueError("Axis must be 'x' or 'y'.")
    if axis == 'x':
        pos = np.array([ax.xaxis.label.get_position() for ax in axs])
        m = pos[:, 1].min()
    else:
        pos = np.array([ax.yaxis.label.get_position() for ax in axs])
        m = pos[:, 0].min()
    # transform into axis units
    inv = axs[0].transAxes.inverted()

    if axis == 'x':
        mt = inv.transform([0, m])[1]
        for ax in axs:
            ax.xaxis.set_label_coords(.5, mt)
    else:
        mt = inv.transform([m, 0])[0]
        for ax in axs:
            ax.yaxis.set_label_coords(mt, .5)


def ecdf_plot(
        data_dict,
        label_dict=None,
        style_dict=None,
        xi=None,
        ax=None,
        xmin=None,
        xmax=None,
        legend=True
):
    """
    Generate a plot showing the ECDF for each of the elements of data_dict.
    :param data: Iterable of pd.Series or np.array objects, each containing data for one sample.
    :param label_dict: If supplied, it is a dictionary containing alternative labels to those used in data_dict.
    Keys must match those in data_dict.
    :param style_dict:
    :param xi: Values to use for the x axis of the ECDF. These are automatically sele
    :param xmin: If specified, data are discarded below this value
    :param ax:
    :return:
    """

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if xmin is None:
        xmin = 1e12
        for v in data_dict.values():
            xmin = min(xmin, min(v))

    if xmax is None:
        xmax = -1e12
        for v in data_dict.values():
            xmax = max(xmax, max(v))

    if xi is None:
        xi = np.linspace(xmin, xmax, 100)

    for k, v in data_dict.items():
        lbl = None
        sd = {}
        lbl = k
        if label_dict is not None:
            lbl = label_dict[k]

        if style_dict is not None:
            sd = style_dict[k]

        this_dat = v.loc[v >= xmin]
        this_ecdf_fun = basic.ecdf_func(this_dat)

        yi = this_ecdf_fun(xi)
        ax.plot(xi, yi, label=lbl, **sd)

    if legend:
        ax.legend(loc='lower right')

    return ax