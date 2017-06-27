import numpy as np

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