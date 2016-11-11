import numpy as np


def axis_border(ax, c=None, lw=None):
    """
    Add a border around the supplied axis.
    :param c: The line colour. If not supplied, the default is used.
    :param lw: The line width. If not supplied, the default is used.
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