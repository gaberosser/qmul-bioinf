from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib_venn import venn2, venn3
from utils import setops


def venn_diagram(*args, **kwargs):
    ax = kwargs.pop('ax', plt.gca())
    n = len(args)
    venn = None
    if n not in {2, 3, 4}:
        raise NotImplementedError("At present, we only support 2 and 3 way Venn diagrams")
    venn_sets, venn_counts = setops.venn_from_arrays(*args, **kwargs)
    if n == 2:
        venn = venn2(subsets=venn_counts, ax=ax, **kwargs)
    elif n == 3:
        venn = venn3(subsets=venn_counts, ax=ax, **kwargs)
    elif n == 4:
        venn = venn4(venn_counts, ax=ax, **kwargs)
    return venn, venn_sets, venn_counts


def venn4(data, set_labels=None, show_names=True, ax=None, **kwds):
    from collections import Iterable
    alignment = {'horizontalalignment': 'center', 'verticalalignment': 'baseline'}

    if (set_labels is None) or (len(set_labels) != 4):
        set_labels = ("set 1", "set 2", "set 3", "set 4")

    if ax is None:

        # set figure size
        if 'figsize' in kwds and len(kwds['figsize']) == 2:
            # if 'figsize' is in kwds, and it is a list or tuple with length of 2
            figsize = kwds['figsize']
        else: # default figure size
            figsize = (10, 10)

        fig = plt.figure(figsize=figsize)  # set figure size
        ax = fig.gca()

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['r', 'g', 'b', 'c']

    # draw ellipse, the coordinates are hard coded in the rest of the function

    patches = []
    width, height = 170, 110  # width and height of the ellipses
    patches.append(Ellipse((170, 170), width, height, -45, color=colors[0], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -45, color=colors[1], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -135, color=colors[2], alpha=0.5))
    patches.append(Ellipse((230, 170), width, height, -135, color=colors[3], alpha=0.5))
    for e in patches:
        ax.add_patch(e)
    ax.set_xlim(80, 320); ax.set_ylim(80, 320)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")

    ### draw text
    # 1
    ax.text(120, 200, data['1000'], **alignment)
    ax.text(280, 200, data['0100'], **alignment)
    ax.text(155, 250, data['0010'], **alignment)
    ax.text(245, 250, data['0001'], **alignment)
    # 2
    ax.text(200, 115, data['1100'], **alignment)
    ax.text(140, 225, data['1010'], **alignment)
    ax.text(145, 155, data['1001'], **alignment)
    ax.text(255, 155, data['0110'], **alignment)
    ax.text(260, 225, data['0101'], **alignment)
    ax.text(200, 240, data['0011'], **alignment)
    # 3
    ax.text(235, 205, data['0111'], **alignment)
    ax.text(165, 205, data['1011'], **alignment)
    ax.text(225, 135, data['1101'], **alignment)
    ax.text(175, 135, data['1110'], **alignment)
    # 4
    ax.text(200, 175, data['1111'], **alignment)
    # names of different groups
    if show_names:
        ax.text(110, 110, set_labels[0], **alignment)
        ax.text(290, 110, set_labels[1], **alignment)
        ax.text(130, 275, set_labels[2], **alignment)
        ax.text(270, 275, set_labels[3], **alignment)

    return ax