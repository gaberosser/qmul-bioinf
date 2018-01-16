from matplotlib import pyplot as plt, gridspec
from matplotlib.patches import Ellipse
from matplotlib_venn import venn2, venn3
import numpy as np
from utils import setops


def venn_diagram(*args, **kwargs):
    ax = kwargs.pop('ax', plt.gca())
    n = len(args)
    venn = None
    if n not in {2, 3, 4}:
        raise NotImplementedError("At present, we only support 2, 3 and 4 way Venn diagrams")
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
        ax = fig.add_subplot(111)

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


def upset_set_size_plot(
        data,
        set_labels,
        order_by_n_members=False,
        include_singletons=False,
        min_size=None,
        n_plot=None,
        bar_width=0.9,
        point_ms=10,
        *args,
        **kwargs
):
    """
    Produce a summary plot showing the set sizes when the number of sets is > 4.
    Inspired / totally copying UpsetR: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
    :param data: Array of iterables containing the full data set of each member.
    :param set_labels: Array of strings giving the name of each member, in the same order as data.
    :param order_by_n_members: If True, order the plot by the number of members participating in each set. This has
    the effect of generating a bar chart that has multiple bunches of descending bars.
    :param include_singletons: If True, singleton sets are included in the main bar. Not really necessary as they are
    also plotted in the lower left bar.
    :param min_size: If not None, this is used to exclude sets falling below the minimum size.
    :param n_plot: If not None, this is used to limit the number of sets plotted.
    :param bar_width: Used for plotting bar charts.
    :param point_ms: Size of the circles in the lower right plot.
    """
    venn_sets, venn_ct = setops.venn_from_arrays(*data, **kwargs)

    lightgrey = '#cecece'

    if include_singletons:
        sort_input = venn_ct
    else:
        # exclude any results with only one set
        sort_input = dict([t for t in venn_ct.items() if len(t[0].replace('0', '')) > 1])

    if min_size is not None:
        sort_input = dict([
                              (k, v) for k, v in sort_input.items() if v > min_size
                          ])

    if order_by_n_members:
        ordered_counts = []
        for i in range(1, len(set_labels) + 1):
            # make a collection of results with this many members
            this_collection = []
            for k in setops.binary_combinations_sum_eq(len(set_labels), i):
                # check whether this entry is present, if not it has already been filtered out
                if k in sort_input:
                    this_collection.append((k, sort_input[k]))
            # sort in descending order and append to list
            ordered_counts.extend(sorted(this_collection, key=lambda x: x[1], reverse=True))
    else:
        ordered_counts = sorted(sort_input.items(), key=lambda x: x[1], reverse=True)

    if n_plot:
        ordered_counts = ordered_counts[:n_plot]

    gs_kw = dict(
        left=0.05,
        right=0.99,
        top=0.99,
        bottom=0.1,
        wspace=0.1,
        hspace=0.01,
        height_ratios=[6, 3],
        width_ratios=[3, 6],
    )

    gs = gridspec.GridSpec(nrows=2, ncols=2, **gs_kw)
    fig = plt.figure(figsize=(9, 6))
    ax_tl = fig.add_subplot(gs[0, 0])
    ax_set_size = fig.add_subplot(gs[1, 0])
    ax_intersect = fig.add_subplot(gs[1, 1], sharey=ax_set_size)
    ax_main = fig.add_subplot(gs[0, 1], sharex=ax_intersect)

    # hide some things
    ax_tl.set_visible(False)
    plt.setp(ax_intersect.get_yticklabels(), visible=False)
    plt.setp(ax_main.get_xticklabels(), visible=False)
    plt.setp(ax_intersect.get_xticklabels(), visible=False)

    # data
    x_arr = np.arange(len(ordered_counts)) + 0.5
    y_arr = np.arange(len(set_labels))

    # main bar chart
    ax_main.bar(x_arr, [t[1] for t in ordered_counts], width=bar_width)
    ax_main.set_ylabel('Number of DE genes in set')

    # bottom right set intersections
    # grey markers everywhere
    for y in y_arr:
        ax_intersect.plot(x_arr, np.ones_like(x_arr) * y, marker='o', mfc=lightgrey, mec='none', ms=point_ms, ls='none')
    # black markers only on sets that are included
    for i, (k, v) in enumerate(ordered_counts):
        x = x_arr[i]
        y = [j for j, u in enumerate(k) if u == '1']
        ax_intersect.plot(x * np.ones(len(y)), y, marker='o', mfc='k', mec='k', ms=point_ms, ls='none')

    # bottom left : set size and singleton (unique) set size
    set_sizes = np.array([len(t) for t in data])
    str_fmt = "{0:0%db}" % len(set_labels)
    # NB the string must be reversed here
    singleton_sizes = np.array([venn_ct[str_fmt.format(2 ** i)[::-1]] for i in range(len(set_labels))])

    this_x = set_sizes - singleton_sizes
    ax_set_size.barh(y_arr + 0.5, this_x, -bar_width, align='edge', label='Non-unique')
    ax_set_size.barh(y_arr + 0.5, singleton_sizes, -bar_width, left=this_x, align='edge', color='#ff8484', label='Unique')
    ax_set_size.invert_xaxis()
    ax_set_size.set_ylim([-.5, len(set_labels) - .5])
    ax_set_size.yaxis.tick_right()
    ax_set_size.set_yticks(y_arr)
    ax_set_size.set_yticklabels(set_labels)
    ax_set_size.set_xlabel("Number of DE genes in single comparison")
    ax_set_size.legend(loc='lower left', fontsize=8, frameon=True, fancybox=True, facecolor='w', framealpha=0.3)

    return {
        'gs': gs,
        'axes': {
            'set_size': ax_set_size,
            'intersection': ax_intersect,
            'main': ax_main,
            'top_left': ax_tl
        },
        'figure': fig
    }

