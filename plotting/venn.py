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
    set_colours=None,
    order_by_n_members=False,
    include_singletons=False,
    min_size=1,
    n_plot=None,
    bar_width=0.9,
    point_ms=10,
    default_colour='#4C72B0',
    **kwargs
):
    """
    Produce a summary plot showing the set sizes when the number of sets is > 4.
    Inspired / totally copying UpsetR: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
    :param data: Array of iterables containing the full data set of each member.
    :param set_labels: Array of strings giving the name of each member, in the same order as data.
    :param set_colours: Dict/list of tuples/OrderedDict giving the name and shading of one or more sets. E.g.
    [
        (group_A, {'sets': ['010', '011'], 'colour': 'red'}),
        (group_B, {'sets': ['110', '001'], 'colour': 'blue'}),
    ]
    The name is used for the legend. It can be 'None' to skip any entry for this group.
    Use the ordered options if order matters in the lower left stacked plot. Use a list to have multiple entries
    with the same group name.
    If supplied, these will be used for shading all three plots. If not, we just shade the singleton sets in the
    lower left set size plot.
    :param order_by_n_members: If True, order the plot by the number of members participating in each set. This has
    the effect of generating a bar chart that has multiple bunches of descending bars.
    :param include_singletons: If True, singleton sets are included in the main bar. Not really necessary as they are
    also plotted in the lower left bar.
    :param min_size: This is used to exclude sets falling below the minimum size. Can be disabled (set to None), but
    this is pointless since it involves plotting empty sets, which cannot be ordered meaningfully.
    :param n_plot: If not None, this is used to limit the number of sets plotted.
    :param bar_width: Used for plotting bar charts.
    :param point_ms: Size of the circles in the lower right plot.
    :param default_colour: The colour used for anything that isn't otherwise shaded.
    """

    n_set = len(set_labels)
    if len(data) != len(set_labels):
        raise AttributeError("Number of supplied data sets (%d) doesn't match the length of set_labels (%d)." % (
            len(data), n_set
        ))
    venn_sets, venn_ct = setops.venn_from_arrays(*data, **kwargs)

    if set_colours is None:
        str_fmt = "{0:0%db}" % n_set
        # NB the string must be reversed here
        singleton_sets = set([str_fmt.format(2 ** i)[::-1] for i in range(n_set)])
        other_sets = set([k for k in venn_ct if k not in singleton_sets])
        set_colours = [
            ('Non-unique', {'sets': other_sets, 'colour': default_colour}),
            ('Unique', {'sets': singleton_sets, 'colour': '#ff8484'}),
        ]
    else:
        try:
            set_colours = set_colours.items()
        except AttributeError:
            set_colours = list(set_colours)
        sets_seen = set()
        for nm, d in set_colours:
            this_sets = d['sets']
            if len(sets_seen.intersection(this_sets)):
                raise ValueError(
                    "Group %s contains one or more sets already contained elsewhere in set_colours" % nm
                )
            sets_seen.update(this_sets)
        sets_remaining = set(venn_ct.keys()).difference(sets_seen)
        if len(sets_remaining) > 0:
            set_colours = [(None, {'sets': sets_remaining, 'colour': default_colour})] + set_colours

    # convenience function to find the colour matching a given set
    def set_lookup(k):
        for t in set_colours:
            if k in t[1]['sets']:
                return t

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

    # set up axis grid
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
    y_arr = np.arange(n_set)

    # main bar chart
    colours = [set_lookup(t[0])[1]['colour'] for t in ordered_counts]
    ax_main.bar(x_arr, [t[1] for t in ordered_counts], width=bar_width, color=colours)
    ax_main.set_ylabel('Number of DE genes in set')

    # bottom right set intersections
    # grey markers everywhere
    for y in y_arr:
        ax_intersect.plot(x_arr, np.ones_like(x_arr) * y, marker='o', mfc=lightgrey, mec='none', ms=point_ms, ls='none')
    # overplot shaded markers on sets that are included
    for i, (k, v) in enumerate(ordered_counts):
        x = x_arr[i]
        y = [j for j, u in enumerate(k) if u == '1']
        c = set_lookup(k)[1]['colour']
        ax_intersect.plot(x * np.ones(len(y)), y, marker='o', mfc=c, mec=c, ms=point_ms, ls='none')

    # bottom left : set size and singleton (unique) set size
    left = np.zeros(n_set)

    set_sizes = []
    for nm, d in set_colours:
        this_ss = np.zeros(n_set)
        for k in d['sets']:
            for i in range(n_set):
                if k[i] == '1':
                    this_ss[i] += venn_ct[k]
        set_sizes.append([nm, this_ss])
        ax_set_size.barh(
            y_arr + 0.5,
            this_ss,
            height=-bar_width,
            left=left,
            align='edge',
            label=nm,
            color=d['colour']
        )
        left += this_ss

    ax_set_size.invert_xaxis()
    ax_set_size.set_ylim([-.5, len(set_labels) - .5])
    ax_set_size.yaxis.tick_right()
    ax_set_size.set_yticks(y_arr)
    ax_set_size.set_yticklabels(set_labels)
    ax_set_size.set_xlabel("Number of DE genes in single comparison")
    ax_set_size.legend(
        loc='lower left',
        # fontsize=8,
        frameon=False,
        facecolor='w',
        # edgecolor='k',
        bbox_to_anchor=(0.05, 1.1),  # place above and outside the axis
    )

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


def expanded_core_sets(venn_set, subgroup_ind):
    """
    Compute the sets that belong to the 'expanded core', which comprises any combination of members that bridges
    multiple subgroups.
    :param venn_set: As returned by setops.venn_from_arrays
    :param subgroup_ind: Dict, keys are the names of the subgroups, values are boolean arrays with the membership for
    each patient. It's important that the ordering here is the same as venn_set.
    :return:
    """
    ecs = []
    for k in venn_set:
        this_k = np.array([t for t in k]).astype(bool)
        nmatch = 0
        for grp, grp_idx in subgroup_ind.items():
            if this_k[grp_idx].any():
                nmatch += 1
        if nmatch > 1:
            # add to the expanded core set
            ecs.append(k)
    return ecs


def upset_plot_with_groups(
        data,
        set_labels,
        subgroup_ind,
        subgroup_colours,
        venn_set=None,
        other_lbl='Expanded core',
        specific_lbl='Specific',
        default_colour='gray',
        **kwargs
):
    """
    Wrapper around the basic upset plotting function. This allows us to highlight sets that fully or partially
    overlap with a pre-defined subgroup.
    :param data: Passed to upset_set_size_plot. Iterable of identifiers used to process venn sets.
    :param set_labels: Iterable of set labels.
    :param subgroup_ind: Dictionary, keys are set_labels, entries are Boolean indexes showing which of set_labels
    are in this subgroup. If ordering is desired, use an OrderedDict.
    :param subgroup_colours: Dict giving the colour for each of the subsets defined in subgroup ind. For each set S,
    two entries are needed, keyed `S full` and `S partial`.
    We can also define two additional colours, which otherwise have default values:
    `Expanded core` (or whatever `other_lbl` is set to) and `Specific`.
    :param venn_set: Output of setops.venn_from_arrays(data). Can supply it to skip recomputing.
    :param other_lbl: Label used to identify those sets that span multiple subgroups.
    :param specific_lbl: Label used to identify those sets that are specific to a single member.
    :param kwargs: Passed to upset_set_size_plot
    :return: Same output as upset plot function.
    """
    # UpsetR attribute plots
    default_colour_other = '#4C72B0'
    default_colour_specific = '#f4e842'

    if venn_set is None:
        venn_set, _ = setops.venn_from_arrays(*data)

    # set colours for UpsetR plot
    sets_full = {}
    sets_partial = {}
    sets_unique = []

    ## TODO: merge this with setops.full_partial_unique_other_sets_from_groups
    for k in venn_set:
        this_k = np.array([t for t in k]).astype(bool)
        if this_k.sum() == 1:
            sets_unique.append(k)
        elif this_k.sum() > 1:
            for grp, grp_idx in subgroup_ind.items():
                n_member = grp_idx.sum()
                # no other matches
                if this_k[~grp_idx].sum() == 0:
                    if this_k[grp_idx].sum() == n_member:
                        sets_full.setdefault(grp, []).append(k)
                    else:
                        sets_partial.setdefault(grp, []).append(k)

    set_colours = []
    for grp_name in subgroup_ind:
        k_full = "%s full" % grp_name
        if grp_name in sets_full:
            set_colours.append(
                (k_full, {'sets': sets_full[grp_name], 'colour': subgroup_colours.get(k_full, default_colour)})
            )

        k_part = "%s partial" % grp_name
        if grp_name in sets_partial:
            set_colours.append(
                (k_part, {'sets': sets_partial[grp_name], 'colour': subgroup_colours.get(k_part, default_colour)})
            )

    set_colours.append(
        (other_lbl, {
            'sets': expanded_core_sets(venn_set, subgroup_ind),
            'colour': subgroup_colours.get(other_lbl, default_colour_other)
        }),
    )

    set_colours.append(
        (specific_lbl, {
            'sets': sets_unique,
            'colour': subgroup_colours.get(specific_lbl, default_colour_specific)
        }),
    )

    return upset_set_size_plot(
        data,
        set_labels,
        set_colours=set_colours,
        default_colour=default_colour,
        **kwargs
    )