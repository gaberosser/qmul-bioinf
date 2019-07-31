from matplotlib import pyplot as plt, gridspec
import numpy as np
from scipy import stats
from utils import setops
import collections
import bisect
import multiprocessing as mp


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


def one_random_perm(set_sizes, N):
    rand_sets = [np.random.choice(N, v) for v in set_sizes.values()]
    _, vc = setops.venn_from_arrays(*rand_sets)
    return vc


def set_permutation_test(data, n_iter=1000, parallel=True):
    K = len(data)
    N = len(setops.reduce_union(*data.values()))

    set_sizes = collections.OrderedDict([(k, len(v)) for k, v in data.items()])
    simulated_sizes = collections.defaultdict(list)

    if parallel:
        pool = mp.Pool()
        jobs = {}
        for i in range(n_iter):
            jobs[i] = pool.apply_async(one_random_perm, args=(set_sizes, N))

        pool.close()
        pool.join()
        for i, j in jobs.items():
            vc = j.get()
            for k, v in vc.items():
                simulated_sizes[k].append(v)
    else:
        for i in range(n_iter):
            vc = one_random_perm(set_sizes, N)
            for k, v in vc.items():
                simulated_sizes[k].append(v)

    _, vc_true = setops.venn_from_arrays(*data.values())

    # to calculate the P value, we EITHER need to specify a single sided test OR decide how to compute a two-sided P
    # Some interesting discussions on this topic:
    # https://stats.stackexchange.com/questions/140107/p-value-in-a-two-tail-test-with-asymmetric-null-distribution
    # https://stats.stackexchange.com/questions/360864/2-tailed-permutation-tests-for-obviously-non-symmetric-data
    # https://stats.stackexchange.com/questions/34052/two-sided-permutation-test-vs-two-one-sided
    # However, a 'Z' value is easier to compute
    z = {}
    p = {}
    for k in simulated_sizes.keys():
        obs = vc_true[k]
        t = stats.percentileofscore(simulated_sizes[k], obs)
        if t <= 50:
            p[k] = 2 * t / 100.
        else:
            p[k] = 2 * (1 - t / 100.)

        z[k] = t - 50.

    return {
        'simulated_set_sizes': simulated_sizes,
        'observed_set_sizes': vc_true,
        'p': p,
        'z': z
    }


if __name__ == "__main__":
    """
    Here I'm trying to assemble a function that automates statistical testing of upset plot intersection sizes against
    a fixed-set-size uniform random null.
    """
    n_iter = 1000
    data = {
        'A': range(5) + range(10, 16),
        'B': range(0, 21, 2),
        'C': range(1, 22, 2),
        'D': range(0, 25, 4)
    }

    K = len(data)
    full = setops.reduce_union(*data.values())
    N = len(full)
    set_sizes = collections.OrderedDict([(k, len(v)) for k, v in data.items()])
    # n_intersections = int(sum([special.comb(K, i, exact=True) for i in range(1, K + 1)]))
    intersections = list(setops.binary_combinations(K))
    simulated_sizes = collections.defaultdict(list)

    pool = mp.Pool()
    jobs = {}
    for i in range(n_iter):
        jobs[i] = pool.apply_async(one_random_perm, args=(set_sizes, N))

    pool.close()
    pool.join()
    for i, j in jobs.items():
        vc = j.get()
        for k, v in vc.items():
            simulated_sizes[k].append(v)

    for i in range(n_iter):

        rand_sets = [np.random.choice(N, v) for v in set_sizes.values()]
        _, vc = setops.venn_from_arrays(*rand_sets)
        for k, v in vc.items():
            simulated_sizes[k].append(v)

    _, vc_true = setops.venn_from_arrays(*data.values())

    # to calculate the P value, we EITHER need to specify a single sided test OR decide how to compute a two-sided P
    # Some interesting discussions on this topic:
    # https://stats.stackexchange.com/questions/140107/p-value-in-a-two-tail-test-with-asymmetric-null-distribution
    # https://stats.stackexchange.com/questions/360864/2-tailed-permutation-tests-for-obviously-non-symmetric-data
    # https://stats.stackexchange.com/questions/34052/two-sided-permutation-test-vs-two-one-sided
    # However, a 'Z' value is easier to compute
    z = {}
    p = {}
    for k in simulated_sizes.keys():
        obs = vc_true[k]
        t = stats.percentileofscore(simulated_sizes[k], obs)
        if t <= 50:
            p[k] = 2 * t / 100.
        else:
            p[k] = 2 * (1 - t / 100.)

        z[k] = t - 50.





    ##### Will (hopefully) turn this next part into a function
    set_labels = data.keys()
    dat = data.values()

    set_colours=None
    order_by_n_members=False
    include_singletons=False
    min_size=1
    n_plot=None
    bar_width=0.9
    point_ms=10
    default_colour='#4C72B0'
    kwargs = {}
    n_iter = 1000
    parallel = True

    perm_res = set_permutation_test(dict(zip(set_labels, dat)), n_iter=n_iter, parallel=parallel)

    n_set = len(set_labels)
    if len(data) != len(set_labels):
        raise AttributeError("Number of supplied data sets (%d) doesn't match the length of set_labels (%d)." % (
            len(data), n_set
        ))
    venn_sets, venn_ct = setops.venn_from_arrays(*dat, **kwargs)

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
        ax_intersect.plot(x_arr, np.ones_like(x_arr) * y, marker='o', mfc=lightgrey, mec='none', ms=point_ms,
                          ls='none')
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