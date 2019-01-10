from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from utils import dictionary, string_manipulation, setops
from plotting import common
import collections


def nested_pie_chart(
    dat,
    facecolours,
    legend_entries='same',
    legend_loc='upper right',
    inner_radius=0.,
    width_per_level=0.5,
    startangle=90,
    ax=None,
    **wedgeprops
):
    """
    Plot a nested pie chart. We ensure (through renormalisation) that the plot represents the nesting faithfully.
    :param dat: Dictionary, number of levels is the number of levels of nested segments.
    In each case, the keys should correspond to an entry of the facecolours dictionary.
    :param facecolours: Dictionary keyed according to dat, with values corresponding to colours.
    Any entry with a None value is omitted from the chart (i.e. a gap is left)
    :param legend_entries: If supplied, this is a dictionary containing the labels to add in a legend. If these are
    identical to the keys used, set this to 'same' (default)
    :param legend_loc: Controls legend placement.
    :param edgecolour:
    :param inner_radius: The radius at which the first level commences. Default is starting at the centre (0).
    :param width_per_level: The width of each level of segments.
    :param startangle:
    :param ax:
    :return:
    """
    if ax is None:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')

    col = {}
    col.update(facecolours)

    dat_flat = dictionary.nested_dict_to_flat(dat, ordered=True)
    all_keys = []
    for k in dat_flat:
        for t in k:
            if t not in all_keys:
                all_keys.append(t)
    # all_keys = setops.reduce_union(*dat_flat.keys())

    # link keys to their level
    key_level = {}
    for k in dat_flat:
        for i, t in enumerate(k):
            if t not in key_level:
                key_level[t] = i

    # generate a random 'null key'
    def generate_null_key(all_keys):
        null_key = all_keys[0]
        while null_key in all_keys:
            null_key = string_manipulation.random_string_generator()
        all_keys.append(null_key)
        return null_key

    # generate a null key
    null_key = generate_null_key(all_keys)

    # the null key shouldn't be plotted
    col[null_key] = None

    n_levels = max(len(k) for k in dat_flat)
    # wherever the levels are unequal, pad to the maximum number
    for k in dat_flat:
        if len(k) < n_levels:
            curr = dat_flat.pop(k)
            new_key = k + tuple([null_key] * (n_levels - len(k)))
            dat_flat[new_key] = curr
            # generate a new null key
            null_key = generate_null_key(all_keys)
            col[null_key] = None

    # reform nested dict
    dat_n = dictionary.flat_dict_to_nested(dat_flat, ordered=True)

    def get_level_totals(d):
        res = dict([(k, 0) for k in d])
        for k, v in d.items():
            if isinstance(v, dict):
                sub_res = get_level_totals(v)
                for k2 in v.keys():
                    res[k] += sub_res[k2]
                res.update(sub_res)
            else:
                res[k] += v
        return res


    # get the total for each level, ascending the hierarchy
    level_totals = get_level_totals(dat_n)

    # separate segments by level
    dat_flat_n = dictionary.nested_dict_to_flat(dat_n, ordered=True)

    wedgeprops = dict(wedgeprops)
    wedgeprops['width'] = width_per_level

    # plot
    patches = []
    texts = []
    curr_radius = inner_radius
    parent_keys = None

    for i in range(0, n_levels):
        if parent_keys is None:
            child_keys = dat_n.keys()
        else:
            child_keys = []
            for k1 in parent_keys:
                for k2 in dat_flat_n.keys():
                    if (k2[i - 1] == k1) and (k2[i] not in child_keys):
                        child_keys.append(k2[i])

        this_values = [level_totals[k] for k in child_keys]
        this_cols = [col[k] for k in child_keys]
        ps, ts = ax.pie(
            this_values,
            colors=this_cols,
            radius=curr_radius + width_per_level,
            wedgeprops=wedgeprops,
            startangle=startangle
        )
        curr_radius += width_per_level
        # set segments invisible where required
        this_patches = {}
        this_texts = {}
        for p, c, k, t in zip(ps, this_cols, child_keys, ts):
            if c is None:
                p.set_visible(False)
            this_patches[k] = p
            this_texts[k] = t

        patches.append(this_patches)
        texts.append(this_texts)

        parent_keys = child_keys

    if legend_entries is not None:
        if legend_entries == 'same':
            legend_entries = collections.OrderedDict([(k, k) for k in all_keys])
        legend_dict = collections.OrderedDict()
        for k, txt in legend_entries.items():
            if k in key_level:
                lvl = key_level[k]
                if k not in patches[lvl]:
                    continue
                ## TODO: is there an easier way to generate these kwargs automatically?
                legend_dict[k] = {
                    'class': 'patch',
                    'facecolor': patches[lvl][k].get_facecolor(),
                    'edgecolor': patches[lvl][k].get_edgecolor(),
                    'linewidth': patches[lvl][k].get_linewidth(),
                }
        legend_dict = {'': legend_dict}
        common.add_custom_legend(ax, legend_dict, loc=legend_loc)

    ax.set_xlim(curr_radius * 1.1 * np.array([-1, 1]))

    return {
        'ax': ax,
        'patches': patches,
        'texts': texts
    }


if __name__ == "__main__":
    dat = {'a': {'b': 1, 'c': 3, 'd': {'d1': 4, 'd2': 7}}, 'e': {'f':10, 'g': 20}, 'h': 3}
    facecolours = {
        'a': 'k',
        'b': 'g',
        'c': 'r',
        'd': 'purple',
        'd1': None,
        'd2': 'gray',
        'e': 'm',
        'f': 'y',
        'g': 'b',
        'h': 'w'
    }
    legend_entries = dict([(k, k) for k in facecolours if facecolours[k] is not None])

    res = nested_pie_chart(dat, facecolours, inner_radius=1., edgecolor='k', linewidth=1.)