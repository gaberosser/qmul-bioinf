from plotting import bar, common, pie
import matplotlib
from matplotlib import pyplot as plt, patches
from matplotlib.colors import Normalize
from matplotlib import cm
import seaborn as sns
import numpy as np
import collections


def polar_heatmap_in_segments(
    density,
    segment_lengths,
    gap_radians=2 * np.pi / 200.,
    inner_r=3,
    delta_r=.25,
    ax=None,
    theta0=0.,
    plot_border=True,
    **kwargs
):
    if set(density.keys()) != set(density.keys()).intersection(segment_lengths.keys()):
        raise ValueError("Each entry of the density input must have an accompanying segment length entry.")

    if ax is None:
        ax = plt.gca(projection='polar')
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location('N')

    sum_segment_length = float(sum([segment_lengths[k] for k in density.keys()]))
    radians_per_unit = (2 * np.pi - gap_radians * len(density)) / sum_segment_length
    outer_r = inner_r + delta_r

    curr_theta = theta0

    # hack required to convert patches correctly to polar coords
    # https://github.com/matplotlib/matplotlib/issues/8521
    ax.bar(0, 1).remove()

    h_pcolor = collections.OrderedDict()
    h_border = collections.OrderedDict()

    for ftr in density.keys():
        th = curr_theta + np.array(density[ftr].index.tolist() + [segment_lengths[ftr]]) * radians_per_unit
        tt = np.array([th, th])
        rr = np.zeros_like(tt) + inner_r
        rr[1] = outer_r

        cc = density[ftr].values[None, :]
        if np.isnan(cc).any():
            cc = np.ma.masked_where(np.isnan(cc), cc)
        h_pcolor[ftr] = ax.pcolor(tt, rr, cc, **kwargs)

        if plot_border:
            this_patch = plt.Rectangle(
                [curr_theta, inner_r],
                width=segment_lengths[ftr] * radians_per_unit,
                height=delta_r,
                edgecolor='k',
                facecolor='none',
                linewidth=1.,
                zorder=999,
            )
            ax.add_artist(this_patch)
            h_border[ftr] = this_patch

        curr_theta = th[-1] + gap_radians

    return {
        'pcolor': h_pcolor,
        'border': h_border
    }

    pass