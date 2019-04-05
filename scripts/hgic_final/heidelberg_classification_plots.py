"""
Added on 4th Sept 2018

Modified version of scripts.publications.cruk_grant_jan_2018.heidelberg_classification_plots, which also loads
raw score data and plots pie charts.

This version just generates a ternary plot for each subtype, showing classification results for each line and
associated FFPE.
"""

from matplotlib import pyplot as plt
import ternary
import numpy as np
import pandas as pd
import collections
import os, sys
from utils import output, setops
from plotting import common
from settings import METHYLATION_ARRAY_DIR
import consts


if __name__ == "__main__":
    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    outdir = output.unique_output_dir(script_name)
    in_fn = os.path.join(METHYLATION_ARRAY_DIR, 'classification.xlsx')

    subgroups = consts.SUBGROUPS
    pids = consts.PIDS
    subgroups_ind = setops.groups_to_ind(pids, subgroups)

    # number of classes to include
    n = 3

    # load summary sheet
    df = pd.read_excel(in_fn, sheet_name='GBM specific', header=0, index_col=None)
    # fix format of PIDs
    df.insert(0, 'pid', df.loc[:, 'Patient ID'].apply(lambda x: "%03.0f" % x))

    # get all classes in top N
    colours = {
        'GBM, MES': '#f4e842',
        'GBM, RTK I': '#0d680f',
        'GBM, RTK II': '#820505',
        'MB, G3': 'orange',
        'PLEX, PED B': '#4C72B0',
        'Other': 'gray'
    }

    # pick out relevant results
    cc_ix = df.pid.isin(consts.PIDS) & (df.Tissue == 'Primary culture')
    ff_ix = df.pid.isin(consts.PIDS) & (df.Tissue == 'FFPE tumour')

    cc = df.loc[cc_ix].set_index('pid').sort_index()
    ff = df.loc[ff_ix].set_index('pid').sort_index()
    cc.insert(1, 'mean_passage', cc.Passage.astype(str).apply(lambda x: np.mean([float(t) for t in x.split(';')])))

    # let's make a TERNARY plot with the three components: RTK I, RTK II, MES
    fig, axs = plt.subplots(ncols=3, figsize=(8, 3))
    taxs = []

    for ax in axs:
        fig, tax = ternary.figure(scale=1., ax=ax)
        taxs.append(tax)
        tax.boundary(linewidth=1.0)
        tax.gridlines(multiple=0.2, color="gray", linewidth=1.0, alpha=0.2)
        tax.clear_matplotlib_ticks()

        ax = tax.get_axes()

        ax.text(-.02, -.02, 'Mesenchymal', horizontalalignment='left', verticalalignment='top')
        ax.text(1.02, -.02, 'RTK I', horizontalalignment='right', verticalalignment='top')
        ax.text(0.5, 0.5 * np.sqrt(3) + 0.02, 'RTK II', horizontalalignment='center', verticalalignment='bottom')

        ax.set_aspect('equal')
        ax.axis('off')

    # bundle them up
    tax_dict = {
        'RTK I': taxs[0],
        'RTK II': taxs[1],
        'MES': taxs[2]
    }

    # each patient is a 'trajectory': FFPE -> early pass -> later passage
    bases = ['GBM_RTK_I', 'GBM_RTK_II', 'GBM_MES']
    # cmap_func = common.continuous_cmap(cc.mean_passage.max(), cmap='Blues')
    cmap = common.get_best_cmap(len(consts.PIDS))
    ff_colour = '#ffffff'

    for p in consts.PIDS:
        this_ff_score = ff.loc[p, bases]
        this_cc = cc.loc[p].sort_values(by='mean_passage')
        this_cc_score = this_cc.loc[:, bases]
        this_cc_pass = this_cc.loc[:, 'mean_passage']
        # this_colours = [ff_colour] + [cmap_func(x - 1) for x in this_cc_pass]
        this_colours = cmap[consts.PIDS.index(p)]
        this_sizes = 20 + this_cc_pass.values ** 2.

        points = np.concatenate([[this_ff_score.values], this_cc_score.values], axis=0)

        ## FIXME: we could just use the hardcoded classification here, eliminating the need to look it up?
        tax = tax_dict[ff.loc[p, 'Result']]
        # here's how to annotate:
        # tax.annotate('GBM%s' % p, points[0], horizontalalignment='left')
        tax.scatter(points[:1], marker='^', c=this_colours, s=30, edgecolor='k', zorder=10, alpha=0.8, label=None)
        tax.scatter(points[1:], marker='o', c=this_colours, s=this_sizes, edgecolor='k', zorder=10, alpha=0.8, label=p)
        for i in range(2):
            tax.line(points[i], points[i+1], c='gray', lw=1., zorder=9)

    # make a custom legend for each subplot
    leg_dicts = {}
    leg_kwargs = {
        'class': 'patch',
        'edgecolor': 'k',
        'linewidth': 1.,
        'alpha': 0.8,
    }
    type_attrs = {
        'class': 'line',
        'linestyle': 'none',
        'markeredgecolor': 'k',
        'markeredgewidth': 1.,
        'markerfacecolor': 'none',
        'markersize': 8
    }

    leg_type_dict = {
        'FFPE': dict(type_attrs),
        'GIC': dict(type_attrs),
    }
    leg_type_dict['FFPE']['marker'] = '^'
    leg_type_dict['GIC']['marker'] = 'o'

    for k in tax_dict:
        this_leg_dict = collections.OrderedDict()
        # FIXME: we use the actual results to determine subgroup earlier, here it's hardcoded
        for pid in subgroups[k]:
            t = dict()
            t.update(leg_kwargs)
            t['facecolor'] = cmap[pids.index(pid)]
            this_leg_dict[pid] = t
        if k == tax_dict.keys()[-1]:
            leg_dicts[k] = collections.OrderedDict()
            leg_dicts[k]['#1'] = this_leg_dict
            leg_dicts[k]['#2'] = leg_type_dict
        else:
            leg_dicts[k] = this_leg_dict
        common.add_custom_legend(
            tax_dict[k].ax,
            leg_dicts[k],
            loc_outside=True,
            loc_outside_horiz='right',
            loc_outside_vert='top',
            bbox_to_anchor=(0.7, 1.),
            frameon=False
        )

    fig.subplots_adjust(left=0.0, right=.95, wspace=0.)
    fig.savefig(os.path.join(outdir, 'ternary_plot_classification.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'ternary_plot_classification.tiff'), dpi=200)
