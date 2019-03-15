import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

from scripts.hgic_final import consts
from plotting import common, bar
from settings import HGIC_LOCAL_DIR
from utils import output, setops


def count_de_by_direction(de_res, pids=consts.PIDS, direction_col='logFC'):
    """
    For each patient, count the regions/probes based on direction of differential methylation
    :param dmr_res: Dictionary, keyed by patient ID. Each value is another dictionary, keyed by probe or cluster ID.
    Each value of those contains an object with the attribute `median_change` that is used to assess direction.
    :return: pd.DataFrame, cols correspond to patients, rows correspond to hypo/hyper
    """
    de_direction = {}

    for pid in pids:
        this_es = de_res[pid][direction_col]
        de_direction[pid] = {
            'Up': (this_es > 0).sum(),
            'Down': (this_es < 0).sum(),
        }
    return pd.DataFrame.from_dict(de_direction)[pids].loc[['Down', 'Up']]


def direction_of_de_bar_plot(
        res,
        pids=consts.PIDS,
        stacked=True,
        as_pct=True,
        legend=True,
        legend_loc='outside',
        **kwargs
):
    """
    Generate bar chart showing the split between hypo- and hypermethylation in the supplied differential methylation
    results.
    :param dmr_res: Dictionary, keyed by patient ID. Each value is another dictionary, keyed by probe or cluster ID.
    Each value of those contains an object with the attribute `median_change` that is used to assess direction.
    :param pids: Useful to specify the ordering of the PIDs.
    :param as_pct: If True (default), values are plotted as percentages.
    :param kwargs: Passed to the plotting function `bar.stacked_bar_chart()`. For example, might include `ax=...` to
    specify axis.
    :return:
    """

    for_plot = count_de_by_direction(res, pids=pids)

    if as_pct:
        for_plot = for_plot.divide(for_plot.sum(), axis=1) * 100.
    colours = pd.Series(
        {
            'Up': consts.METHYLATION_DIRECTION_COLOURS['hyper'],
            'Down': consts.METHYLATION_DIRECTION_COLOURS['hypo']
        }
    )

    with sns.axes_style('whitegrid'):
        if stacked:
            fig, ax = bar.stacked_bar_chart(for_plot, colours=colours, width=0.8, legend=legend, **kwargs)
        else:
            ax = for_plot.transpose().plot.bar(
                color=[colours[t] for t in for_plot.index],
                legend=legend,
                **kwargs
            )
            fig = ax.figure
    if as_pct:
        ax.set_ylim([0, 100])
    ax.set_xlim(np.array([0, len(pids)]) - 0.5)
    if legend:
        if legend_loc == 'outside':
            # shrink main axis and put legend on RHS
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
            common.legend_outside_axes(ax)
        else:
            ax.legend(loc=legend_loc)

    return fig, ax


def bar_plot(res, keys=None):

    if keys is None:
        keys = sorted(res.keys())

    # bar plot showing balance of DMR direction
    fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(5.5, 5.5))

    direction_of_de_bar_plot(
        res,
        pids=keys,
        as_pct=True,
        ax=axs[0],
        legend=False
    )
    axs[0].set_ylabel('% DEs')
    direction_of_de_bar_plot(
        res,
        pids=keys,
        as_pct=False,
        ax=axs[1],
        legend=False
    )
    axs[1].set_ylabel('Number DEs')
    plt.setp(axs[1].xaxis.get_ticklabels(), rotation=90)

    return {
        'axs': axs,
        'fig': fig
    }


if __name__ == '__main__':
    """
    Apply the same kind of analysis approach used on the DMR results to the DE results.
    """
    outdir = output.unique_output_dir()
    de_res_fn = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq', 'full_de_syngeneic_only.xlsx')
    pids = consts.PIDS

    de_res_wide = pd.read_excel(de_res_fn)

    # create dictionary from wideform df
    de_res = {}
    for pid in pids:
        this_res = de_res_wide.loc[de_res_wide[pid] == 'Y', de_res_wide.columns.str.contains('%s_' % pid)]
        this_res.columns = this_res.columns.str.replace('%s_' % pid, '')
        this_res.insert(0, 'Gene Symbol', de_res_wide.loc[this_res.index, 'Gene Symbol'])
        de_res[pid] = this_res

    de_by_direction = count_de_by_direction(de_res)

    plt_dict = bar_plot(de_res, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "syngeneic_full_list_directions.png"), dpi=200)

    # specific list
    spec_ix = setops.specific_features(*[de_res[pid].index for pid in pids])
    for_plot = dict([
        (pid, de_res[pid].loc[s]) for pid, s in zip(pids, spec_ix)
    ])
    plt_dict = bar_plot(for_plot, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "syngeneic_specific_list_directions.png"), dpi=200)

    # We don't see the phenotype