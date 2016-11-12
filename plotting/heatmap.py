from matplotlib import pyplot as plt, gridspec
import seaborn as sns
import numpy as np
from plotting import utils


def grouped_expression_heatmap(
        groups,
        data,
        vmax=None,
        vmin=None,
        cbar=True,
        fig_kwargs=None,
        gs_kwargs=None,
        heatmap_kwargs=None,
        orientation='horizontal'
):
    if orientation not in ('horizontal', 'vertical'):
        raise ValueError("Unsupported orientation %s. Options are horizontal and vertical", orientation)
    horiz =  (orientation == 'horizontal')

    # set GridSpec dims inputs based on options
    if horiz:
        ncol = len(groups)
        if cbar:
            nrow = 2
        else:
            nrow = 1
    else:
        nrow = len(groups)
        if cbar:
            ncol = 2
        else:
            ncol = 1

    if fig_kwargs is None:
        fig_kwargs = {}
    if gs_kwargs is None:
        gs_kwargs = {}
    if heatmap_kwargs is None:
        heatmap_kwargs = {}

    heatmap_kwargs.setdefault('vmin', vmin)
    heatmap_kwargs.setdefault('vmax', vmax)
    heatmap_kwargs.setdefault('square', True)
    if cbar:
        default_cbar_kws = {'orientation': orientation}
        if vmax is not None and vmin is not None:
            default_cbar_kws['ticks'] = [vmin, vmax]

        heatmap_kwargs.setdefault('cmap', 'RdBu_r')
        heatmap_kwargs.setdefault('cbar_kws', default_cbar_kws)

    if horiz:
        gs_kw = dict(
            left=0.2,
            right=0.95,
            top=0.9,
            bottom=0.1,
            wspace=0.1,
            hspace=0.05
        )
        gs_kw.update(gs_kwargs)

        if cbar:
            gs_kw.setdefault('width_ratios', [len(arr) for _, arr in groups])
            gs_kw.setdefault('height_ratios', [1, 16])

        else:
            gs_kw.setdefault('width_ratios', [len(arr) for _, arr in groups])

    else:
        gs_kw = dict(
            left=0.2,
            right=0.9,
            top=0.98,
            bottom=0.1,
            wspace=0.05,
            hspace=0.05
        )
        gs_kw.update(gs_kwargs)
        if cbar:
            gs_kw.setdefault('height_ratios', [len(arr) for _, arr in groups])
            gs_kw.setdefault('width_ratios', [16, 1])
        else:
            gs_kw.setdefault('height_ratios', [len(arr) for _, arr in groups])

    gs = gridspec.GridSpec(nrows=nrow, ncols=ncol, **gs_kw)
    fig = plt.figure(**fig_kwargs)

    axs = []

    for i, (grp, arr) in enumerate(groups):
        if horiz:
            if cbar:
                ax = fig.add_subplot(gs[1:, i])
            else:
                ax = fig.add_subplot(gs[0, i])
        else:
            if cbar:
                ax = fig.add_subplot(gs[i, :-1])
            else:
                ax = fig.add_subplot(gs[i, 0])

        axs.append(ax)
        cbar_ax = None
        if i == (len(groups) - 1):
            this_cbar = bool(cbar)
            if this_cbar:
                if horiz:
                    cbar_ax = fig.add_subplot(gs[0, 0])
                else:
                    cbar_ax = fig.add_subplot(gs[0, 1])
        else:
            this_cbar = False

        this_data = data.loc[arr, :]
        if horiz:
            this_data = this_data.transpose()
        sns.heatmap(
            this_data,
            ax=ax,
            cbar=this_cbar,
            cbar_ax=cbar_ax,
            **heatmap_kwargs
        )
        if horiz:
            ax.set_xticklabels(arr, rotation=90)
            ax.set_xlabel(grp)
            if i != 0:
                # only left-most axis needs yticklabels
                ax.set_yticklabels([])

        else:
            ax.set_yticklabels(arr[::-1], rotation=0)  # NB: reverse array ordering due to reversed y axis
            ax.set_ylabel(grp)
            if i != len(groups) - 1:
                # only bottom-most axis needs xticklabels
                ax.set_xticklabels([])

        for tick in ax.get_xticklabels():
            tick.set_rotation(90)
        for tick in ax.get_yticklabels():
            tick.set_rotation(0)

    # align all xlabels
    fig.canvas.draw()
    if horiz:
        utils.align_labels(axs, 'x')
    else:
        # pass
        utils.align_labels(axs, 'y')

    # add axis borders
    for ax in axs:
        utils.axis_border(ax, c='0.3')
    if cbar:
        utils.axis_border(cbar_ax, c='0.3')

    return fig, axs, cbar_ax, gs
