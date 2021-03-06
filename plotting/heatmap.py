from matplotlib import pyplot as plt, gridspec
import seaborn as sns
import numpy as np
from plotting import common


def grouped_expression_heatmap(
        groups,
        data,
        vmax=None,
        vmin=None,
        cbar=True,
        fig_kwargs=None,
        gs_kwargs=None,
        heatmap_kwargs=None,
        orientation='horizontal',
        drop_na=True,
):
    """
    Plot a heatmap with subaxes for different groups. Typically, a 'group' is a group of genes.
    :param groups: Iterable of length 2 tuples: (group_name, group_keys). The keys index the ROWS of data.
    :param data: Pandas dataframe. The index (rows) must correspond to keys in groups.
    :param drop_na: If True, any rows that are all NaN will be removed before plotting. Set False if comparing multiple
    plots where some but not all results are NaN (to avoid the number of rows/columns changing).
    """
    if orientation not in ('horizontal', 'vertical'):
        raise ValueError("Unsupported orientation %s. Options are horizontal and vertical", orientation)
    horiz =  (orientation == 'horizontal')

    # we need to define the colour bar boundaries, otherwise the sale could be inconsistent between groups
    if vmin is None:
        vmin = data.values.flatten().min()
    if vmax is None:
        vmax = data.values.flatten().max()

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
    if cbar:
        default_cbar_kws = {'orientation': orientation}
        default_cbar_kws['ticks'] = [vmin, vmax]

        cmap = heatmap_kwargs.setdefault('cmap', 'RdBu_r')
        heatmap_kwargs.setdefault('cbar_kws', default_cbar_kws)

        # create manual scalarmappable for the colour bar
        norm = plt.cm.colors.Normalize(vmin=vmin, vmax=vmax)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(cmap))
        sm.set_array([vmin, vmax])

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
                    cbar_ax = fig.add_subplot(gs[0, :])
                else:
                    cbar_ax = fig.add_subplot(gs[:, 1])
        else:
            this_cbar = False

        this_data = data.loc[arr, :]
        if drop_na:
            this_data = this_data.dropna(axis=0, how='all')
            arr = this_data.index
        if horiz:
            this_data = this_data.transpose()

        sns.heatmap(
            this_data,
            ax=ax,
            cbar=False,
            **heatmap_kwargs
        )
        if horiz:
            ax.set_xticklabels(arr, rotation=90)
            ax.set_xlabel(grp)
            if i != 0:
                # only left-most axis needs yticklabels and label
                ax.set_yticklabels([])
                ax.yaxis.label.set_visible(False)

        else:
            ax.set_yticklabels(arr[::-1], rotation=0)  # NB: reverse array ordering due to reversed y axis
            ax.set_ylabel(grp)
            if i != len(groups) - 1:
                # only bottom-most axis needs xticklabels and label
                ax.set_xticklabels([])
                ax.xaxis.label.set_visible(False)

        for tick in ax.get_xticklabels():
            tick.set_rotation(90)
        for tick in ax.get_yticklabels():
            tick.set_rotation(0)

    if cbar:
        cb = fig.colorbar(sm, cax=cbar_ax, **heatmap_kwargs['cbar_kws'])
        if horiz:
            cbar_ax.xaxis.set_ticks_position('top')
        else:
            cbar_ax.yaxis.set_ticks_position('right')

    # align all xlabels
    fig.canvas.draw()
    if horiz:
        common.align_labels(axs, 'x')
    else:
        # pass
        common.align_labels(axs, 'y')

    # add axis borders
    for ax in axs:
        common.axis_border(ax, c='0.3')
    if cbar:
        common.axis_border(cbar_ax, c='0.3')

    return fig, axs, cbar_ax, gs


def single_heatmap(
        data,
        orientation='horizontal',
        ax=None,
        cbar=True,
        vmin=None,
        vmax=None,
        fig_kwargs=None,
        **heatmap_kwargs
):
    if fig_kwargs is None:
        fig_kwargs = {}
    if ax is None:
        fig = plt.figure(**fig_kwargs)
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    horiz = (orientation == 'horizontal')

    if horiz:
        data = data.transpose()

    sns.heatmap(
        data,
        ax=ax,
        cbar=False,
        vmin=vmin,
        vmax=vmax,
        **heatmap_kwargs
    )

    if cbar:
        quadmesh = ax.collections[-1]  # should only be one item in the array
        cax = fig.colorbar(quadmesh)
        common.axis_border(cax.ax, c='0.3', lw=1.)
    else:
        cax = None

    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    for tick in ax.get_yticklabels():
        tick.set_rotation(0)

    common.axis_border(ax, c='0.3', lw=1.)
    return ax, cax