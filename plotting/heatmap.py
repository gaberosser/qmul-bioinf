from matplotlib import pyplot as plt, gridspec
import seaborn as sns
import numpy as np


def align_labels(axs, axis='x'):
    """
    :param axs:
    :param axis: 'x' or 'y'
    :return:
    """
    if axis not in ('x', 'y'):
        raise ValueError("Axis must be 'x' or 'y'.")
    if axis == 'x':
        pos = np.array([ax.xaxis.label.get_position() for ax in axs])
        m = pos[:, 1].min()
    else:
        pos = np.array([ax.yaxis.label.get_position() for ax in axs])
        m = pos[:, 0].min()
    # transform into axis units
    inv = axs[0].transAxes.inverted()

    if axis == 'x':
        mt = inv.transform([0, m])[1]
        for ax in axs:
            ax.xaxis.set_label_coords(.5, mt)
    else:
        mt = inv.transform([m, 0])[0]
        for ax in axs:
            ax.yaxis.set_label_coords(mt, .5)


def grouped_expression_heatmap(
        groups,
        data,
        vmax=None,
        cbar=True,
        fig_kwargs=None,
        gs_kwargs=None,
        heatmap_kwargs=None,
        orientation='horizontal'
):
    if orientation not in ('horizontal', 'vertical'):
        raise ValueError("Unsupported orientation %s. Options are horizontal and vertical", orientation)
    horiz =  (orientation == 'horizontal')

    if vmax is None:
        vmin = None
    else:
        vmin = -vmax

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
        heatmap_kwargs.setdefault('cmap', 'RdBu_r')
        heatmap_kwargs.setdefault('cbar_kws', {'orientation': orientation})

    if horiz:
        if cbar:
            gs_kwargs.setdefault('width_ratios', [len(arr) for _, arr in groups])
            gs_kwargs.setdefault('height_ratios', [1, 16])
            gs = gridspec.GridSpec(2, len(groups), **gs_kwargs)
        else:
            gs_kwargs.setdefault('width_ratios', [len(arr) for _, arr in groups])
            gs = gridspec.GridSpec(1, len(groups), **gs_kwargs)

    else:
        if cbar:
            gs_kwargs.setdefault('height_ratios', [len(arr) for _, arr in groups])
            gs_kwargs.setdefault('width_ratios', [16, 1])
            gs = gridspec.GridSpec(len(groups), 2, **gs_kwargs)
        else:
            gs_kwargs.setdefault('height_ratios', [len(arr) for _, arr in groups])
            gs = gridspec.GridSpec(len(groups), 1, **gs_kwargs)

    fig = plt.figure(**fig_kwargs)

    if horiz:
        gs.update(
            left=0.2,
            right=0.95,
            top=0.9,
            bottom=0.1,
            wspace=0.1,
            hspace=0.)
    else:
        gs.update(
            left=0.2,
            right=0.9,
            top=0.98,
            bottom=0.1,
            wspace=0.,
            hspace=0.1)

    axs = []

    for i, (grp, arr) in enumerate(groups):
        if horiz:
            ax = fig.add_subplot(gs[1:, i])
        else:
            ax = fig.add_subplot(gs[i, 0:])
        axs.append(ax)
        cbar_ax = None
        if i == (len(groups) - 1):
            this_cbar = True if cbar else False
            if this_cbar:
                if horiz:
                    cbar_ax = fig.add_subplot(gs[0, :])
                else:
                    cbar_ax = fig.add_subplot(gs[:, 1])
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
        # plt.yticks(rotation=0)
        # plt.xticks(rotation=90)

    # align all xlabels
    fig.canvas.draw()
    if horiz:
        align_labels(axs, 'x')
    else:
        # pass
        align_labels(axs, 'y')

    return fig, axs, cbar_ax
