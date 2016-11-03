import pandas as pd
import numpy as np


def plot_correlation_coefficient_array(
        pd_arr,
        cmap='Reds',
        vmin=0,
        vmax=1,
        fmin=None,
        fmax=None,
        fig_kwargs=None,
        plot_kwargs=None,
):
    """
    Plot a right-angled triangular array of coloured squares showing the pairwise correlation between each column in
    the supplied array.
    All columns are used and the column names are applied to the plot.
    """
    from matplotlib import pyplot as plt, rc
    plt.interactive(True)
    import seaborn as sns
    sns.set_style('white')
    # rc('text', usetex=True)

    if fig_kwargs is None:
        fig_kwargs = {}

    if plot_kwargs is None:
        plot_kwargs = {}

    origin = plot_kwargs.pop('origin', 'lower')

    n = pd_arr.shape[1]
    # generate correlation matrix
    c = pd_arr.corr()

    # get min/max if required
    vals = c.values[np.triu_indices_from(c, 1)]
    vals.sort()
    if fmin is not None:
        vmin = vals[int(len(vals) * fmin)]
    if fmax is not None:
        vmax = vals[int(len(vals) * fmax)]

    # zero upper triangle
    c.values[np.triu_indices_from(c, 1)] = np.nan

    # plot
    fig = plt.figure(**fig_kwargs)
    ax = fig.add_subplot(111)
    h = ax.matshow(c, cmap=cmap, vmin=vmin, vmax=vmax, origin=origin, **plot_kwargs)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(pd_arr.columns, rotation=30)
    ax.set_yticklabels(pd_arr.columns)
    fig.colorbar(h)

