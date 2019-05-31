import numpy as np
from stats import basic
from matplotlib import pyplot as plt
import seaborn as sns


def log_cpm_ecdf_plot(
        data_dict,
        units='counts',
        log_cpm_lookup_values=None,
        label_dict=None,
        style_dict=None,
        min_cpm=1,
        ax=None,
        legend=True
):
    """
    Generate a plot showing the ECDF of log2(CPM) values for each of the elements of data_dict.
    :param data_dict: Dictionary of pd.Series objects, each containing RAW COUNTS for one sample.
    :param units: If 'counts' (default), we start by calculating CPM values. If 'cpm', we skip this step.
    :param label_dict:
    :param style_dict:
    :param min_cpm:
    :param ax:
    :return:
    """
    lbl_dict = {
        'counts': 'log2(CPM)',
        'cpm': 'log2(CPM)',
        'tpm': 'log2(TPM)',
    }
    units = units.lower()

    if units not in lbl_dict.keys():
        raise NotImplementedError("Unrecognised units '%s'" % units)

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    transformed_dat = {}
    for c, col in data_dict.iteritems():

        this_dat = col.loc[col > 0]
        if units == 'counts':
            this_dat = this_dat / this_dat.sum() * 1e6

        if min_cpm is not None:
            this_dat = this_dat.loc[this_dat >= min_cpm]

        transformed_dat[c] = np.log2(this_dat)

    if log_cpm_lookup_values is None:
        log_cpm_lookup_values = np.linspace(
            np.floor(np.log2(min_cpm) if min_cpm > 0 else -3),
            np.ceil(max([dat.max() for dat in transformed_dat.values()])),
            100
        )

    for c, col in transformed_dat.iteritems():

        if label_dict is not None:
            lbl = label_dict.get(c, c)
        else:
            lbl = c

        sd = {}
        if style_dict is not None:
            sd = style_dict[c]

        this_ecdf_fun = basic.ecdf_func(col)
        this_log_cpm_ecdf = this_ecdf_fun(log_cpm_lookup_values)
        ax.plot(log_cpm_lookup_values, this_log_cpm_ecdf, label=lbl, **sd)

    if legend:
        ax.legend(loc='lower right')

    ax.set_xlabel(lbl_dict[units])
    ax.set_ylabel("Empirical CDF")

    return ax