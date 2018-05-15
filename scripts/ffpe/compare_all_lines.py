import numpy as np
import pandas as pd
from rnaseq import loader, filter
from stats import basic, transformations
from matplotlib import pyplot as plt
from plotting import common
import seaborn as sns


if __name__ == "__main__":
    min_cpm = 1
    obj = loader.load_by_patient('all', type='ffpe')

    pids = loader.PATIENT_LOOKUP_FFPE.keys()
    batch_order = sorted(obj.meta.batch.unique())

    colour_by_patient = dict([
        (p, common.COLOUR_BREWERS[len(pids)][i]) for i, p in enumerate(pids)
    ])

    # filter
    dat = filter.filter_by_cpm(obj.data, min_n_samples=1, min_cpm=min_cpm)

    # ECDF
    x_cdf = np.linspace(-5, 15, 500)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    log_cpm_ecdf = {}

    for s, col in dat.iteritems():
        this_dat = col.loc[col >= min_cpm] + 1
        this_cpm = this_dat.divide(this_dat.sum()) * 1e6
        this_log_cpm = np.log2(this_cpm)
        this_ecdf_fun = basic.ecdf_func(this_log_cpm)
        this_log_cpm_ecdf = this_ecdf_fun(x_cdf)
        log_cpm_ecdf[s] = this_log_cpm_ecdf
        p = obj.meta.loc[s, 'reference_id'].replace('GBM', '')
        c = colour_by_patient[p]
        ls = '-'

        this_batch = obj.meta.loc[s, 'batch']
        match_ix = (obj.meta.reference_id == obj.meta.loc[s, 'reference_id']) & (obj.meta.index != s)

        if match_ix.sum() > 0:
            other_batch_inds = [batch_order.index(t) for t in obj.meta.loc[match_ix, 'batch']]
            this_batch_ix = batch_order.index(this_batch)
            if this_batch_ix > max(other_batch_inds):
                ls = '--'

        ax.plot(x_cdf, this_log_cpm_ecdf, ls=ls, c=c, label='%s (%s)' % (s, this_batch))
    ax.legend(loc='upper left')


    # take logCPM values
    logcpm = np.log10((dat + 1).divide(dat.sum() + 1, axis=1) * 1e6)
