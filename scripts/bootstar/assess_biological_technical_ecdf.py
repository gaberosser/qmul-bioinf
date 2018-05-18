from rnaseq import loader, differential_expression, filter, general
from plotting import common, clustering
from stats import transformations, basic
import pandas as pd
import numpy as np
from scipy import stats
import math
from matplotlib import pyplot as plt
import seaborn as sns
from mpltools import color
from adjustText import adjust_text
from utils import output, setops
import references
import os


def log_cpm(dat, base=2, offset=1.):
    dat = dat + offset
    if len(dat.shape) == 2:
        cpm = dat.divide(dat.sum(), axis=1) * 1e6
    else:
        cpm = dat.divide(dat.sum()) * 1e6
    return np.log(cpm) / np.log(base)


if __name__ == '__main__':
    min_cpm = 1

    outdir = output.unique_output_dir("biological_technical_ecdf")

    # all our patient data

    our_patient_obj = loader.load_by_patient('all', source='star')

    # mouse data

    our_mouse_obj = loader.load_references(['wtchg_p170506', 'wtchg_p170390'], tax_id=10090, strandedness='r')
    # eliminate unneeded samples
    our_mouse_obj.meta = our_mouse_obj.meta.loc[
        ~(our_mouse_obj.meta.index.str.contains('_1') | our_mouse_obj.meta.index.str.contains('CRL3034'))
    ]
    our_mouse_obj.data = our_mouse_obj.data[our_mouse_obj.meta.index.tolist()]

    # SS2 data
    ss2_obj = loader.load_references('wtchg_p180059', strandedness='u')

    all_data = dict([(k, se) for k, se in our_patient_obj.data.iteritems()])
    all_data.update(dict([(k, se) for k, se in our_mouse_obj.data.iteritems()]))

    labels = {}
    for i, k in enumerate(our_patient_obj.meta.index):
        labels[k] = 'Human poly(A)' if i == 0 else None
    for i, k in enumerate(our_mouse_obj.meta.index):
        labels[k] = 'Mouse poly(A)' if i == 0 else None

    plot_style = dict([
        (k, {'c': 'r', 'alpha': 0.4}) for k in our_patient_obj.meta.index
    ])
    plot_style.update(
        dict([
            (k, {'c': 'b', 'alpha': 0.4}) for k in our_mouse_obj.meta.index
        ])
    )

    x_cdf = np.linspace(-5, 15, 500)

    # ECDFs
    # 1: poly(A) mouse and human

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for c, col in all_data.iteritems():
        this_dat = col.loc[col >= min_cpm] + 1
        this_cpm = this_dat.divide(this_dat.sum()) * 1e6
        this_log_cpm = np.log2(this_cpm)
        this_ecdf_fun = basic.ecdf_func(this_log_cpm)
        this_log_cpm_ecdf = this_ecdf_fun(x_cdf)

        ax.plot(x_cdf, this_log_cpm_ecdf, label=labels[c], **plot_style[c])
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")
    # fig.savefig(os.path.join(outdir, "ecdf_by_library_prep.png"), dpi=200)

    # 2: same, just for mouse samples

    # set the colour cycle to distinguish the samples more easily
    color.cycle_cmap(our_mouse_obj.meta.shape[0], cmap='jet')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    e = references.gene_symbol_to_ensembl('Gapdh', tax_id=10090)

    for c, col in our_mouse_obj.data.iteritems():
        this_dat = col.loc[col >= min_cpm] + 1
        this_cpm = this_dat.divide(this_dat.sum()) * 1e6
        this_log_cpm = np.log2(this_cpm)
        this_ecdf_fun = basic.ecdf_func(this_log_cpm)

        this_log_cpm_ecdf = this_ecdf_fun(x_cdf)

        ens_ecdf = this_ecdf_fun(this_log_cpm[e])

        ax.plot(x_cdf, this_log_cpm_ecdf, label=c)
        # optional - plot a single gene on each line
        ax.plot(this_log_cpm[e], ens_ecdf, 'ko')

    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")

    # 3: Add OPCs and NSC (SS2) in
    all_data.update(dict([(k, se) for k, se in ss2_obj.data.iteritems()]))
    opc_lbl = False
    nsc_lbl = False
    for i, k in enumerate(ss2_obj.meta.index):
        if 'NSC' in k:
            plot_style[k] = {'c': 'g', 'alpha': 0.6}
            if not nsc_lbl:
                labels[k] = 'SmartSeq2 NSC'
                nsc_lbl = True
            else:
                labels[k] = None
        if 'OPC' in k:
            plot_style[k] = {'c': 'k', 'alpha': 0.6}
            if not opc_lbl:
                labels[k] = 'SmartSeq2 OPC'
                opc_lbl = True
            else:
                labels[k] = None

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for c, col in all_data.iteritems():
        this_dat = col.loc[col >= min_cpm] + 1
        this_cpm = this_dat.divide(this_dat.sum()) * 1e6
        this_log_cpm = np.log2(this_cpm)
        this_ecdf_fun = basic.ecdf_func(this_log_cpm)
        this_log_cpm_ecdf = this_ecdf_fun(x_cdf)

        ax.plot(x_cdf, this_log_cpm_ecdf, label=labels[c], **plot_style[c])
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")