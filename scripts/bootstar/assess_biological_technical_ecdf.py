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
import collections
from plotting.rnaseq import log_cpm_ecdf_plot
import os


def log_cpm(dat, base=2, offset=1.):
    dat = dat + offset
    if len(dat.shape) == 2:
        cpm = dat.divide(dat.sum(), axis=1) * 1e6
    else:
        cpm = dat.divide(dat.sum()) * 1e6
    return np.log(cpm) / np.log(base)


if __name__ == '__main__':
    min_cpm = 0.01

    outdir = output.unique_output_dir("biological_technical_ecdf")

    # all our patient data (cell culture)

    our_patient_obj = loader.load_by_patient('all', source='star')

    # all our patient data (FFPE culture)
    ffpe_samples = [
        'NH15_1661DEF2C',
        'NH15_1877_SP1C',
        'NH15_2101_DEF1A',
        'NH16_270_DEF1Ereplacement',
        'NH16_616DEF1B',
        'NH16_677_SP1A',
        'NH16_1574DEF1A',
        'NH16_1976_DEF1Areplacement',
        'NH16_2063_DEF1Areplacement',
        'NH16_2214DEF1A',
        'NH16_2255DEF1B2',
        'NH16_2806DEF3A1',
    ]
    our_ffpe_obj = loader.load_by_patient('all', type='ffpe')
    our_ffpe_obj.meta = our_ffpe_obj.meta.loc[ffpe_samples]
    our_ffpe_obj.data = our_ffpe_obj.data[our_ffpe_obj.meta.index]

    # mouse data

    our_mouse_obj = loader.load_references(['wtchg_p170506', 'wtchg_p170390'], tax_id=10090, strandedness='r')
    # eliminate unneeded samples
    our_mouse_obj.meta = our_mouse_obj.meta.loc[
        ~(our_mouse_obj.meta.index.str.contains('_1') | our_mouse_obj.meta.index.str.contains('CRL3034'))
    ]
    our_mouse_obj.data = our_mouse_obj.data[our_mouse_obj.meta.index.tolist()]

    # SS2 data
    ss2_obj = loader.load_references('wtchg_p180059', strandedness='u')
    ss2_obj.meta.index = ["%s_SS2" % t for t in ss2_obj.meta.index]
    ss2_obj.data.columns = ss2_obj.meta.index


    # define counts seprately for convenience
    polya_human_counts = our_patient_obj.data
    polya_mouse_counts = our_mouse_obj.data
    ffpe_counts = our_ffpe_obj.data
    ss2_counts = ss2_obj.data

    first_label = {}
    for i, k in enumerate(our_patient_obj.meta.index):
        first_label[k] = 'Human poly(A)' if i == 0 else None
    for i, k in enumerate(our_mouse_obj.meta.index):
        first_label[k] = 'Mouse poly(A)' if i == 0 else None
    for i, k in enumerate(our_ffpe_obj.meta.index):
        first_label[k] = 'Human FFPE Ribozero' if i == 0 else None
    opc_lbl = False
    nsc_lbl = False
    for i, k in enumerate(ss2_obj.meta.index):
        if 'NSC' in k:
            if not nsc_lbl:
                first_label[k] = 'SmartSeq2 NSC'
                nsc_lbl = True
            else:
                first_label[k] = None
        if 'OPC' in k:
            if not opc_lbl:
                first_label[k] = 'SmartSeq2 OPC'
                opc_lbl = True
            else:
                first_label[k] = None


    plot_style = dict([
        (k, {'c': 'r', 'alpha': 0.6}) for k in our_patient_obj.meta.index
    ])
    plot_style.update(
        dict([
            (k, {'c': 'b', 'alpha': 0.6}) for k in our_mouse_obj.meta.index
        ])
    )
    plot_style.update(
        dict([
            (k, {'c': 'y', 'alpha': 0.7}) for k in our_ffpe_obj.meta.index
        ])
    )
    for i, k in enumerate(ss2_obj.meta.index):
        if 'NSC' in k:
            plot_style[k] = {'c': 'g', 'alpha': 0.6}
        if 'OPC' in k:
            plot_style[k] = {'c': 'k', 'alpha': 0.6}

    x_cdf = np.linspace(-6, 12.5, 500)

    # ECDFs
    # -1: FFPE, poly(A) human, poly(A) mouse
    this_data = dict(our_patient_obj.data.iteritems())
    this_data.update(our_ffpe_obj.data.iteritems())
    this_data.update(our_mouse_obj.data.iteritems())

    ax = log_cpm_ecdf_plot(
        this_data,
        label_dict=first_label,
        log_cpm_lookup_values=x_cdf,
        style_dict=plot_style,
        min_cpm=min_cpm
    )

    # again with norming
    polya_human_cpm_n = transformations.edger_tmm_normalisation_cpm(our_patient_obj.data)
    polya_mouse_cpm_n = transformations.edger_tmm_normalisation_cpm(our_mouse_obj.data)
    ffpe_cpm_n = transformations.edger_tmm_normalisation_cpm(our_ffpe_obj.data)
    ss2_cpm_n = transformations.edger_tmm_normalisation_cpm(ss2_obj.data)

    this_data = dict(polya_human_cpm_n.iteritems())
    this_data.update(polya_mouse_cpm_n.iteritems())
    this_data.update(ffpe_cpm_n.iteritems())
    this_data.update(ss2_cpm_n.iteritems())

    ax = log_cpm_ecdf_plot(
        this_data,
        units='cpm',
        label_dict=first_label,
        log_cpm_lookup_values=x_cdf,
        style_dict=plot_style,
        min_cpm=min_cpm
    )

    # 1: poly(A) mouse and human
    this_data = dict(polya_human_counts.iteritems())
    this_data.update(dict(polya_mouse_counts.iteritems()))

    ax = log_cpm_ecdf_plot(
        this_data,
        log_cpm_lookup_values=x_cdf,
        label_dict=first_label,
        style_dict=plot_style,
        min_cpm=min_cpm
    )
    # fig.savefig(os.path.join(outdir, "ecdf_by_library_prep.png"), dpi=200)

    # 2: same, just for mouse samples

    # set the colour cycle to distinguish the samples more easily
    color.cycle_cmap(our_mouse_obj.meta.shape[0], cmap='jet')

    this_data = collections.OrderedDict(sorted(polya_mouse_counts.iteritems()))
    ax = log_cpm_ecdf_plot(
        this_data,
        log_cpm_lookup_values=x_cdf,
        label_dict=dict([(k, k) for k in this_data]),
        # style_dict=plot_style,
        min_cpm=min_cpm,
    )

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    #
    # e = references.gene_symbol_to_ensembl('Gapdh', tax_id=10090)
    #
    # for c, col in our_mouse_obj.data.iteritems():
    #     this_dat = col.loc[col >= min_cpm] + 1
    #     this_cpm = this_dat.divide(this_dat.sum()) * 1e6
    #     this_log_cpm = np.log2(this_cpm)
    #     this_ecdf_fun = basic.ecdf_func(this_log_cpm)
    #
    #     this_log_cpm_ecdf = this_ecdf_fun(x_cdf)
    #
    #     ens_ecdf = this_ecdf_fun(this_log_cpm[e])
    #
    #     ax.plot(x_cdf, this_log_cpm_ecdf, label=c)
    #     # optional - plot a single gene on each line
    #     ax.plot(this_log_cpm[e], ens_ecdf, 'ko')
    #
    # ax.legend(loc='lower right')
    # ax.set_xlabel("log2(CPM)")
    # ax.set_ylabel("Empirical CDF")

    # 3: Add OPCs and NSC (SS2) in
    this_data.update(dict(ss2_counts.iteritems()))
    ax = log_cpm_ecdf_plot(
        this_data,
        log_cpm_lookup_values=x_cdf,
        label_dict=first_label,
        style_dict=plot_style
    )

    ## TODO: update from here using transformations

    """

    # So how does EgdeR cope with this? Can the TMM normalisation collapse these lines?
    ss2_dat = ss2_obj.data.copy()
    human_data = pd.concat(
        (our_patient_obj.data, ss2_dat), axis=1
    )

    differential_expression.robjects.packages.importr('edgeR')
    rdata = differential_expression.pandas2ri.py2ri(human_data)

    y = differential_expression.r("DGEList")(rdata)

    yn = differential_expression.r("calcNormFactors")(y)
    cpm_2 = differential_expression.pandas2ri.ri2py_dataframe(differential_expression.r('cpm')(yn))
    cpm_2.index = human_data.index
    cpm_2.columns = human_data.columns

    # FWIW, edgeR agrees with my own calculation to within 1e-10 numerical error:
    # cpm_1 = differential_expression.pandas2ri.ri2py_dataframe(differential_expression.r('cpm')(y))
    # cpm_1.index = human_data.index
    # cpm_1.columns = human_data.columns
    # cpm_mine = human_data.divide(human_data.sum(), axis=1) * 1e6
    # ((cpm_mine - cpm_1).abs() > 1e-10).sum().sum() == 0

    normed_data = dict(cpm_2.iteritems())

    opc_lbl = False
    nsc_lbl = False
    for i, k in enumerate(ss2_obj.meta.index):
        k2 = "%s_SS2" % k
        if 'NSC' in k:
            plot_style[k2] = {'c': 'g', 'alpha': 0.6}
            if not nsc_lbl:
                first_label[k2] = 'SmartSeq2 NSC'
                nsc_lbl = True
            else:
                first_label[k2] = None
        if 'OPC' in k:
            plot_style[k2] = {'c': 'k', 'alpha': 0.6}
            if not opc_lbl:
                first_label[k2] = 'SmartSeq2 OPC'
                opc_lbl = True
            else:
                first_label[k2] = None

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for c, col in normed_data.iteritems():
        this_dat = col.loc[col >= min_cpm] + 1
        this_cpm = this_dat.divide(this_dat.sum()) * 1e6
        this_log_cpm = np.log2(this_cpm)
        this_ecdf_fun = basic.ecdf_func(this_log_cpm)
        this_log_cpm_ecdf = this_ecdf_fun(x_cdf)

        ax.plot(x_cdf, this_log_cpm_ecdf, label=first_label[c], **plot_style[c])
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")

    """