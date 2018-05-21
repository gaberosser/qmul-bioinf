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
    ss2_nsc_counts = ss2_obj.data.loc[:, ss2_obj.data.columns.str.contains('NSC')]
    ss2_opc_counts = ss2_obj.data.loc[:, ss2_obj.data.columns.str.contains('OPC')]

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
        (k, {'c': 'r', 'alpha': 0.5}) for k in our_patient_obj.meta.index
    ])
    plot_style.update(
        dict([
            (k, {'c': 'b', 'alpha': 0.6, 'zorder': 99}) for k in our_mouse_obj.meta.index
        ])
    )
    plot_style.update(
        dict([
            (k, {'c': 'y', 'alpha': 0.7}) for k in our_ffpe_obj.meta.index
        ])
    )
    for i, k in enumerate(ss2_obj.meta.index):
        if 'NSC' in k:
            plot_style[k] = {'c': 'g', 'alpha': 0.6, 'zorder': 100}
        if 'OPC' in k:
            plot_style[k] = {'c': 'k', 'alpha': 0.6, 'zorder': 101}

    x_cdf = np.linspace(-6, 12.5, 500)

    # ECDFs
    # 1a) FFPE, poly(A) human, poly(A) mouse
    this_data = dict(our_patient_obj.data.iteritems())
    this_data.update(our_ffpe_obj.data.iteritems())
    this_data.update(our_mouse_obj.data.iteritems())

    ax = log_cpm_ecdf_plot(
        this_data,
        label_dict=first_label,
        # log_cpm_lookup_values=x_cdf,
        style_dict=plot_style,
        min_cpm=min_cpm
    )
    ax.figure.savefig(os.path.join(outdir, "ecdf_polya_ffpe.png"), dpi=200)

    # 1b) again with TMM norming
    polya_human_cpm_n = transformations.edger_tmm_normalisation_cpm(our_patient_obj.data)
    polya_mouse_cpm_n = transformations.edger_tmm_normalisation_cpm(our_mouse_obj.data)
    ffpe_cpm_n = transformations.edger_tmm_normalisation_cpm(our_ffpe_obj.data)
    ss2_cpm_n = transformations.edger_tmm_normalisation_cpm(ss2_obj.data)
    ss2_nsc_cpm_n = ss2_cpm_n.loc[:, ss2_cpm_n.columns.str.contains('NSC')]
    ss2_opc_cpm_n = ss2_cpm_n.loc[:, ss2_cpm_n.columns.str.contains('OPC')]

    this_data = dict(polya_human_cpm_n.iteritems())
    this_data.update(polya_mouse_cpm_n.iteritems())
    this_data.update(ffpe_cpm_n.iteritems())

    ax = log_cpm_ecdf_plot(
        this_data,
        units='cpm',
        label_dict=first_label,
        style_dict=plot_style,
        min_cpm=min_cpm
    )
    ax.figure.savefig(os.path.join(outdir, "ecdf_polya_ffpe_tmm.png"), dpi=200)

    # 2a) same, just for mouse samples

    # set the colour cycle to distinguish the samples more easily
    color.cycle_cmap(our_mouse_obj.meta.shape[0], cmap='jet')

    # set the order because it makes it easier to distinguish samples
    this_data = collections.OrderedDict([
        (c, polya_mouse_counts[c]) for c in our_mouse_obj.meta.index
    ])

    ax = log_cpm_ecdf_plot(
        this_data,
        label_dict=dict([(k, k) for k in this_data]),
        min_cpm=min_cpm,
    )

    # overlay human samples in grey
    this_data = dict(our_patient_obj.data.iteritems())
    this_plot_style = dict([
        (k, {'c': 'gray', 'alpha': 0.4, 'zorder': 1}) for k in our_patient_obj.meta.index
    ])
    ax = log_cpm_ecdf_plot(
        this_data,
        label_dict=first_label,
        style_dict=this_plot_style,
        min_cpm=min_cpm,
        ax=ax
    )
    ax.figure.savefig(os.path.join(outdir, "ecdf_mouse_samples.png"), dpi=200)

    # 2b) Again with TMM norming
    this_data = collections.OrderedDict([
        (c, polya_mouse_cpm_n[c]) for c in our_mouse_obj.meta.index
    ])

    ax = log_cpm_ecdf_plot(
        this_data,
        label_dict=dict([(k, k) for k in this_data]),
        min_cpm=min_cpm,
        units='cpm'
    )

    # overlay human samples in grey
    this_data = dict(polya_human_cpm_n.iteritems())
    this_plot_style = dict([
        (k, {'c': 'gray', 'alpha': 0.4, 'zorder': 1}) for k in our_patient_obj.meta.index
    ])
    ax = log_cpm_ecdf_plot(
        this_data,
        units='cpm',
        label_dict=first_label,
        style_dict=this_plot_style,
        min_cpm=min_cpm,
        ax=ax
    )
    ax.figure.savefig(os.path.join(outdir, "ecdf_mouse_samples_tmm.png"), dpi=200)

    # 3a) NSC (SS2) and NSC poly(A)
    this_data = dict(polya_human_counts.iteritems())
    this_data.update(dict(ss2_nsc_counts.iteritems()))
    ax = log_cpm_ecdf_plot(
        this_data,
        label_dict=first_label,
        style_dict=plot_style,
        min_cpm=min_cpm
    )
    ax.figure.savefig(os.path.join(outdir, "ecdf_nsc_two_preps.png"), dpi=200)

    # 3b) with TMM norming
    this_data = dict(polya_human_cpm_n.iteritems())
    this_data.update(dict(ss2_nsc_cpm_n.iteritems()))
    ax = log_cpm_ecdf_plot(
        this_data,
        units='cpm',
        label_dict=first_label,
        style_dict=plot_style,
        min_cpm=min_cpm
    )
    ax.figure.savefig(os.path.join(outdir, "ecdf_nsc_two_preps_tmm.png"), dpi=200)

    # 4a) NSC (SS2) and NSC poly(A) and OPC (SS2)
    this_data = dict(polya_human_counts.iteritems())
    this_data.update(dict(ss2_nsc_counts.iteritems()))
    this_data.update(dict(ss2_opc_counts.iteritems()))
    ax = log_cpm_ecdf_plot(
        this_data,
        label_dict=first_label,
        style_dict=plot_style,
        min_cpm=min_cpm
    )
    ax.figure.savefig(os.path.join(outdir, "ecdf_nsc_opc_two_preps.png"), dpi=200)

    # 4b) with TMM norming
    this_data = dict(polya_human_cpm_n.iteritems())
    this_data.update(dict(ss2_nsc_cpm_n.iteritems()))
    this_data.update(dict(ss2_opc_cpm_n.iteritems()))
    ax = log_cpm_ecdf_plot(
        this_data,
        units='cpm',
        label_dict=first_label,
        style_dict=plot_style,
        min_cpm=min_cpm
    )
    ax.figure.savefig(os.path.join(outdir, "ecdf_nsc_opc_two_preps_tmm.png"), dpi=200)