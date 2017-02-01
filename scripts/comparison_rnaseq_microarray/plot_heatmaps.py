from scripts.comparison_rnaseq_microarray import load_illumina_data, load_rnaseq_data, load_references, comparisons, consts
from plotting import heatmap
from microarray.process import aggregate_by_probe_set
import os
import pandas as pd
import numpy as np
from matplotlib import rc, pyplot as plt
import seaborn as sns
plt.interactive(True)
sns.set_style('white')

SAVE_PLOTS = True
# SAVE_PLOTS = False

if SAVE_PLOTS:
    OUTDIR = 'temp_results.0'
    i = 1
    while os.path.exists(OUTDIR):
        OUTDIR = 'temp_results.%d' % i
        i += 1
    print "Creating temp output dir %s" % OUTDIR
    os.makedirs(OUTDIR)

NORTHCOTT_C_D = consts.NORTHCOTT_GENES[2:]

# standardised vmin
ZMAX = 5.

# addition for logging
eps = 1e-12
# aggregation method for microarray
AGGR_METHOD = 'median'

# create some standard heatmap kwargs for reuse
heatmap_kwargs = {
    'vmin': -ZMAX,
    'vmax': ZMAX,
    'cbar': True,
    'orientation': 'vertical',
}

# load RNA-Seq healthy and MB
he_ct, he_tpm, he_meta = load_references.load_cerebellum_rnaseq_reference_data()
mb_tpm = load_rnaseq_data.load_rnaseq_cufflinks_gene_count_data(unit='tpm')

# change one gene name for ease of comparison
new_ix = np.array(he_tpm.index)
new_ix[new_ix == 'EYS'] = 'EGFL11'
he_tpm.index = new_ix
he_ct.index = new_ix

new_ix = np.array(mb_tpm.index)
new_ix[new_ix == 'EYS'] = 'EGFL11'
mb_tpm.index = new_ix

# log
mb_tpm_log = np.log2(mb_tpm + eps)
he_tpm_log = np.log2(he_tpm + eps)

# concatenate
all_tpm = pd.concat((mb_tpm, he_tpm), axis=1, join='inner')
all_tpm_log = np.log2(all_tpm + eps)

# standardize using mean from pool of ALL data
# all_tpm_n = all_tpm.subtract(all_tpm.mean(axis=1), axis=0).divide(all_tpm.std(axis=1), axis=0)
# all_tpm_nlog = all_tpm_log.subtract(all_tpm_log.mean(axis=1), axis=0).divide(all_tpm_log.std(axis=1), axis=0)

# standardize using mean from pool of HEALTHY data
m = he_tpm.mean(axis=1)
s = he_tpm.std(axis=1)
all_tpm_n = all_tpm.subtract(m, axis=0).divide(s, axis=0)

mlog = he_tpm_log.mean(axis=1)
slog = he_tpm_log.std(axis=1)
all_tpm_nlog = all_tpm_log.subtract(mlog, axis=0).divide(slog, axis=0)

# all_northcott = []
# [all_northcott.extend(v) for _, v in consts.NORTHCOTT_GENES]

# g = sns.clustermap(all_tpm_nlog.loc[all_northcott, :].dropna(), vmin=-3, vmax=3, row_cluster=False, col_cluster=True)
# g.cax.set_visible(False)
# plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
# plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

if SAVE_PLOTS:
    # Plot: RNA-Seq scramble vs Allen HBA, nanostring
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        consts.NANOSTRING_GENES,
        all_tpm_nlog,
        gs_kwargs={'bottom': 0.16},
        fig_kwargs={'figsize': [5.5, 8.5]},
        **heatmap_kwargs
    )
    fig.savefig(os.path.join(OUTDIR, "rnaseq_all-ahba_nanostring.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_all-ahba_nanostring.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_all-ahba_nanostring.pdf"))

if SAVE_PLOTS:
    # Plot: RNA-Seq scramble vs Allen HBA, full Northcott
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        consts.NORTHCOTT_GENES,
        all_tpm_nlog,
        fig_kwargs={'figsize': [5.5, 12]},
        heatmap_kwargs={'square': False},
        gs_kwargs={'left': 0.25},
        **heatmap_kwargs
    )
    # reduce y label font size
    for ax in axs:
        plt.setp(ax.yaxis.get_ticklabels(), fontsize=8.5)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_all-ahba_ncott.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_all-ahba_ncott.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_all-ahba_ncott.pdf"))

if SAVE_PLOTS:
    # Plot: RNA-Seq scramble vs Allen HBA, Northcott C and D only
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        NORTHCOTT_C_D,
        all_tpm_nlog,
        fig_kwargs={'figsize': [5.5, 8.5]},
        heatmap_kwargs={'square': False},
        gs_kwargs={'left': 0.25, 'bottom': 0.16},
        **heatmap_kwargs
    )
    fig.savefig(os.path.join(OUTDIR, "rnaseq_all-ahba_ncottcd.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_all-ahba_ncottcd.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_all-ahba_ncottcd.pdf"))


# repeat with the 2 scramble samples only
scr_tpm_nlog = all_tpm_nlog.loc[:, [
    u'Scramble.1', u'Scramble.2',
    u'9861_111', u'9861_112', u'9861_113', u'9861_114',
    u'10021_116', u'10021_117', u'10021_118', u'10021_119', u'10021_120'
]]

if SAVE_PLOTS:
    # Plot: RNA-Seq scramble vs Allen HBA, nanostring
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        consts.NANOSTRING_GENES,
        scr_tpm_nlog,
        gs_kwargs={'bottom': 0.16},
        fig_kwargs={'figsize': [5.5, 8.5]},
        **heatmap_kwargs
    )
    fig.savefig(os.path.join(OUTDIR, "rnaseq_scr-ahba_nanostring.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_scr-ahba_nanostring.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_scr-ahba_nanostring.pdf"))

if SAVE_PLOTS:
    # Plot: RNA-Seq scramble vs Allen HBA, full Northcott
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        consts.NORTHCOTT_GENES,
        scr_tpm_nlog,
        fig_kwargs={'figsize': [5.5, 12]},
        heatmap_kwargs={'square': False},
        gs_kwargs={'left': 0.25},
        **heatmap_kwargs
    )
    # reduce y label font size
    for ax in axs:
        plt.setp(ax.yaxis.get_ticklabels(), fontsize=8.5)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_scr-ahba_ncott.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_scr-ahba_ncott.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_scr-ahba_ncott.pdf"))

if SAVE_PLOTS:
    # Plot: RNA-Seq scramble vs Allen HBA, Northcott C and D only
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        NORTHCOTT_C_D,
        scr_tpm_nlog,
        fig_kwargs={'figsize': [5.5, 8.5]},
        heatmap_kwargs={'square': False},
        gs_kwargs={'left': 0.25, 'bottom': 0.16},
        **heatmap_kwargs
    )
    fig.savefig(os.path.join(OUTDIR, "rnaseq_scr-ahba_ncottcd.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_scr-ahba_ncottcd.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "rnaseq_scr-ahba_ncottcd.pdf"))


"""
There are currently 2 related but slightly different ways to load the Zhao data.
Here we can switch between them.
The results are qualitatively similar, but there are a few changes.
Main differences.
Method 1:
- Aggr raw before taking logs.
- Mean of repeats before taking logs.
- Add a constant to each sample's raw values to bring minimum up to 0. Log2 (x + eps).

Method 2:
- Replace any value <1 with 1. Log2(x).
- Aggregate and take mean over repeats AFTER taking logs.
"""

# METHOD 1: the old loader
LOAD_METHOD = 3

if LOAD_METHOD == 1:
    marray_data, pvals = load_illumina_data.load_normed_microarray_data(pval=None, return_pvals=True)
    # add constant to each array to force non-negative values
    marray_data = marray_data.subtract(marray_data.min(axis=0))

    probe_set = load_illumina_data.load_illumina_array_library()
    marray_ann = load_illumina_data.add_gene_symbol_column(marray_data, probe_set)
    marray_all = aggregate_by_probe_set(marray_ann, method=AGGR_METHOD)

    # take mean over repeats
    for sn in load_illumina_data.SAMPLE_NAMES:
        marray_all.loc[:, sn] = marray_all.loc[:, [sn, sn + '-R']].mean(axis=1)

    # log2
    marray_all_log = np.log2(marray_all + eps)

elif LOAD_METHOD == 2:

    # METHOD 2: the new loader
    from load_data import microarray_data
    marray_all_log, marray_meta = microarray_data.load_annotated_gse28192(aggr_field='SYMBOL', aggr_method=AGGR_METHOD)

    # take mean over repeats
    for sn in load_illumina_data.SAMPLE_NAMES:
        marray_all_log.loc[:, sn] = marray_all_log.loc[:, [sn, sn + '-R']].mean(axis=1)

elif LOAD_METHOD == 3:
    from load_data import microarray_data
    from microarray import process
    marray_raw, marray_meta = microarray_data.load_annotated_gse28192(log2=False)
    marray_all_log = process.variance_stabilizing_transform(marray_raw)
    # NB: this isn't really logged data, but it's playing the same role
    marray_all_log = microarray_data.annotate_and_aggregate_gse28192(marray_all_log, aggr_field='SYMBOL', aggr_method=AGGR_METHOD)



# extract healthy data for standardisation
HEALTHY_SAMPLES = dict(consts.SAMPLE_GROUPS_ZHAO)['Healthy cerebellum']
# marray_he = marray_all.loc[:, HEALTHY_SAMPLES]
marray_he_log = marray_all_log.loc[:, HEALTHY_SAMPLES]


# 1299 + HEALTHY

MB_SAMPLES = [
    'Pt1299',
    'ICb1299-I',
    'ICb1299-III',
    'ICb1299-IV',
]
# order here affects the order of the columns
# healthy, MB:
# keep_cols = list(HEALTHY_SAMPLES) + MB_SAMPLES
# MB, healthy:
keep_cols = MB_SAMPLES + list(HEALTHY_SAMPLES)

marr_log = marray_all_log.loc[:, keep_cols].dropna(axis=0, how='all')

# standardise using pool of HEALTHY data
marr_nlog = marr_log.subtract(marray_he_log.mean(axis=1), axis=0).divide(marray_he_log.std(axis=1), axis=0)

if SAVE_PLOTS:
    # Plot: marr MB vs healthy, nanostring
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        consts.NANOSTRING_GENES,
        marr_nlog,
        fig_kwargs={'figsize': [5, 8.5]},
        **heatmap_kwargs
    )
    fig.savefig(os.path.join(OUTDIR, "marr_1299-cer_nanostring.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_1299-cer_nanostring.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_1299-cer_nanostring.pdf"))

if SAVE_PLOTS:
    # Plot: marr MB vs healthy, Northcott full
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
            consts.NORTHCOTT_GENES,
            marr_nlog,
            fig_kwargs={'figsize': [5, 11]},
            heatmap_kwargs={'square': False},
            gs_kwargs={'left': 0.25},
            **heatmap_kwargs
        )
    # reduce y label font size
    for ax in axs:
        plt.setp(ax.yaxis.get_ticklabels(), fontsize=8.5)
    fig.savefig(os.path.join(OUTDIR, "marr_1299-cer_ncott.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_1299-cer_ncott.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_1299-cer_ncott.pdf"))

if SAVE_PLOTS:
    # Plot: marr MB vs healthy, Northcott C and D only
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        NORTHCOTT_C_D,
        marr_nlog,
        fig_kwargs={'figsize': [5, 8.5]},
        heatmap_kwargs={'square': False},
        gs_kwargs={'left': 0.25},
        **heatmap_kwargs
    )
    fig.savefig(os.path.join(OUTDIR, "marr_1299-cer_ncottcd.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_1299-cer_ncottcd.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_1299-cer_ncottcd.pdf"))


# various + HEALTHY

MB_SAMPLES = [
    'Pt1078',
    'ICb1078-I',
    'ICb1078-III',
    'ICb1078-V',
    'Pt1299',
    'ICb1299-I',
    'ICb1299-III',
    'ICb1299-IV',
    'Pt1487',
    'ICb1487-I',
    'ICb1487-III',
    'Pt1595',
    'ICb1595-I',
    'ICb1595-III',
    'Pt1338',
    'ICb1338-I',
    'ICb1338-III',
]
# order here affects the order of the columns
# healthy, MB:
# keep_cols = list(HEALTHY_SAMPLES) + MB_SAMPLES
# MB, healthy:
keep_cols = MB_SAMPLES + list(HEALTHY_SAMPLES)

marr_log = marray_all_log.loc[:, keep_cols].dropna(axis=0, how='all')

# standardise using pool of HEALTHY data
marr_nlog = marr_log.subtract(marray_he_log.mean(axis=1), axis=0).divide(marray_he_log.std(axis=1), axis=0)

if SAVE_PLOTS:
    # Plot: marr MB vs healthy, nanostring
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        consts.NANOSTRING_GENES,
        marr_nlog,
        fig_kwargs={'figsize': [5, 8.5]},
        heatmap_kwargs = {'square': False},
        **heatmap_kwargs
    )
    #  add dividing lines
    xbreaks = [4, 8, 11, 14, 17]
    for ax in axs:
        for t in xbreaks:
            ax.axvline(t, color='w', linewidth=2.5)
            ax.axvline(t, color='0.4', linewidth=1.)
    fig.savefig(os.path.join(OUTDIR, "marr_mbgrp2-cer_nanostring.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_mbgrp2-cer_nanostring.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_mbgrp2-cer_nanostring.pdf"))

if SAVE_PLOTS:
    # Plot: marr MB vs healthy, Northcott full
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
            consts.NORTHCOTT_GENES,
            marr_nlog,
            fig_kwargs={'figsize': [5, 12]},
            heatmap_kwargs={'square': False},
            gs_kwargs={'left': 0.25},
            **heatmap_kwargs
        )
    # reduce y label font size and add dividing lines
    xbreaks = [4, 8, 11, 14, 17]
    for ax in axs:
        plt.setp(ax.yaxis.get_ticklabels(), fontsize=8.5)
        for t in xbreaks:
            ax.axvline(t, color='w', linewidth=2.5)
            ax.axvline(t, color='0.4', linewidth=1.)
    fig.savefig(os.path.join(OUTDIR, "marr_mbgrp2-cer_ncott.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_mbgrp2-cer_ncott.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_mbgrp2-cer_ncott.pdf"))



# various + HEALTHY

MB_SAMPLES = [
    # 'Pt984',
    'ICb984-I',
    'ICb984-III',
    'ICb984-V',
    'Pt1338',
    'ICb1338-I',
    'ICb1338-III',
    'Pt1595',
    'ICb1595-I',
    'ICb1595-III',
    'Pt1299',
    'ICb1299-I',
    'ICb1299-III',
    'ICb1299-IV',
    'Pt1487',
    'ICb1487-I',
    'ICb1487-III',
]
# order here affects the order of the columns
# healthy, MB:
# keep_cols = list(HEALTHY_SAMPLES) + MB_SAMPLES
# MB, healthy:
keep_cols = MB_SAMPLES + list(HEALTHY_SAMPLES)

marr_log = marray_all_log.loc[:, keep_cols].dropna(axis=0, how='all')

# standardise using pool of HEALTHY data
marr_nlog = marr_log.subtract(marray_he_log.mean(axis=1), axis=0).divide(marray_he_log.std(axis=1), axis=0)


if SAVE_PLOTS:
    # Plot: marr MB vs healthy, nanostring
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
        consts.NANOSTRING_GENES,
        marr_nlog,
        fig_kwargs={'figsize': [5, 8.5]},
        heatmap_kwargs = {'square': False},
        **heatmap_kwargs
    )
    #  add dividing lines
    xbreaks = [3, 6, 9, 13, 16]
    for ax in axs:
        for t in xbreaks:
            ax.axvline(t, color='w', linewidth=2.5)
            ax.axvline(t, color='0.4', linewidth=1.)
    fig.savefig(os.path.join(OUTDIR, "marr_mb_various-cer_nanostring.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_mb_various-cer_nanostring.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_mb_various-cer_nanostring.pdf"))

if SAVE_PLOTS:
    # Plot: marr MB vs healthy, Northcott full
    fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
            consts.NORTHCOTT_GENES,
            marr_nlog,
            fig_kwargs={'figsize': [5, 12]},
            heatmap_kwargs={'square': False},
            gs_kwargs={'left': 0.25},
            **heatmap_kwargs
        )
    # reduce y label font size and add dividing lines
    xbreaks = [3, 6, 9, 13, 16]
    for ax in axs:
        plt.setp(ax.yaxis.get_ticklabels(), fontsize=8.5)
        for t in xbreaks:
            ax.axvline(t, color='w', linewidth=2.5)
            ax.axvline(t, color='0.4', linewidth=1.)
    fig.savefig(os.path.join(OUTDIR, "marr_mb_various-cer_ncott.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_mb_various-cer_ncott.tiff"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "marr_mb_various-cer_ncott.pdf"))
