from scripts.comparison_rnaseq_microarray import load_illumina_data, load_rnaseq_data, load_references, comparisons
import references
import pandas as pd
import numpy as np
from scipy import stats
import collections
from matplotlib import rc, pyplot as plt, gridspec as gridspec
import seaborn as sns
plt.interactive(True)
sns.set_style('white')
# rc('text', usetex=True)


ALL_NORTHCOTT = [
    'WIF1',
    'TNC',
    'GAD1',
    'DKK2',
    'EMX2',
    'ADAMTSL1',
    'NKD1',
    'PDE11A',
    'EPHA7',
    'RUNX2',
    'C20orf103',
    'ABHD12B',
    'PAX3',
    'LGR5',
    'EPHA3',
    'PGM5',
    'TPH1',
    'TNFRSF19',
    'TMEM51',
    'C9orf94',
    'GABRE',
    'TNFRSF11B',
    'FREM2',
    'LEF1',
    'GABRG3',
    'PDLIM3',
    'EYA1',
    'HHIP',
    'ATOH1',
    'SFRP1',
    'PPP2R2C',
    'SOX2',
    'DMRTA1',
    'NDP',
    'GNG3',
    'NDST3',
    'GRIA4',
    'SATB2',
    'CXCR4',
    'CYYR1',
    'SCG5',
    'ALDH1A3',
    'KIF26A',
    'C6orf117',
    'BOC',
    'PRLR',
    'C4orf18',
    'CA14',
    'POU3F2',
    'SEMA6A',
    'IMPG2',
    'GABRA5',
    'EGFL11',
    'NRL',
    'MAB21L2',
    'TTR',
    'NPR3',
    'TBR1',
    'FSTL5',
    'TMEM16B',
    'C5orf23',
    'GNB3',
    'DMD',
    'PCDH21',
    'USH2A',
    'RCVRN',
    'PDE6H',
    'RASGRF2',
    'FOXG1B',
    'SAMD3',
    'TSHZ3',
    'MAK',
    'PPP2R2B',
    'RD3',
    'FAM19A4',
    'KCNA1',
    'EOMES',
    'KHDRBS2',
    'RBM24',
    'UNC5D',
    'OAS1',
    'GRM8',
    'CDH18',
    'LOC138046',
    'SNCAIP',
    'MPP3',
    'CHL1',
    'PEX5L',
    'LMX1A',
    'GPR12',
    'FBXL21',
    'SH3GL3',
    'NID2',
    'LINGO2',
    'PTHLH',
    'CA4',
    'PRL',
    'KCNIP4',
    'NEUROD2',
    'ST18',
    'OTX2',  # not Northcott, requested by SB
]

REF_GROUPS = (
    ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
    ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
    ('Group C', ('IMPG2', 'GABRA5', 'EGFL11', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
    ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1', 'OTX2')),  # OTX2 added by SB
)

MB_GROUPS = (
    ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
    ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
    ('Group C', ('IMPG2', 'GABRA5', 'EYS', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
    ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1', 'OTX2')),  # OTX2 added by SB
)

SAMPLE_GROUPS = (
    ('WNT', ('Pt1140', 'ICb1140-II', 'ICb1140-III', 'Pt1192', 'ICb1192-I', 'ICb1192-III', 'ICb1192-V')),
    ('SSH', ('Pt1338', 'ICb1338-I', 'ICb1338-III', 'ICb984-I', 'ICb984-III', 'ICb984-V')),
    ('Group C', (
        'ICb1197-I',
        'ICb1197-III',
        'Pt1494',
        'ICb1494-I',
        'ICb1494-III',
        'ICb1494-V',
        'Pt1572',
        'ICb1572-I',
        'ICb1572-III',
        'ICb1572-V',
        'Pt1595',
        'ICb1595-I',
        'ICb1595-III',
    )),
    ('Group D', (
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
    )),
)

from plotting import bar

plt = bar.plt

# addition for logging
eps = 1e-12

# standardised vmin

VMIN_N = -2.

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

# also try log version
mb_tpm_log = np.log10(mb_tpm + eps)
he_tpm_log = np.log10(he_tpm + eps)

# standardize using mean from pool of ALL data
all_tpm = mb_tpm.copy()
all_tpm = pd.concat((all_tpm, he_tpm), axis=1, join='inner')
all_tpm_log = np.log10(all_tpm + eps)

all_tpm_n = all_tpm.subtract(all_tpm.mean(axis=1), axis=0).divide(all_tpm.std(axis=1), axis=0)
mb_tpm_n = mb_tpm.subtract(all_tpm.mean(axis=1), axis=0).divide(all_tpm.std(axis=1), axis=0)

all_tpm_nlog = all_tpm_log.subtract(all_tpm_log.mean(axis=1), axis=0).divide(all_tpm_log.std(axis=1), axis=0)
mb_tpm_nlog = mb_tpm_log.subtract(all_tpm_log.mean(axis=1), axis=0).divide(all_tpm_log.std(axis=1), axis=0)

# plot standardized gene scores for Northcott groups

fig = plt.figure(figsize=[5, 4])
gs = gridspec.GridSpec(2, len(REF_GROUPS),
                       height_ratios=[1, 16],
                       width_ratios=[len(arr) for _, arr in REF_GROUPS])
gs.update(
    left=0.2,
    right=0.95,
    top=0.9,
    bottom=0.1,
    wspace=0.1,
    hspace=0.)
cbar_kws = {"orientation": "horizontal"}


for i, (grp, arr) in enumerate(REF_GROUPS):
    ax = fig.add_subplot(gs[1:, i])
    if i == (len(REF_GROUPS) - 1):
        cbar = True
        cbar_ax = fig.add_subplot(gs[0, :])
    else:
        cbar = False
        cbar_ax = None
    sns.heatmap(
        all_tpm_n.loc[arr, :].transpose(),
        ax=ax,
        vmin=VMIN_N,
        vmax=-VMIN_N,
        square=True,
        cmap='RdBu_r',
        cbar=cbar,
        cbar_ax=cbar_ax,
        cbar_kws=cbar_kws
    )
    ax.set_xticklabels(arr, rotation=90)
    if i == 0:
        plt.yticks(rotation=0)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel(grp)
    ax.xaxis.set_label_coords(.5, -.48)
cbar_ax.set_title('Standardised score by gene')

fig.savefig("rnaseq_mb_standardised_by_gene_activity_heatmap.png", dpi=200)
fig.savefig("rnaseq_mb_standardised_by_gene_activity_heatmap.pdf", dpi=200)

# plot standardized LOG gene scores for Northcott groups

fig = plt.figure(figsize=[5, 4])
gs = gridspec.GridSpec(2, len(REF_GROUPS),
                       height_ratios=[1, 16],
                       width_ratios=[len(arr) for _, arr in REF_GROUPS])
gs.update(
    left=0.2,
    right=0.95,
    top=0.9,
    bottom=0.1,
    wspace=0.1,
    hspace=0.)
cbar_kws = {"orientation": "horizontal"}

for i, (grp, arr) in enumerate(REF_GROUPS):
    ax = fig.add_subplot(gs[1:, i])
    if i == (len(REF_GROUPS) - 1):
        cbar = True
        cbar_ax = fig.add_subplot(gs[0, :])
    else:
        cbar = False
        cbar_ax = None
    sns.heatmap(
        all_tpm_nlog.loc[arr, :].transpose(),
        ax=ax,
        vmin=VMIN_N,
        vmax=-VMIN_N,
        square=True,
        cmap='RdBu_r',
        cbar=cbar,
        cbar_ax=cbar_ax,
        cbar_kws=cbar_kws
    )
    ax.set_xticklabels(arr, rotation=90)
    if i == 0:
        plt.yticks(rotation=0)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel(grp)
    ax.xaxis.set_label_coords(.5, -.48)
cbar_ax.set_title('Log10 standardised score by gene')

fig.savefig("rnaseq_mb_log_standardised_by_gene_activity_heatmap.png", dpi=200)
fig.savefig("rnaseq_mb_log_standardised_by_gene_activity_heatmap.pdf", dpi=200)

METHOD = 'median'
HEALTHY_SAMPLE_NAMES = [
    'NT1197',
    'NCb1',
    'NCb2',
    'A911105',
    'A508112',
    'A508285',
]

# load full microarray data
# marray_data = load_illumina_data.load_normed_microarray_data(pval=0.05)
marray_data, pvals = load_illumina_data.load_normed_microarray_data(pval=None, return_pvals=True)

# fill NaN with zeros
# marray_data.fillna(value=0., inplace=True)

# OR load raw data
# marray_data, meta = load_illumina_data.load_raw_microarray_data()


probe_set = load_illumina_data.load_illumina_array_library()
marray_ann = load_illumina_data.add_gene_symbol_column(marray_data, probe_set)
marray_by_gene = load_illumina_data.aggregate_by_probe_set(marray_ann, method=METHOD)

# take mean over repeats
for sn in load_illumina_data.SAMPLE_NAMES:
    marray_by_gene.loc[:, sn] = marray_by_gene.loc[:, [sn, sn + '-R']].mean(axis=1)
marray_by_gene = marray_by_gene.loc[:, load_illumina_data.SAMPLE_NAMES]

marray_by_gene_log = np.log10(marray_by_gene + eps)

## All Northcott genes

all_northcott_patients = []
for grp, arr in SAMPLE_GROUPS:
    all_northcott_patients.extend(arr)

if False:

    he = marray_by_gene.loc[:, HEALTHY_SAMPLE_NAMES].mean(axis=1)
    a = marray_by_gene.loc[ALL_NORTHCOTT, all_northcott_patients].divide(he.loc[ALL_NORTHCOTT], axis=0)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.heatmap(
        a - 1.,
        vmin=-1.,
        vmax=1.,
        square=True,
        cmap='RdBu_r',
        ax=ax
    )
    # xt = ax.get_xticklabels()
    # ax.set_xticklabels(xt, rotation=90)
    # yt = ax.get_yticklabels()
    # ax.set_yticklabels(yt, rotation=0)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)

## WEEK 5

all_mb_genes = []
for _, arr in REF_GROUPS:
    all_mb_genes.extend(arr)

# standardised scores by gene

# relative to ALL
# marray_by_gene_stand = (
#     marray_by_gene.subtract(marray_by_gene.mean(axis=1), axis=0)
#         .divide(marray_by_gene.std(axis=1), axis=0)
# )
# marray_by_gene_stand_log = (
#     marray_by_gene_log.subtract(marray_by_gene_log.mean(axis=1), axis=0)
#         .divide(marray_by_gene_log.std(axis=1), axis=0)
# )

# relative to HEALTHY
marray_by_gene_stand = (
    marray_by_gene.subtract(marray_by_gene.loc[:, HEALTHY_SAMPLE_NAMES].mean(axis=1), axis=0)
        .divide(marray_by_gene.loc[:, HEALTHY_SAMPLE_NAMES].std(axis=1), axis=0)
)
marray_by_gene_stand_log = (
    marray_by_gene_log.subtract(marray_by_gene_log.loc[:, HEALTHY_SAMPLE_NAMES].mean(axis=1), axis=0)
        .divide(marray_by_gene_log.loc[:, HEALTHY_SAMPLE_NAMES].std(axis=1), axis=0)
)

# v1: absolute values
if False:
    VMAX = 15000

    fig = plt.figure(figsize=[5, 8])
    gs = gridspec.GridSpec(2, len(REF_GROUPS),
                           height_ratios=[1, 12],
                           width_ratios=[len(arr) for _, arr in REF_GROUPS])
    gs.update(
        left=0.2,
        right=0.95,
        top=0.95,
        bottom=0.15,
        wspace=0.,
        hspace=0.1)
    cbar_kws = {"orientation": "horizontal"}


    for i, (grp, arr) in enumerate(REF_GROUPS):
        ax = fig.add_subplot(gs[1:, i])
        if i == (len(REF_GROUPS) - 1):
            cbar = True
            cbar_ax = fig.add_subplot(gs[0, :])
        else:
            cbar = False
            cbar_ax = None
        sns.heatmap(
            marray_by_gene.loc[arr, all_northcott_patients + HEALTHY_SAMPLE_NAMES].transpose(),
            ax=ax,
            vmin=0,
            vmax=VMAX,
            square=True,
            cmap='Reds',
            cbar=cbar,
            cbar_ax=cbar_ax,
            cbar_kws=cbar_kws
        )
        ax.set_xticklabels(arr, rotation=90)
        if i == 0:
            plt.yticks(rotation=0)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel(grp)
        ax.xaxis.set_label_coords(.5, -.15)
    cbar_ax.set_title('$\log_2$(Normalised intensity)')

    fig.savefig("marray_all_samples_mb_gene_activity_heatmap.png", dpi=200)
    fig.savefig("marray_all_samples_mb_gene_activity_heatmap.pdf", dpi=200)

# v2: log2(absolute) values

if False:
    fig = plt.figure(figsize=[5, 8])
    gs = gridspec.GridSpec(2, len(REF_GROUPS),
                           height_ratios=[1, 12],
                           width_ratios=[len(arr) for _, arr in REF_GROUPS])
    gs.update(
        left=0.2,
        right=0.95,
        top=0.95,
        bottom=0.15,
        wspace=0.,
        hspace=0.1)
    cbar_kws = {"orientation": "horizontal"}


    for i, (grp, arr) in enumerate(REF_GROUPS):
        ax = fig.add_subplot(gs[1:, i])
        if i == (len(REF_GROUPS) - 1):
            cbar = True
            cbar_ax = fig.add_subplot(gs[0, :])
        else:
            cbar = False
            cbar_ax = None
        sns.heatmap(
            np.log2(marray_by_gene.loc[arr, all_northcott_patients + HEALTHY_SAMPLE_NAMES].transpose()),
            ax=ax,
            vmin=0,
            vmax=np.ceil(np.log2(VMAX)),
            square=True,
            cmap='Reds',
            cbar=cbar,
            cbar_ax=cbar_ax,
            cbar_kws=cbar_kws
        )
        ax.set_xticklabels(arr, rotation=90)
        if i == 0:
            plt.yticks(rotation=0)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel(grp)
        ax.xaxis.set_label_coords(.5, -.15)
    cbar_ax.set_title('$\log_2$(Normalised intensity)')

    fig.savefig("marray_all_samples_mb_gene_log_activity_heatmap.png", dpi=200)
    fig.savefig("marray_all_samples_mb_gene_log_activity_heatmap.pdf", dpi=200)

# v3: standardised by score values

fig = plt.figure(figsize=[5, 8])
gs = gridspec.GridSpec(2, len(REF_GROUPS),
                       height_ratios=[1, 12],
                       width_ratios=[len(arr) for _, arr in REF_GROUPS])
gs.update(
    left=0.2,
    right=0.95,
    top=0.95,
    bottom=0.15,
    wspace=0.,
    hspace=0.1)
cbar_kws = {"orientation": "horizontal"}


for i, (grp, arr) in enumerate(REF_GROUPS):
    ax = fig.add_subplot(gs[1:, i])
    if i == (len(REF_GROUPS) - 1):
        cbar = True
        cbar_ax = fig.add_subplot(gs[0, :])
    else:
        cbar = False
        cbar_ax = None
    sns.heatmap(
        marray_by_gene_stand.loc[arr, all_northcott_patients + HEALTHY_SAMPLE_NAMES].transpose(),
        ax=ax,
        vmin=-2,
        vmax=2.,
        square=True,
        cmap='RdBu_r',
        cbar=cbar,
        cbar_ax=cbar_ax,
        cbar_kws=cbar_kws
    )
    ax.set_xticklabels(arr, rotation=90)
    if i == 0:
        plt.yticks(rotation=0)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel(grp)
    ax.xaxis.set_label_coords(.5, -.15)
cbar_ax.set_title('Standardised score by gene')

fig.savefig("marray_all_samples_mb_standardised_by_gene_activity_heatmap.png", dpi=200)
fig.savefig("marray_all_samples_mb_standardised_by_gene_activity_heatmap.pdf", dpi=200)


# v3b: as v3 but logged data

fig = plt.figure(figsize=[5, 8])
gs = gridspec.GridSpec(2, len(REF_GROUPS),
                       height_ratios=[1, 12],
                       width_ratios=[len(arr) for _, arr in REF_GROUPS])
gs.update(
    left=0.2,
    right=0.95,
    top=0.95,
    bottom=0.15,
    wspace=0.,
    hspace=0.1)
cbar_kws = {"orientation": "horizontal"}


for i, (grp, arr) in enumerate(REF_GROUPS):
    ax = fig.add_subplot(gs[1:, i])
    if i == (len(REF_GROUPS) - 1):
        cbar = True
        cbar_ax = fig.add_subplot(gs[0, :])
    else:
        cbar = False
        cbar_ax = None
    sns.heatmap(
        marray_by_gene_stand_log.loc[arr, all_northcott_patients + HEALTHY_SAMPLE_NAMES].transpose(),
        ax=ax,
        vmin=-2,
        vmax=2.,
        square=True,
        cmap='RdBu_r',
        cbar=cbar,
        cbar_ax=cbar_ax,
        cbar_kws=cbar_kws
    )
    ax.set_xticklabels(arr, rotation=90)
    if i == 0:
        plt.yticks(rotation=0)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel(grp)
    ax.xaxis.set_label_coords(.5, -.15)
cbar_ax.set_title('Log10 standardised score by gene')

fig.savefig("marray_all_samples_mb_log_standardised_by_gene_activity_heatmap.png", dpi=200)
fig.savefig("marray_all_samples_mb_log_standardised_by_gene_activity_heatmap.pdf", dpi=200)

# v4: combine with RNA-Seq

fig = plt.figure(figsize=[5, 8])
gs = gridspec.GridSpec(2, len(REF_GROUPS),
                       height_ratios=[1, 12],
                       width_ratios=[len(arr) for _, arr in REF_GROUPS])
gs.update(
    left=0.2,
    right=0.95,
    top=0.95,
    bottom=0.15,
    wspace=0.,
    hspace=0.1)
cbar_kws = {"orientation": "horizontal"}


for i, (grp, arr) in enumerate(REF_GROUPS):
    m = marray_by_gene_stand.loc[arr, all_northcott_patients + HEALTHY_SAMPLE_NAMES]
    r = mb_tpm_n.loc[arr, :]
    a = pd.concat((m, r), axis=1)

    ax = fig.add_subplot(gs[1:, i])
    if i == (len(REF_GROUPS) - 1):
        cbar = True
        cbar_ax = fig.add_subplot(gs[0, :])
    else:
        cbar = False
        cbar_ax = None
    sns.heatmap(
        a.transpose(),
        ax=ax,
        vmin=-2,
        vmax=2.,
        square=True,
        cmap='RdBu_r',
        cbar=cbar,
        cbar_ax=cbar_ax,
        cbar_kws=cbar_kws
    )
    ax.set_xticklabels(arr, rotation=90)
    if i == 0:
        plt.yticks(rotation=0)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel(grp)
    ax.xaxis.set_label_coords(.5, -.15)
cbar_ax.set_title('Standardised score by gene')

fig.savefig("marray_and_rnaseq_all_samples_mb_standardised_by_gene_activity_heatmap.png", dpi=200)
fig.savefig("marray_and_rnaseq_all_samples_mb_standardised_by_gene_activity_heatmap.pdf", dpi=200)


# v4b: as v4 but log

fig = plt.figure(figsize=[5, 8])
gs = gridspec.GridSpec(2, len(REF_GROUPS),
                       height_ratios=[1, 12],
                       width_ratios=[len(arr) for _, arr in REF_GROUPS])
gs.update(
    left=0.2,
    right=0.95,
    top=0.95,
    bottom=0.15,
    wspace=0.,
    hspace=0.1)
cbar_kws = {"orientation": "horizontal"}


for i, (grp, arr) in enumerate(REF_GROUPS):
    m = marray_by_gene_stand_log.loc[arr, all_northcott_patients + HEALTHY_SAMPLE_NAMES]
    r = mb_tpm_n.loc[arr, :]
    a = pd.concat((m, r), axis=1)

    ax = fig.add_subplot(gs[1:, i])
    if i == (len(REF_GROUPS) - 1):
        cbar = True
        cbar_ax = fig.add_subplot(gs[0, :])
    else:
        cbar = False
        cbar_ax = None
    sns.heatmap(
        a.transpose(),
        ax=ax,
        vmin=-2,
        vmax=2.,
        square=True,
        cmap='RdBu_r',
        cbar=cbar,
        cbar_ax=cbar_ax,
        cbar_kws=cbar_kws
    )
    ax.set_xticklabels(arr, rotation=90)
    if i == 0:
        plt.yticks(rotation=0)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel(grp)
    ax.xaxis.set_label_coords(.5, -.15)
cbar_ax.set_title('Log10 standardised score by gene')

fig.savefig("marray_and_rnaseq_all_samples_mb_log_standardised_by_gene_activity_heatmap.png", dpi=200)
fig.savefig("marray_and_rnaseq_all_samples_mb_log_standardised_by_gene_activity_heatmap.pdf", dpi=200)


def box_and_scatter(a, b1, b2, ylabel=None):
    fig = plt.figure(figsize=(5, 7))
    ax = fig.add_subplot(111)
    sns.boxplot(data=a, ax=ax, flierprops={'markersize': 0})
    sns.stripplot(data=b1, size=6, jitter=True, edgecolor='k', lw=1, ax=ax)
    sns.stripplot(data=b2, size=6, jitter=True, color='w', edgecolor='k', lw=1, ax=ax)
    ax.set_xlabel(grp)
    if ylabel:
        ax.set_ylabel(ylabel)
    plt.tight_layout()
    return fig, ax


def box_and_scatter_activity(grp, arr):
    a = all_tpm.loc[arr, :].transpose()
    b1 = he_tpm.loc[arr, :].transpose()
    b2 = mb_tpm.loc[arr, :].transpose()
    fig, ax = box_and_scatter(a, b1, b2, '# reads')

    fig.savefig("marray_and_rnaseq_box_and_scatter_activity_%s.png" % grp, dpi=200)
    fig.savefig("marray_and_rnaseq_box_and_scatter_activity_%s.pdf" % grp, dpi=200)


def box_and_scatter_log_activity(grp, arr):
    a = np.log10(all_tpm.loc[arr, :] + 1e-12).transpose()
    b1 = np.log10(he_tpm.loc[arr, :] + 1e-12).transpose()
    b2 = np.log10(mb_tpm.loc[arr, :] + 1e-12).transpose()
    fig, ax = box_and_scatter(a, b1, b2, '$\log_{10}($# reads)')

    fig.savefig("marray_and_rnaseq_box_and_scatter_log10_activity_%s.png" % grp, dpi=200)
    fig.savefig("marray_and_rnaseq_box_and_scatter_log10_activity_%s.pdf" % grp, dpi=200)

sns.set_style('whitegrid')
for grp, arr in REF_GROUPS:
    box_and_scatter_activity(grp, arr)
    box_and_scatter_log_activity(grp, arr)


## WEEK 4

if False:

    mb_samples = ('ICb1299-III', 'ICb1299-IV')
    mb_sample_names = list(mb_samples) + [t + '-R' for t in mb_samples]

    # pick out samples and aggregate
    mb = marray_by_gene.loc[:, mb_sample_names].mean(axis=1)
    he = marray_by_gene.loc[:, HEALTHY_SAMPLE_NAMES].mean(axis=1)

    # load reference RNA-Seq
    rna_ref, meta = load_references.load_cerebellum_rnaseq_reference_data()


    # lower threshold by percentile (not including zeros)
    # fmin_marr = 0.001

    # lower threshold by percentile
    # he_min = he[he != 0].sort_values().ix[int(he.size * fmin_marr)]

    # standardise healthy microarray results
    st_he = he[he != 0].sort_values()
    # st_he[st_he < he_min] = he_min
    st_he = np.log2(st_he)
    st_he /= st_he.max()

    # standardise

    # take mean of reference samples and sort by ascending value
    rna_ref_mean = rna_ref.mean(axis=1).sort_values()
    # lower threshold by percentile
    # eps = rna_ref_mean[rna_ref_mean != 0].ix[int(rna_ref_mean.size * fmin_rna)]

    # standardise RNA-Seq
    # either ignore all zeros, because that's what we did for microarray:
    # st_rna = rna_ref_mean.loc[rna_ref_mean != 0.].copy()
    # or keep them:
    st_rna = rna_ref_mean.copy()  # Zhao et al keep them then threshold

    # floor the values, either by setting a fixed threshold:
    # rna_floor = 5e-8
    # or by using a percentile lookup on the non-zero values:
    fmin_rna = 0.01
    rna_floor = rna_ref_mean[rna_ref_mean != 0].ix[int(rna_ref_mean.size * fmin_rna)]
    # or by using a percentile lookup on ALL values:
    # fmin_rna = 0.3  # approx value used by Zhao et al.
    # rna_floor = rna_ref_mean.ix[int(rna_ref_mean.size * fmin_rna)]


    # threshold all values below the floor
    st_rna[st_rna < rna_floor] = rna_floor

    # standardise
    rna_min = st_rna.min()

    # ensure all values > 1
    # st_rna[st_rna < rna_min] = rna_min
    st_rna = np.log2(st_rna / rna_min)
    st_rna /= st_rna.max()


    # plot 'dynamic range'
    rele_he = np.linspace(0, 1, st_he.size)
    rele_rna = np.linspace(0, 1, st_rna.size)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(rele_he, st_he, 'k', label='Microarray')
    ax.plot(rele_rna, st_rna, 'r', label='RNA-Seq')
    ax.set_xlabel('Percentile rank')
    ax.set_ylabel('Standardised log2 expression')
    ax.legend(loc='upper left')

    # plot scatter
    joint_ind = st_he.index.intersection(st_rna.index)
    obj = stats.linregress(st_he.loc[joint_ind], st_rna.loc[joint_ind])

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    plt.scatter(st_he.loc[joint_ind], st_rna.loc[joint_ind], c='b', alpha=0.4)
    plt.plot([0, 1], [obj.intercept, obj.intercept + obj.slope], 'r', label='$R = %.3f$' % obj.rvalue)
    ax.set_xlim([-1e-2, 1 + 1e-2])
    ax.set_ylim([-1e-2, 1 + 1e-2])
    ax.set_aspect('equal')
    ax.set_xlabel('Standardised microarray expression')
    ax.set_ylabel('Standardised RNA-Seq expression')
    ax.legend(fontsize=14, loc='upper left')


    ## comparing by only a subregion, selected by defining bounds on the rank percentiles
    fmin = 0.25
    fmax = 0.75

    a = rna_ref_mean.sort_values().ix[int(fmin * rna_ref_mean.size):int(fmax * rna_ref_mean.size)]
    b = he.sort_values().ix[int(fmin * he.size):int(fmax * he.size)]

    # joint index
    joint_ind = a.index.intersection(b.index)
    a = a.loc[joint_ind]
    b = b.loc[joint_ind]

    # fig = plt.figure()
    # plt.scatter(a, b)
