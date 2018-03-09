import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from scipy import stats
from statsmodels import api as sm
from sklearn.preprocessing import Imputer

from settings import GIT_LFS_DATA_DIR
from utils.output import unique_output_dir
from plotting import common

BIOMARKERS = [
    'Hb',
    'Plt',
    'Neutrophil',
    'Lymphocyte',
    # 'PT',  # PT test and INR are highly related
    'INR',
    # 'APTT',  # too many missing values
    'Urea',
    'Creatinine',
    'ALT',
    'ALP',
    'CRP'
]
BIOMARKER_PEAK_COLS = ["%s peak" % t for t in BIOMARKERS]
BIOMARKER_PEAK_AGE_COLS = ["%s peak age" % t for t in BIOMARKERS]
BIOMARKER_TROUGH_COLS = ['Plt trough']
BIOMARKER_TROUGH_AGE_COLS = ['Plt trough age']

OUTCOME_COL = 'Outcome'


def load_cleaned_data():
    # fn = os.path.join(GIT_LFS_DATA_DIR, 'divyen_shah', 'cleaned_data_jan_2018.csv')
    fn = os.path.join(GIT_LFS_DATA_DIR, 'divyen_shah', 'cleaned_data_feb_2018.csv')
    dat = pd.read_csv(fn, header=0, na_values='-', index_col=0)
    dat.loc[:, 'batch'] = [t[:2] for t in dat.index]
    return dat


def load_full_data():
    fn = os.path.join(GIT_LFS_DATA_DIR, 'divyen_shah', 'cleaned_data_full_cohort_feb_2018.csv')
    dat = pd.read_csv(fn, header=0, na_values='-', index_col=0)
    dat.loc[:, 'batch'] = [t[:2] for t in dat.index]
    return dat


def impute_missing(data, strategy='median'):
    X = data.copy()
    imp = Imputer(missing_values='NaN', strategy=strategy, axis=0)
    X = imp.fit_transform(X)
    return pd.DataFrame(X, index=data.index, columns=data.columns)


def scatter_two_vars_two_categories(dat1, dat2, cls1, cls2, ecs, cs):

    i1, lbl1 = pd.factorize(cls1)
    i2, lbl2 = pd.factorize(cls2)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for (j1, l1) in enumerate(lbl1):
        for (j2, l2) in enumerate(lbl2):
            idx = (i1 == j1) & (i2 == j2)
            ax.scatter(
                dat1.loc[idx],
                dat2.loc[idx],
                color=cs[l1],
                edgecolor=ecs[l2],
                linewidths=1.5,
                label='%s, %s' % (l1, l2)
            )

    ax.set_xlabel(dat1.name)
    ax.set_ylabel(dat2.name)
    ax.legend(loc='upper right')
    fig.tight_layout()

    return fig, ax


def get_stat_asterisks(pval):
    if pval <= 0.0001:
        return "****"
    if pval <= 0.001:
        return "***"
    if pval <= 0.01:
        return "**"
    if pval <= 0.05:
        return "*"

if __name__ == "__main__":

    outdir = unique_output_dir("hie_peak_times", reuse_empty=True)
    dat = load_cleaned_data()

    biomarkers = dat.loc[:, (
        BIOMARKER_PEAK_COLS
        + BIOMARKER_TROUGH_COLS
        + BIOMARKER_PEAK_AGE_COLS
        + BIOMARKER_TROUGH_AGE_COLS
    )]
    outcomes = dat.loc[:, OUTCOME_COL]
    peaks_dat = dat.loc[:, BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS]
    nvar = peaks_dat.shape[1]

    platelet_idx = dat.loc[:, 'Blood product'].str.contains('Platelets').fillna(False)

    lightblue = '#5aa1c4'
    lightgreen = '#64ba51'

    ncols = nvar / 2 + nvar % 2
    fig, axs = plt.subplots(nrows=2, ncols=ncols, sharex=True, figsize=(12, 7))
    axs = axs.flat

    age_pval_mwu = {}

    for i, c in enumerate(peaks_dat.columns):
        age_col = "%s age" % c
        ax = axs[i]

        this_dat = dat.loc[:, [age_col, 'Outcome']]
        age_pval_mwu[age_col] = stats.mannwhitneyu(this_dat.loc[this_dat.Outcome == 2, age_col].values,
                                             this_dat.loc[this_dat.Outcome == 1, age_col].values).pvalue

        this_dat.loc[this_dat.loc[:, 'Outcome'] == 1, 'Outcome'] = 'Unfavourable'
        this_dat.loc[this_dat.loc[:, 'Outcome'] == 2, 'Outcome'] = 'Favourable'

        sns.boxplot(data=this_dat, x='Outcome', y=age_col, ax=ax, color='w', fliersize=0)
        sns.swarmplot(data=this_dat, x='Outcome', y=age_col, ax=ax, alpha=0.8)

        ax.xaxis.label.set_visible(False)
        ax.yaxis.label.set_visible(False)

        ttl = age_col
        stars = get_stat_asterisks(age_pval_mwu[age_col])
        if stars is not None:
            ttl = "%s (%s)" % (age_col, stars)

        ax.set_title(ttl)
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
        # ax.set_xticklabels(['Fav', 'Unfav'], rotation=45)
        ylim = list(ax.get_ylim())
        ylim[0] = -5
        ax.set_ylim(ylim)
        # if ylim[0] < 0:
        #     ylim[0] = 0.
        #     ax.set_ylim(ylim)

    fig.tight_layout()

    if nvar % 2 == 1:
        axs[-1].set_visible(False)

    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter_peak_time.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter_peak_time.pdf'))

    # replot the Plt peak case with platelet idx marked
    fig, axs = plt.subplots(1, 2, sharey=True)
    c = 'Plt peak'
    age_col = "%s age" % c

    this_dat = dat.loc[:, [age_col, 'Outcome']]

    p_with = age_pval_mwu[age_col]
    p_without = stats.mannwhitneyu(this_dat.loc[(this_dat.Outcome == 2) & (~platelet_idx), age_col].values,
                                   this_dat.loc[(this_dat.Outcome == 1) & (~platelet_idx), age_col].values).pvalue

    this_dat.loc[this_dat.loc[:, 'Outcome'] == 1, 'Outcome'] = 'Unfavourable'
    this_dat.loc[this_dat.loc[:, 'Outcome'] == 2, 'Outcome'] = 'Favourable'

    sns.boxplot(data=this_dat, x='Outcome', y=age_col, ax=axs[0], color='w', fliersize=0)
    sns.swarmplot(data=this_dat.loc[~platelet_idx], x='Outcome', y=age_col, ax=axs[0], alpha=0.8)
    sns.swarmplot(data=this_dat.loc[platelet_idx], x='Outcome', y=age_col, edgecolor='k', linewidth=1.5, ax=axs[0],
                  alpha=0.8)

    sns.boxplot(data=this_dat.loc[~platelet_idx], x='Outcome', y=age_col, ax=axs[1], color='w', fliersize=0)
    sns.swarmplot(data=this_dat.loc[~platelet_idx], x='Outcome', y=age_col, ax=axs[1], alpha=0.8)

    ttl_with = "%s, all (p=%.3f)" % (age_col, p_with)
    ttl_without = "%s, not given platelets (p=%.3f)" % (age_col, p_without)

    axs[0].set_title(ttl_with)
    axs[1].set_title(ttl_without)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'plt_peak_time_with_without.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'plt_peak_time_with_without.pdf'))

    # t test and MWU test
    # look for differences in the distribution of individual variables between outcomes
    p_ttest = {}
    p_mwu_test = {}

    for col in peaks_dat.columns:
        this_dat = peaks_dat.loc[:, col].groupby(outcomes).apply(lambda x: x.dropna().values)
        args = []
        this_dropped = []
        for k, t in this_dat.iteritems():
            if len(t) == 0:
                print "One of the 2 outcome groups is missing variable %s. Skipping" % col
                continue
            else:
                args.append(t)
        this_res = stats.ttest_ind(*args)
        p_ttest[col] = this_res.pvalue
        this_res = stats.mannwhitneyu(*args)
        p_mwu_test[col] = this_res.pvalue

    # a few different representations
    o = pd.Series('Unfavourable', index=peaks_dat.index)
    o.loc[outcomes == 2] = 'Favourable'

    p = pd.Series('not given platelets', index=peaks_dat.index)
    p[platelet_idx] = 'given platelets'

    ecs = {
        'given platelets': 'k',
        'not given platelets': 'none',
    }

    cs = {
        'Favourable': lightblue,
        'Unfavourable': lightgreen
    }

    fig, ax = scatter_two_vars_two_categories(
        peaks_dat.loc[:, 'Plt peak'],
        peaks_dat.loc[:, 'CRP peak'],
        o,
        p,
        ecs=ecs,
        cs=cs
    )
    ax.legend(loc='upper right', frameon=True, facecolor='w', edgecolor='k')
    fig.savefig(os.path.join(outdir, 'plt_peak_vs_crp_peak.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'plt_peak_vs_crp_peak.pdf'))

    plt_ratio = peaks_dat.loc[:, 'Plt peak'] / peaks_dat.loc[:, 'Plt trough']
    plt_ratio.name = 'Plt peak / Plt trough'
    fig, ax = scatter_two_vars_two_categories(
        plt_ratio,
        peaks_dat.loc[:, 'CRP peak'],
        o,
        p,
        ecs=ecs,
        cs=cs
    )
    ax.legend(loc='lower right', frameon=True, facecolor='w', edgecolor='k')
    fig.savefig(os.path.join(outdir, 'plt_ratio_vs_crp_peak.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'plt_ratio_vs_crp_peak.pdf'))

    fig, ax = scatter_two_vars_two_categories(
        plt_ratio[~platelet_idx],
        peaks_dat.loc[~platelet_idx, 'CRP peak'],
        o[~platelet_idx],
        p[~platelet_idx],
        ecs=ecs,
        cs=cs
    )
    ax.legend(loc='lower right', frameon=True, facecolor='w', edgecolor='k')
    fig.savefig(os.path.join(outdir, 'plt_ratio_vs_crp_peak_no_platelets.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'plt_ratio_vs_crp_peak_no_platelets.pdf'))

    plt_diff = peaks_dat.loc[:, 'Plt peak'] - peaks_dat.loc[:, 'Plt trough']
    plt_diff.name = 'Plt peak - Plt trough'
    fig, ax = scatter_two_vars_two_categories(
        plt_diff,
        peaks_dat.loc[:, 'CRP peak'],
        o,
        p,
        ecs=ecs,
        cs=cs
    )
    ax.legend(loc='upper right', frameon=True, facecolor='w', edgecolor='k')
    fig.savefig(os.path.join(outdir, 'plt_diff_vs_crp_peak.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'plt_diff_vs_crp_peak.pdf'))

    # Peak values again

    ncols = nvar / 2 + nvar % 2
    fig, axs = plt.subplots(nrows=2, ncols=ncols, sharex=True, figsize=(12, 7))
    axs = axs.flat

    for i, c in enumerate(peaks_dat.columns):
        ax = axs[i]

        this_dat = dat.loc[:, [c, 'Outcome']]
        pval = p_mwu_test[c]

        this_dat.loc[this_dat.loc[:, 'Outcome'] == 1, 'Outcome'] = 'Unfavourable'
        this_dat.loc[this_dat.loc[:, 'Outcome'] == 2, 'Outcome'] = 'Favourable'

        sns.boxplot(data=this_dat, x='Outcome', y=c, ax=ax, color='w', fliersize=0)
        sns.swarmplot(data=this_dat, x='Outcome', y=c, ax=ax, alpha=0.8)

        ax.xaxis.label.set_visible(False)
        ax.yaxis.label.set_visible(False)

        ttl = c
        stars = get_stat_asterisks(pval)
        if stars is not None:
            ttl = "%s (%s)" % (c, stars)

        ax.set_title(ttl)
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
        ylim = list(ax.get_ylim())
        if ylim[0] <= 0:
            ylim[0] = -5.
            ax.set_ylim(ylim)

    fig.tight_layout()

    if nvar % 2 == 1:
        axs[-1].set_visible(False)

    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter_peak_values.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter_peak_values.pdf'))