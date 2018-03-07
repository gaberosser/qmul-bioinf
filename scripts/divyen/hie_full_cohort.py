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
    'INR',
    'APTT',
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
    fn = os.path.join(GIT_LFS_DATA_DIR, 'divyen_shah', 'cleaned_data_full_cohort_feb_2018.csv')
    return pd.read_csv(fn, header=0, na_values='-', index_col=0)


def standardise(data, axis=0):
    if axis == 0:
        return data.subtract(data.mean(axis=0), axis=1).divide(data.std(axis=0), axis=1)
    elif axis == 1:
        return data.subtract(data.mean(axis=1), axis=0).divide(data.std(axis=1), axis=0)
    else:
        raise AttributeError("Axis must be 0 (norm by col) or 1 (norm by row)")


def impute_missing(data, strategy='median'):
    X = data.copy()
    imp = Imputer(missing_values='NaN', strategy=strategy, axis=0)
    X = imp.fit_transform(X)
    return pd.DataFrame(X, index=data.index, columns=data.columns)


if __name__ == "__main__":
    outdir = unique_output_dir("hie_full_cohort_results", reuse_empty=True)

    dat = load_cleaned_data()
    dat.loc[:, 'batch'] = [t[:2] for t in dat.index]
    biomarkers = dat.loc[:, (
        BIOMARKER_PEAK_COLS
        + BIOMARKER_TROUGH_COLS
        + BIOMARKER_PEAK_AGE_COLS
        + BIOMARKER_TROUGH_AGE_COLS
    )]
    outcomes = dat.loc[:, OUTCOME_COL]
    peaks_dat = dat.loc[:, BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS]
    nvar = peaks_dat.shape[1]
    X = impute_missing(peaks_dat, strategy='median')

    meconium_idx = dat.loc[:, 'Meconium Aspiration'] == 'Y'
    culture_idx = dat.loc[:, 'Culture'] == 'Y'
    mec_index = meconium_idx | culture_idx

    # scatter + boxplot showing distribution of variables between mec / nonmec
    if nvar % 2 == 0:
        fig, axs = plt.subplots(nrows=2, ncols=nvar / 2, sharex=True, figsize=(12, 8))
    else:
        fig, axs = plt.subplots(nrows=2, ncols=(nvar + 1) / 2, sharex=True, figsize=(12, 8))
    axs = axs.flat

    for i, c in enumerate(peaks_dat.columns):
        ax = axs[i]

        this_dat = peaks_dat.loc[:, [c]]
        this_dat.insert(1, 'Mec', meconium_idx | culture_idx)

        # this_dat = dat.loc[:, [c, 'Outcome']]
        idx = this_dat.loc[:, 'Mec']
        this_dat.loc[idx, 'Mec'] = 'Meconium / culture'
        this_dat.loc[~idx, 'Mec'] = 'Normal'

        sns.boxplot(data=this_dat, x='Mec', y=c, ax=ax, color='w')
        sns.swarmplot(data=this_dat, x='Mec', y=c, ax=ax)

        # ax.scatter(scat[:, 0], scat[:, 1], c=scat[:, 2], cmap='RdBu')

        ax.xaxis.label.set_visible(False)
        ax.yaxis.label.set_visible(False)
        ax.set_title(c)
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
        # ax.set_xticklabels(['Fav', 'Unfav'], rotation=45)
        ylim = list(ax.get_ylim())
        if ylim[0] < 0:
            ylim[0] = 0.
            ax.set_ylim(ylim)
    # fig.subplots_adjust(left=0.025, right=0.99, wspace=0.4, bottom=0.2)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter.pdf'))

    # t test
    # look for differences in the distribution of individual variables between outcomes
    p_ttest = {}

    for col in peaks_dat.columns:
        this_dat = peaks_dat.loc[:, col].groupby(mec_index).apply(lambda x: list(x.dropna()))
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

    # to what extent do biomarkers predict meconium aspiration / infection?
    mc_outcome = (meconium_idx | culture_idx).astype(int)
    logit_model = sm.Logit(mc_outcome, sm.add_constant(X))  # outcome is mec/cult positive
    result = logit_model.fit()
    print "Predicting MA/cult. status using biomarkers:"
    print result.summary()
    ci = result.conf_int()
    ci.columns = ["2.5%", "97.5%"]
    ci.insert(0, "Odds ratio", result.params)
    ci = np.exp(ci)
    ci.insert(0, "P value", result.pvalues)

    # CRP peak is significant: plot the box and whisker
    jitter = 0.2
    lbl = 'Meconium aspiration or culture'
    crp_pk = peaks_dat.loc[:, ['CRP peak']]
    crp_pk.insert(0, lbl, (meconium_idx | culture_idx))

    crp_pk.loc[(meconium_idx | culture_idx), lbl] = 'Positive'
    crp_pk.loc[~(meconium_idx | culture_idx), lbl] = 'Negative'

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.boxplot(data=crp_pk, x=lbl, y='CRP peak', ax=ax, color='w')
    sns.swarmplot(data=crp_pk, x=lbl, y='CRP peak', ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "crp_peak_vs_meconium_or_culture.png"), dpi=200)

    # one way ANOVA
    # look for sign differences in the distribution of individual variables between batches
    p_anova = {}
    for col in peaks_dat.columns:
        this_dat = peaks_dat.loc[:, col].groupby(dat.batch).apply(lambda x: list(x.dropna()))
        args = []
        this_dropped = []
        for k, t in this_dat.iteritems():
            if len(t) > 0:
                args.append(t)
            else:
                this_dropped.append(k)
        if len(this_dropped):
            # find the missing ones
            print "Warning: variable %s is entirely missing from batches: %s" % (col, ", ".join(this_dropped))
            print "Running with remaining batches"
        this_res = stats.f_oneway(*args)
        p_anova[col] = this_res.pvalue
        if this_res.pvalue < 0.05:
            print "Variable %s has a significant one-way ANOVA result (p<%.3e)" % (col, this_res.pvalue)

    # for each significantly different biomarker, plot a histogram
    b_sign = False
    for col, p in p_anova.items():
        if p < 0.05:
            b_sign = True
            print "Variable %s has a significant one-way ANOVA result (p<%.3e)" % (col, p)
            this_dat = pd.concat((peaks_dat.loc[:, col], dat.batch), axis=1)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sns.boxplot(data=this_dat, x='batch', y=col, ax=ax, color='w')
            sns.swarmplot(data=this_dat, y=col, x='batch', ax=ax)
            fig.savefig(os.path.join(outdir, "sign_diff_by_anova_%s.png" % col), dpi=200)

    if not b_sign:
        print "No significant differences detected between batches (one-way ANOVA)"

    # t test
    # look for differences in the distribution of individual variables between outcomes
    p_ttest = {}

    for col in peaks_dat.columns:
        this_dat = peaks_dat.loc[:, col].groupby(outcomes).apply(lambda x: list(x.dropna()))
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


    # for each variable, plot a box and whisker for the two outcome classes
    jitter = 0.2
    if nvar % 2 == 0:
        fig, axs = plt.subplots(nrows=2, ncols=nvar / 2, sharex=True, figsize=(12, 8))
        axs = axs.flat
    else:
        fig, axs = plt.subplots(nrows=1, ncols=nvar, sharex=True, figsize=(15, 5))

    for i, c in enumerate(peaks_dat.columns):
        ax = axs[i]

        this_dat = dat.loc[:, [c, 'Outcome']]
        this_dat.loc[this_dat.loc[:, 'Outcome'] == 1, 'Outcome'] = 'Unfavourable'
        this_dat.loc[this_dat.loc[:, 'Outcome'] == 2, 'Outcome'] = 'Favourable'

        sns.boxplot(data=this_dat, x='Outcome', y=c, ax=ax, color='w')
        sns.swarmplot(data=this_dat, x='Outcome', y=c, ax=ax)

        # ax.scatter(scat[:, 0], scat[:, 1], c=scat[:, 2], cmap='RdBu')

        ax.xaxis.label.set_visible(False)
        ax.yaxis.label.set_visible(False)
        ax.set_title(c)
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
        # ax.set_xticklabels(['Fav', 'Unfav'], rotation=45)
        ylim = list(ax.get_ylim())
        if ylim[0] < 0:
            ylim[0] = 0.
            ax.set_ylim(ylim)
    # fig.subplots_adjust(left=0.025, right=0.99, wspace=0.4, bottom=0.2)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter.pdf'))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for i, c in enumerate(peaks_dat.columns):
        ttl = c.lower().replace(' ', '_')
        ax.cla()

        fav = dat.loc[dat.Outcome == 2, c].values
        unfav = dat.loc[dat.Outcome == 1, c].values

        fav_age = dat.loc[dat.Outcome == 2, "%s age" % c].values
        unfav_age = dat.loc[dat.Outcome == 1, "%s age" % c].values

        ax.scatter(fav_age, fav, c='b', label='Fav')
        ax.scatter(unfav_age, unfav, c='g', label='Unfav')
        ax.set_title(c)

        ax.set_xlabel("Age (hours)")
        ax.set_ylabel(c)

        fig.savefig(os.path.join(outdir, "scatter_marker_%s_vs_age.png" % ttl), dpi=200)
        fig.savefig(os.path.join(outdir, "scatter_marker_%s_vs_age.pdf" % ttl))


    # pair plot

    gs = plt.GridSpec(nvar - 1, nvar - 1)
    fig = plt.figure(figsize=(12, 12))
    ms = 10

    row_axs = []
    col_axs = []
    for i, row_lbl in enumerate(peaks_dat.columns):
        if i == 0:
            continue
        rax = fig.add_subplot(gs[i - 1, 0])
        rax.set_ylabel(row_lbl.replace(' peak', ''))
        row_axs.append(rax)
    for j, col_lbl in enumerate(peaks_dat.columns):
        if j == (nvar - 1):
            continue
        cax = fig.add_subplot(gs[nvar - 2, j])
        cax.set_xlabel(col_lbl.replace(' peak', ''))
        col_axs.append(cax)

    for i, row_lbl in enumerate(peaks_dat.columns):
        for j, col_lbl in enumerate(peaks_dat.columns):
            if i == 0:
                continue
            if j == (nvar - 1):
                continue
            cax = col_axs[j]
            rax = row_axs[i - 1]
            ax = None
            sharex = None
            sharey = None

            if i == (nvar - 1):
                ax = cax
                if j != 0:
                    ax.yaxis.set_visible(False)
            else:
                sharex = cax

            if j == 0:
                ax = rax
                if i != (nvar - 1):
                    ax.xaxis.set_visible(False)
            else:
                sharey = rax

            if ax is None:
                ax = fig.add_subplot(gs[i - 1, j], sharex=sharex, sharey=sharey)
                ax.xaxis.set_visible(False)
                ax.yaxis.set_visible(False)

            if j < i:
                x = peaks_dat.loc[:, col_lbl]
                y = peaks_dat.loc[:, row_lbl]
                z = outcomes
                ax.scatter(x[z == 2], y[z == 2], color='b', s=ms, alpha=0.6, label='Favourable')
                ax.scatter(x[z == 1], y[z == 1], color='g', s=ms, alpha=0.6, label='Unfavourable')
                ax.grid(False)
            else:
                ax.grid(False)
                ax.set_facecolor('none')

    gs.update(bottom=0.07, left=0.06, top=.99, right=.99, hspace=0.02, wspace=0.02)
    # FIXME: this need to be run manually, presumably due to some buffer not flushing until the plot is fully plotted?
    # common.align_labels(row_axs, axis='y')
    fig.savefig(os.path.join(outdir, "pairplot.png"), dpi=200)

    # based on this, it's worth plotting Plt peak vs Plt trough on its own
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(peaks_dat.loc[outcomes==2, 'Plt peak'], peaks_dat.loc[outcomes==2, 'Plt trough'], c='b', label='Favourable')
    ax.scatter(peaks_dat.loc[outcomes==1, 'Plt peak'], peaks_dat.loc[outcomes==1, 'Plt trough'], c='g', label='Unfavourable')
    ax.legend(loc='upper right')
    ax.set_xlabel("Plt peak")
    ax.set_ylabel("Plt trough")
    fig.savefig(os.path.join(outdir, "plt_peak_vs_trough.png"), dpi=200)

    # scatter plot with platelet transfusion included
    platelet_idx = dat.loc[:, 'Blood product'].str.contains('Platelets').fillna(False)
    X_platelet = X.copy()
    X_platelet.insert(0, 'platelets', platelet_idx.astype(int))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # platelets
    x = X_platelet.loc[platelet_idx, 'Plt trough']
    y = X_platelet.loc[platelet_idx, 'Plt peak']
    x_fav = x[outcomes.loc[platelet_idx] == 2]
    y_fav = y[outcomes.loc[platelet_idx] == 2]
    x_unfav = x[outcomes.loc[platelet_idx] == 1]
    y_unfav = y[outcomes.loc[platelet_idx] == 1]
    ax.scatter(x_fav, y_fav, facecolor='b', marker='s', label='Favourable, given platelets')
    ax.scatter(x_unfav, y_unfav, facecolor='g', marker='s', label='Unfavourable, given platelets')

    # not platelets
    x = X_platelet.loc[~platelet_idx, 'Plt trough']
    y = X_platelet.loc[~platelet_idx, 'Plt peak']
    x_fav = x[outcomes.loc[~platelet_idx] == 2]
    y_fav = y[outcomes.loc[~platelet_idx] == 2]
    x_unfav = x[outcomes.loc[~platelet_idx] == 1]
    y_unfav = y[outcomes.loc[~platelet_idx] == 1]
    ax.scatter(x_fav, y_fav, facecolor='b', marker='o', label='Favourable, not given platelets')
    ax.scatter(x_unfav, y_unfav, facecolor='g', marker='o', label='Unfavourable, not given platelets')

    # not platelet trend line
    lr = stats.linregress(x, y)
    x = np.array([X.loc[:, 'Plt trough'].min(), X.loc[:, 'Plt trough'].max()])
    y = lr.intercept + lr.slope * x
    ax.plot(x, y, 'k--')

    ax.legend(loc='upper right')
    ax.set_xlabel('Plt trough')
    ax.set_ylabel('Plt peak')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "scatter_plt_platelet_products.png"), dpi=200)

    # build a simple logistic regression model

    # add a constant and fit
    logit_model = sm.Logit(outcomes == 2, sm.add_constant(X))  # outcome is FAVOURABLE
    result = logit_model.fit()
    print "Predicting outcome using biomarkers:"
    print result.summary()
    ci = result.conf_int()
    ci.columns = ["2.5%", "97.5%"]
    ci.insert(0, "Odds ratio", result.params)
    ci = np.exp(ci)
    ci.insert(0, "P value", result.pvalues)

    # plot the decision boundary for model with Plt peak and Plt trough (only)
    # plt_dat = X.loc[:, ['Plt peak', 'Plt trough']]
    # logit_model_plt = sm.Logit(outcomes == 2, sm.add_constant(plt_dat))
    # result_plt = logit_model_plt.fit()
    # coeff = result_plt.params
    # intercept = -coeff['const'] / coeff['Plt peak']
    # slope = -coeff['Plt trough'] / coeff['Plt peak']
    # fit_x = np.linspace(
    #     peaks_dat.loc[:, 'Plt trough'].min(),
    #     peaks_dat.loc[:, 'Plt trough'].max(),
    #     20
    # )
    # fit_y = intercept + slope * fit_x
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.scatter(peaks_dat.loc[outcomes==2, 'Plt trough'], peaks_dat.loc[outcomes==2, 'Plt peak'], c='b', label='Favourable')
    # ax.scatter(peaks_dat.loc[outcomes==1, 'Plt trough'], peaks_dat.loc[outcomes==1, 'Plt peak'], c='g', label='Unfavourable')
    # ax.plot(fit_x, fit_y, 'r--', label='Model decision boundary')
    # ax.legend(loc='upper left')
    # ax.set_xlabel('Plt trough')
    # ax.set_ylabel('Plt peak')
    # fig.tight_layout()
    # fig.savefig(os.path.join(outdir, "logit_plt_peak_trough_decision_boundary.png"), dpi=200)

    # repeat logit but with platelet variable included

    logit_model = sm.Logit(outcomes == 2, sm.add_constant(X_platelet))  # outcome is FAVOURABLE
    result_platelet = logit_model.fit()
    print "Predicting outcome using biomarkers and platelets given:"
    print result_platelet.summary()
    ci_platelet = result_platelet.conf_int()
    ci_platelet.columns = ["2.5%", "97.5%"]
    ci_platelet.insert(0, "Odds ratio", result_platelet.params)
    ci_platelet = np.exp(ci_platelet)
    ci_platelet.insert(0, "P value", result_platelet.pvalues)
