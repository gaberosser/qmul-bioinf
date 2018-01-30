import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from scipy import stats
from statsmodels import api as sm

from settings import GIT_LFS_DATA_DIR
from utils.output import unique_output_dir

BIOMARKERS = [
    'Hb',
    'Plt',
    'Neutrophil',
    'Lymphocyte',
    'PT',  # PT test and INR are highly related
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
    # fn = os.path.join(DATA_DIR, 'divyen_shah', 'cleaned_data_nov_2016.csv')
    fn = os.path.join(GIT_LFS_DATA_DIR, 'divyen_shah', 'cleaned_data_jan_2018.csv')
    return pd.read_csv(fn, header=0, na_values='-', index_col=0)


def standardise(data, axis=0):
    if axis == 0:
        return data.subtract(data.mean(axis=0), axis=1).divide(data.std(axis=0), axis=1)
    elif axis == 1:
        return data.subtract(data.mean(axis=1), axis=0).divide(data.std(axis=1), axis=0)
    else:
        raise AttributeError("Axis must be 0 (norm by col) or 1 (norm by row)")


if __name__ == "__main__":
    outdir = unique_output_dir("hie_results", reuse_empty=True)

    dat = load_cleaned_data()
    dat.loc[:, 'batch'] = [t[:2] for t in dat.index]
    biomarkers = dat.loc[:, (
        BIOMARKER_PEAK_COLS
        + BIOMARKER_TROUGH_COLS
        + BIOMARKER_PEAK_AGE_COLS
        + BIOMARKER_TROUGH_AGE_COLS
    )]
    outcomes = dat.loc[:, OUTCOME_COL]

    ## 1) Correlation between variables
    peaks_dat = dat.loc[:, BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS]
    corr = peaks_dat.corr(method='spearman')
    fig = plt.figure(figsize=(6.7, 5.5))
    ax = fig.add_subplot(111)
    sns.heatmap(corr, ax=ax, cmap='RdBu_r')
    ax.set_aspect('equal')
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "spearman_correlation_plot.png"), dpi=200)

    # CONCLUSION: PT and INR peaks are highly correlated -> remove PT
    peaks_dat = peaks_dat.loc[:, peaks_dat.columns != 'PT peak']
    nvar = peaks_dat.shape[1]

    # standardize - necessary when using array plots to keep the range the same
    peaks_dat = dat.loc[:, BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS]
    bm_peaks = standardise(peaks_dat)
    # bm_peaks = peaks_dat.subtract(peaks_dat.mean(axis=0), axis=1).divide(peaks_dat.std(axis=0), axis=1)
    bm_peaks = pd.concat((bm_peaks, outcomes), axis=1)
    bm_peaks.columns = bm_peaks.columns.str.replace(' peak', '')

    ## TODO: pair plot and pairwise correlation
    ax = sns.heatmap(peaks_dat.corr())
    # FIXME: rotate axes ticklabels and tidy up

    # pairplot by batch
    xx = pd.concat((peaks_dat, dat.loc[:, 'batch']), axis=1)

    # Look for batch effects using PCA on biomarker values
    # this requires us to fill NA values with column mean

    data_for_pca = peaks_dat.copy()

    # use this code to fill missing values with the mean:
    # for col, colval in data_for_pca.iteritems():
    #     if colval.isnull().any():
    #         data_for_pca.loc[:, col] = colval.fillna(colval.mean())

    # use this code to discard any variables with missing data:
    data_for_pca = peaks_dat.loc[:, peaks_dat.isnull().sum(axis=0) == 0]

    pca = PCA(n_components=6)
    pca.fit(data_for_pca)
    pp = pca.transform(data_for_pca)
    # scatter, shading by component
    fig = plt.figure()
    ax = fig.add_subplot(111)

    batches = dat.batch.factorize()
    colours = ['r', 'b', 'g', 'y', 'k']
    for i, bt in enumerate(batches[1]):
        idx = batches[0] == i
        ax.scatter(pp[idx, 0], pp[idx, 1], c=colours[i], label=bt)
    ax.legend(loc='best', fontsize=14, frameon=True, facecolor='w', fancybox=True)
    ax.set_xlabel("Principle component 1")
    ax.set_ylabel("Principle component 2")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "pca_by_batch.pdf"))
    fig.savefig(os.path.join(outdir, "pca_by_batch.png"), dpi=200)

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


    # does peak/trough value correlate with age?
    val_age_corr = {
        'all': {},
        'fav': {},
        'unfav': {}
    }
    p_threshold = 0.05
    fig = plt.figure()
    ax = fig.add_subplot(111)

    cols = BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS

    for i, c in enumerate(cols):
        age_col = "%s age" % c

        all_values = dat.loc[:, [c, age_col]].dropna().values
        fav = dat.loc[dat.Outcome == 2, [c, age_col]].dropna().values
        unfav = dat.loc[dat.Outcome == 1, [c, age_col]].dropna().values

        val_age_corr['all'][c] = stats.linregress(*all_values.transpose())
        val_age_corr['fav'][c] = stats.linregress(*fav.transpose())
        val_age_corr['unfav'][c] = stats.linregress(*unfav.transpose())

        # plot values
        colour = 'r' if val_age_corr['all'][c].pvalue < p_threshold else 'gray'
        ax.plot(i, val_age_corr['all'][c].slope, ls='none', marker='o', c=colour, label='All' if i==0 else None)
        colour = 'r' if val_age_corr['fav'][c].pvalue < p_threshold else 'gray'
        ax.plot(i, val_age_corr['fav'][c].slope, ls='none', marker='s', c=colour, label='Favourable' if i==0 else None)
        colour = 'r' if val_age_corr['unfav'][c].pvalue < p_threshold else 'gray'
        ax.plot(i, val_age_corr['unfav'][c].slope, ls='none', marker='^', c=colour, label='Unfavourable' if i==0 else None)

    ax.legend()
    ax.set_ylabel("Regression slope")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'linregress_value_with_age.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'linregress_value_with_age.pdf'))

    if False:
        # we need to set the range so that the hist function doesn't complain about NaNs
        hist_kws = {'range': (-2, 6)}
        g = sns.pairplot(bm_peaks, hue='Outcome', diag_kws=hist_kws)
        for ax in g.fig.get_axes():
            ax.set_xticks([])
            ax.set_yticks([])
        g.fig.subplots_adjust(left=0.04, bottom=0.04)
        g.fig.legends[0].set_visible(False)


    from scripts.mb_subgroup_classifier.shrunken_centroids import run_validation

    deltas = np.linspace(0., 1.5, 30)
    # use all for training
    training_data = bm_peaks.loc[:, bm_peaks.columns.difference(['Outcome'])].transpose()
    train_err = run_validation(deltas, training_data, bm_peaks.Outcome)

    # box and whisker for each variable by outcome
    bwf = []
    bwu = []
    ttls = []
    for c in BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS:
        bwf.append(dat.loc[dat.loc[:, OUTCOME_COL] == 2, c].dropna())  # favourable outcome
        bwu.append(dat.loc[dat.loc[:, OUTCOME_COL] == 1, c].dropna())  # unfavourable outcome
        ttls.append(c.replace(' peak', ''))

    jitter = 0.2
    if nvar % 2 == 0:
        fig, axs = plt.subplots(nrows=2, ncols=nvar / 2, sharex=True, figsize=(12, 8))
        axs = axs.flat
    else:
        fig, axs = plt.subplots(nrows=1, ncols=nvar, sharex=True, figsize=(15, 5))

    for i, c in enumerate(BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS):
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
    for i, c in enumerate(BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS):
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

    # TODO: scatterplots of peak/trough vs age for each variable, coloured by outcome
    # TODO: add statistical significance of difference between groups in boxplot

    # PCA

    X = peaks_dat.copy()
    from sklearn.preprocessing import Imputer
    imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
    imp.fit(X)
    X = imp.transform(X)

    pca = PCA(n_components=10)
    pca.fit(X)
    y = pca.transform(X)

    fav_idx = np.where(dat.Outcome==2)[0]
    unfav_idx = np.where(dat.Outcome==1)[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(y[fav_idx, 0], y[fav_idx, 1], c='b')
    ax.scatter(y[unfav_idx, 0], y[unfav_idx, 1], c='g')

    # statsmodels logistic regression
    try:
        peaks_dat_const = sm.add_constant(peaks_dat)
        est = sm.OLS(outcomes, peaks_dat_const)
        est2 = est.fit()
        print est2.summary()
    except Exception as exc:
        print "Failed to fit statsmodels logistic regression model"
        print repr(exc)

    # decision tree classifier
    from sklearn import tree
    import graphviz
    clf = tree.DecisionTreeClassifier()
    clf = clf.fit(peaks_dat, outcomes)
    dot_data = tree.export_graphviz(clf, out_file=None, feature_names=peaks_dat.columns, class_names=['Unfav', 'Fav'],
                                    filled=True, rounded=True)
    graph = graphviz.Source(dot_data)
    graph.render(os.path.join(outdir, 'decision_tree_graph'))# TODO


    # decision tree regression

    from sklearn.linear_model import LogisticRegression
    logistic = LogisticRegression(C=1)  # what is C?
    inferred = []
    true = []

    for i in range(1000):

        idx_shuffle = np.random.permutation(X.shape[0])
        cutoff = int(X.shape[0] * 0.9)
        idx_train = idx_shuffle[:cutoff]
        idx_test = idx_shuffle[cutoff:]
        X_train = X[idx_train, :]
        y_train = dat.Outcome.values[idx_train]
        X_test = X[idx_test, :]
        y_test = dat.Outcome.values[idx_test]

        logistic.fit(X_train, y_train)
        y_test_inf = logistic.predict(X_test)

        inferred.append(y_test_inf)
        true.append(y_test)

    # performance
    fp = tp = fn = tn = 0
    n = 0
    for inf, tr in zip(inferred, true):
        tn += ((inf == 2) & (tr == 2)).sum()
        tp += ((inf == 1) & (tr == 1)).sum()
        fn += ((inf == 2) & (tr == 1)).sum()
        fp += ((inf == 1) & (tr == 2)).sum()
    n = tn + tp + fn + fp

    # now fit to the full dataset and analyse the coefficients
    logistic.fit(X, dat.Outcome.values)
    odds_coef = np.exp(pd.Series(logistic.coef_.flatten(), index=peaks_dat.columns))