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
from plotting import common

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

    # TODO: scatterplots of peak/trough vs age for each variable, coloured by outcome

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

    # build a simple logistic regression model
    from statsmodels import api as sm
    from sklearn.preprocessing import Imputer
    # fill nan with median for that variable
    X = peaks_dat.copy()
    imp = Imputer(missing_values='NaN', strategy='median', axis=0)
    X = imp.fit_transform(X)
    X = pd.DataFrame(X, index=peaks_dat.index, columns=peaks_dat.columns)

    # add a constant and fit
    logit_model = sm.Logit(outcomes == 2, sm.add_constant(X))  # outcome is FAVOURABLE
    result = logit_model.fit()
    print "statsmodels.Logit model fitted:"
    print result.summary()
    ci = result.conf_int()
    ci.columns = ["2.5%", "97.5%"]
    ci.insert(0, "Odds ratio", result.params)
    ci = np.exp(ci)
    ci.insert(0, "P value", result.pvalues)

    # plot the decision boundary for model with Plt peak and Plt trough (only)
    plt_dat = X.loc[:, ['Plt peak', 'Plt trough']]
    logit_model_plt = sm.Logit(outcomes == 2, sm.add_constant(plt_dat))
    result_plt = logit_model_plt.fit()
    coeff = result_plt.params
    intercept = -coeff['const'] / coeff['Plt peak']
    slope = -coeff['Plt trough'] / coeff['Plt peak']
    fit_x = np.linspace(
        peaks_dat.loc[:, 'Plt trough'].min(),
        peaks_dat.loc[:, 'Plt trough'].max(),
        20
    )
    fit_y = intercept + slope * fit_x

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(peaks_dat.loc[outcomes==2, 'Plt trough'], peaks_dat.loc[outcomes==2, 'Plt peak'], c='b', label='Favourable')
    ax.scatter(peaks_dat.loc[outcomes==1, 'Plt trough'], peaks_dat.loc[outcomes==1, 'Plt peak'], c='g', label='Unfavourable')
    ax.plot(fit_x, fit_y, 'r--', label='Model decision boundary')
    ax.legend(loc='upper left')
    ax.set_xlabel('Plt trough')
    ax.set_ylabel('Plt peak')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "logit_plt_peak_trough_decision_boundary.png"), dpi=200)



    # TODO: add statistical significance of difference between groups in boxplot

    # PCA

    X = peaks_dat.copy()

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