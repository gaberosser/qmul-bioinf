import pandas as pd
import os
from settings import DATA_DIR
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from scripts.output import unique_output_dir
from sklearn.decomposition import PCA


BIOMARKERS = [
    'Hb',
    'Plt',
    'Neutrophil',
    'Lymphocyte',
    'PT',
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


if __name__ == "__main__":
    outdir = unique_output_dir("hie_results", reuse_empty=True)

    fn = os.path.join(DATA_DIR, 'divyen_shah', 'cleaned_data_nov_2016.csv')
    dat = pd.read_csv(fn, header=0, na_values='-', index_col=0)
    biomarkers = dat.loc[:, (
        BIOMARKER_PEAK_COLS
        + BIOMARKER_TROUGH_COLS
        + BIOMARKER_PEAK_AGE_COLS
        + BIOMARKER_TROUGH_AGE_COLS
        + [OUTCOME_COL]
    )]

    # standardize - necessary when using array plots to keep the range the same
    peaks_dat = dat.loc[:, BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS]
    bm_peaks = peaks_dat.subtract(peaks_dat.mean(axis=0), axis=1).divide(peaks_dat.std(axis=0), axis=1)
    bm_peaks = pd.concat((bm_peaks, dat.loc[:, OUTCOME_COL]), axis=1)
    bm_peaks.columns = bm_peaks.columns.str.replace(' peak', '')

    if False:
        # we need to set the range so that the hist function doesn't complain about NaNs
        hist_kws = {'range': (-2, 6)}
        g = sns.pairplot(bm_peaks, hue='Outcome', diag_kws=hist_kws)
        for ax in g.fig.get_axes():
            ax.set_xticks([])
            ax.set_yticks([])
        g.fig.subplots_adjust(left=0.04, bottom=0.04)
        g.fig.legends[0].set_visible(False)


    from scripts.mb_subgroup_classifier.shrunken_centroids import run_validation, NearestCentroidClassifier
    deltas = np.linspace(0., 1.5, 30)
    # use all for training
    training_data = bm_peaks.loc[:, bm_peaks.columns.difference(['Outcome'])].transpose()
    train_err = run_validation(deltas, training_data, bm_peaks.Outcome)

    # box and whisker for each variable by outcome
    bwf = []
    bwu = []
    ttls = []
    for c in BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS:
        bwf.append(dat.loc[dat.loc[:, OUTCOME_COL] == 'favourable', c].dropna())
        bwu.append(dat.loc[dat.loc[:, OUTCOME_COL] == 'unfavourable', c].dropna())
        ttls.append(c.replace(' peak', ''))

    # set positions
    xf = range(1, len(bwf) * 3, 3)
    xu = range(2, len(bwf) * 3, 3)

    fig, axs = plt.subplots(nrows=1, ncols=len(BIOMARKER_PEAK_COLS) + len(BIOMARKER_TROUGH_COLS), sharex=True, figsize=(15, 5))
    for i in range(len(bwf)):
        ax = axs[i]
        ax.boxplot([bwf[i], bwu[i]], widths=0.5)
        ax.set_xticklabels(['Fav', 'Unfav'], rotation=45)
        ax.set_title(ttls[i])

    fig.subplots_adjust(left=0.05, right=0.99)
    fig.savefig(os.path.join(outdir, 'box_whisker_basic.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'box_whisker_basic.pdf'))

    jitter = 0.2
    fig, axs = plt.subplots(nrows=1, ncols=len(BIOMARKER_PEAK_COLS) + len(BIOMARKER_TROUGH_COLS), sharex=True,
                            figsize=(15, 5))
    for i, c in enumerate(BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS):
        ax = axs[i]

        # # extract data for scatterplot
        # fav = dat.loc[dat.Outcome == 'favourable', c].values
        # unfav = dat.loc[dat.Outcome == 'unfavourable', c].values
        # scat = np.zeros((len(fav) + len(unfav), 3))
        #
        # # add peak age
        # fav_age = dat.loc[dat.Outcome == 'favourable', "%s age" % c].values
        # unfav_age = dat.loc[dat.Outcome == 'unfavourable', "%s age" % c].values
        #
        #
        # scat[:len(fav), 0] = np.random.normal(scale=jitter, size=len(fav))
        # scat[:len(fav), 1] = fav
        # scat[:len(fav), 2] = fav_age
        # scat[len(fav):, 0] = np.random.normal(loc=1, scale=jitter, size=len(unfav))
        # scat[len(fav):, 1] = unfav
        # scat[len(fav):, 2] = unfav_age

        sns.boxplot(data=dat.loc[:, [c, 'Outcome']], x='Outcome', y=c, ax=ax)
        sns.swarmplot(data=dat.loc[:, [c, 'Outcome']], x='Outcome', y=c, ax=ax)

        # ax.scatter(scat[:, 0], scat[:, 1], c=scat[:, 2], cmap='RdBu')

        ax.xaxis.label.set_visible(False)
        ax.yaxis.label.set_visible(False)
        ax.set_title(c)
        ax.set_xticklabels(['Fav', 'Unfav'], rotation=45)
        ylim = list(ax.get_ylim())
        if ylim[0] < 0:
            ylim[0] = 0.
            ax.set_ylim(ylim)
    fig.subplots_adjust(left=0.02, right=0.99, wspace=0.4)
    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'box_whisker_plus_scatter.pdf'))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, c in enumerate(BIOMARKER_PEAK_COLS + BIOMARKER_TROUGH_COLS):
        ttl = c.lower().replace(' ', '_')
        ax.cla()

        fav = dat.loc[dat.Outcome == 'favourable', c].values
        unfav = dat.loc[dat.Outcome == 'unfavourable', c].values

        fav_age = dat.loc[dat.Outcome == 'favourable', "%s age" % c].values
        unfav_age = dat.loc[dat.Outcome == 'unfavourable', "%s age" % c].values

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

    fav_idx = np.where(dat.Outcome=='favourable')[0]
    unfav_idx = np.where(dat.Outcome=='unfavourable')[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(y[fav_idx, 0], y[fav_idx, 1], c='b')
    ax.scatter(y[unfav_idx, 0], y[unfav_idx, 1], c='g')

    from sklearn.linear_model import LogisticRegression
    logistic = LogisticRegression(C=1e5)  # what is C?
    inferred = []

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
