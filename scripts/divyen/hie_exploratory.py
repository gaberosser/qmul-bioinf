import pandas as pd
import os
from settings import DATA_DIR
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

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
    plt.show()

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


    # TODO: scatterplots of peak/trough vs age for each variable, coloured by outcome
    # TODO: add statistical significance of difference between groups in boxplot