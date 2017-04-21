from load_data import rnaseq_data
from microarray import process
from utils.output import unique_output_dir
import os
import pandas as pd
from scipy import stats
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


def extract_present_genes(x):
    data = x.loc[obj.data.index.str.contains('ENSG')]
    data = data.loc[data.any(axis=1), :]
    return data


def norm(x):
    return (x + 1).divide(x.sum(axis=0) + 1, axis=1)


if __name__ == '__main__':

    obj = rnaseq_data.gbm_ribozero_samples_loader(annotate_by='Ensembl Gene ID')

    # few plots about number assigned
    assgn = obj.data.loc[obj.data.index.str.contains('N_')]
    s = obj.data.sum(axis=0)
    s.name = 'N_unique_assignment'
    assgn = assgn.append(s).transpose()
    ax = assgn.plot.bar()
    ax.set_position([0.05, 0.17, 0.7, 0.8])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.setp(ax.xaxis.get_ticklabels(), rotation=45)

    data = obj.data.loc[obj.data.index.str.contains('ENSG')]
    data = data.loc[data.any(axis=1), :]
    data_n = norm(data)

    # extract the paired samples for further work
    a = np.log2(data_n.loc[:, 'GBM031 FF'])
    b = np.log2(data_n.loc[:, 'GBM031 FFPE'])

    # scatter plot of the paired samples
    lr = stats.linregress(a, b)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(a, b, label=None)

    ax.plot(
        [min(a), max(a)], [min(a) * lr.slope + lr.intercept, max(a) * lr.slope + lr.intercept],
        'k--',
        label='r = %.3f' % lr.rvalue
    )
    ax.legend(loc='upper left')
    ax.set_xlabel('GBM031 flash frozen expression')
    ax.set_ylabel('GBM031 FFPE expression')

    # try comparing to the cell line results
    # NB: this means Ribo-Zero vs poly(A)

    obj_polya = rnaseq_data.gbm_paired_samples_loader(annotate_by='Ensembl Gene ID', source='star')
    data_polya = obj_polya.data.loc[:, obj_polya.data.columns.str.contains('GBM')]

    data_all = pd.concat((obj.data, data_polya), axis=1)
    data_all = extract_present_genes(data_all)
    data_all_n = norm(data_all)
    data_all_yg = process.yugene_transform(data_all)
    data_all_log_yg = process.yugene_transform(np.log2(data_all_n))

    # generate a DataFrame with pairwise comparison
    samples = (
        '018', '019', '026', '031'
    )

    pearson_df = pd.DataFrame(columns=['a', 'b', 'slope', 'intercept', 'rvalue', 'pvalue', 'stderr'])
    pearson_matched = pearson_df.copy()

    def add_one_result(a, b, lr, df):
        d = {'a': a, 'b': b}
        d.update(lr.__dict__)
        s = pd.Series(d)
        return df.append(s, ignore_index=True)

    # define a variable here that allows us to choose which data norming method we use
    comp_data = np.log2(data_all_n)

    for i in samples:
        b1 = False
        b2 = False
        t = 'GBM%s FFPE' % i
        if t in data_all_yg.columns:
            b1 = True
            lr = stats.linregress(
                comp_data.loc[:, t],
                comp_data.loc[:, 'GBM%s' % i],
            )
            pearson_df = add_one_result(t, 'GBM%s' % i, lr, pearson_df)

        t = 'GBM%s FF' % i
        if t in data_all_yg.columns:
            b2 = True
            lr = stats.linregress(
                comp_data.loc[:, t],
                comp_data.loc[:, 'GBM%s' % i],
            )
            pearson_df = add_one_result(t, 'GBM%s' % i, lr, pearson_df)

        if b1 and b2:
            lr = stats.linregress(
                comp_data.loc[:, 'GBM%s FFPE' % i],
                comp_data.loc[:, 'GBM%s FF' % i],
            )
            pearson_matched = add_one_result('GBM%s FFPE' % i, 'GBM%s FF' % i, lr, pearson_matched)

    pearson_df.columns = ['ribozero', 'cell_line', 'slope', 'intercept', 'rvalue', 'pvalue', 'stderr']
    pearson_matched.columns = ['ribozero_1', 'ribozero_2', 'slope', 'intercept', 'rvalue', 'pvalue', 'stderr']

