from load_data import rnaseq_data
from microarray import process
from stats import transformations
from utils.output import unique_output_dir
import os
import pandas as pd
from scipy import stats
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from plotting.pca import pca_plot_by_group_2d
from plotting import clustering
from sklearn.decomposition import PCA


def extract_present_genes(x):
    data = x.loc[obj.data.index.str.contains('ENSG')]
    data = data.loc[data.any(axis=1), :]
    return data


def norm(x):
    return (x + 1).divide(x.sum(axis=0) + 1, axis=1)


def get_assigned_proportions(dat):
    assgn = dat.loc[dat.index.str.contains('N_')]
    s = dat.sum(axis=0)
    s.name = 'N_unique_assignment'
    assgn = assgn.append(s).transpose()
    return assgn


def compute_pairwise_corr(cell_line_dat, tissue_dat, n_genes=None):
    """
    Used
    :param _dat: Two pd DataFrames containing cell line and tissue data
    :param n_genes: If supplied, this is an integer, Genes are ranked by MAD
    :return:
    """
    if n_genes is not None:
        mad = transformations.median_absolute_deviation(
            pd.concat((cell_line_dat, tissue_dat), axis=1)
        ).sort_values(ascending=False)
        g = mad.index[:n_genes]
        cell_line_dat = cell_line_dat.loc[g]
        tissue_dat = tissue_dat.loc[g]

    return cell_line_dat.apply(lambda x: tissue_dat.corrwith(x), axis=0)


if __name__ == '__main__':
    outdir = unique_output_dir('ribozero_vs_polya')

    obj = rnaseq_data.gbm_ribozero_samples_loader(annotate_by='Ensembl Gene ID')

    data = obj.data.loc[obj.data.index.str.contains('ENSG')]
    data = data.loc[data.any(axis=1), :]
    data_n = norm(data)

    # load relevant poly(A) data, too

    obj_polya = rnaseq_data.gbm_paired_samples_loader(annotate_by='Ensembl Gene ID', source='star')
    data_polya = obj_polya.data.loc[:, obj_polya.data.columns.str.contains('GBM')]

    data_all = pd.concat((obj.data, data_polya), axis=1)
    data_all = extract_present_genes(data_all)
    data_all_n = norm(data_all)
    data_all_yg = process.yugene_transform(data_all, resolve_ties=False)
    data_all_log_yg = process.yugene_transform(np.log2(data_all_n), resolve_ties=False)


    # number assigned bar chart

    assgn = get_assigned_proportions(obj.data)
    assgn_polya = get_assigned_proportions(obj_polya.data).loc['GBM031']
    assgn = assgn.append(assgn_polya)

    ax = assgn.plot.bar()
    ax.set_position([0.05, 0.17, 0.7, 0.8])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.setp(ax.xaxis.get_ticklabels(), rotation=45)
    ax.set_ylabel('# reads')

    ax.figure.savefig(os.path.join(outdir, 'assignment_proportions_bar.png'), dpi=200)
    ax.figure.savefig(os.path.join(outdir, 'assignment_proportions_bar.pdf'))


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

    fig.savefig(os.path.join(outdir, 'gbm031_matched_scatter.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'gbm031_matched_scatter.pdf'))

    # try comparing to the cell line results
    # NB: this means Ribo-Zero vs poly(A)

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

    # get correlation between (FFPE/frozen) and (primary culture)
    tissue_dat = comp_data.iloc[:, comp_data.columns.str.contains('FF')]
    cell_line_dat = comp_data.iloc[:, ~comp_data.columns.str.contains('FF')]

    # use all genes
    corrmat = compute_pairwise_corr(cell_line_dat, tissue_dat)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = sns.heatmap(
        corrmat,
        square=True,
        annot=True,
        fmt='.3f',
        cmap='RdBu_r',
        ax=ax
    )
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "a_vs_b_pairwise_correlation_allgenes.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "a_vs_b_pairwise_correlation_allgenes.pdf"))

    # use top 1500 genes
    corrmat = compute_pairwise_corr(cell_line_dat, tissue_dat, n_genes=1500)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = sns.heatmap(
        corrmat,
        square=True,
        annot=True,
        fmt='.3f',
        cmap='RdBu_r',
        ax=ax
    )
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "a_vs_b_pairwise_correlation_1500genes.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "a_vs_b_pairwise_correlation_1500genes.pdf"))


    # PCA

    # define a variable here that allows us to choose which data norming method we use
    comp_data = np.log2(data_all_n)

    # by library type
    pca = PCA(n_components=8)
    pca.fit(comp_data.transpose())
    y = pca.transform(comp_data.transpose())
    subgroups = pd.Series(['Ribo-Zero FFPE'] * 3 + ['Ribo-Zero frozen'] * 2 + ['Poly(A)'] * 5)
    ax = pca_plot_by_group_2d(y, subgroups=subgroups, ellipses=False, auto_scale=False)
    ax.legend(loc='upper left', frameon=True, facecolor='w', edgecolor='b')
    plt.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "pca_ribozero_polya.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "pca_ribozero_polya.pdf"))

    # by cell line
    subgroups = data_all_n.columns.str.replace(r'[^0-9]*', '')
    ax = pca_plot_by_group_2d(y, subgroups=subgroups, ellipses=False, auto_scale=False)
    ax.legend(loc='upper left', frameon=True, facecolor='w', edgecolor='b')
    plt.tight_layout()
    ax.figure.savefig(os.path.join(outdir, "pca_ribozero_polya_byline.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "pca_ribozero_polya_byline.pdf"))

    # hierarchical clustering
    subgroups = pd.DataFrame(
        ['b'] * 3 + ['r'] * 2 + ['g'] * 5,
        index=data_all_n.columns,
        columns=['Prep type']
    )
    legend_lbl = {'FFPE': 'b', 'frozen': 'r', 'Poly(A)': 'g'}
    res = clustering.dendrogram_with_colours(
        comp_data,
        subgroups,
        legend_labels=legend_lbl,
        metric='correlation',
        method='average'
    )
    res['fig'].savefig(os.path.join(outdir, 'clustering_dendrogram.png'), dpi=200)
    res['fig'].savefig(os.path.join(outdir, 'clustering_dendrogram.pdf'))

    # can add n_gene kwarg here to pick top N genes by MAD:
    cg = clustering.plot_correlation_clustermap(comp_data)
    cg.savefig(os.path.join(outdir, 'clustering_corr_map.png'), dpi=200)
    cg.savefig(os.path.join(outdir, 'clustering_corr_map.pdf'))