from load_data import rnaseq_data
from plotting import corr, clustering
from stats import transformations
from matplotlib import pyplot as plt
from utils.output import unique_output_dir
import os
import numpy as np
import pandas as pd
import seaborn as sns


if __name__ == "__main__":
    outdir = unique_output_dir('mouse_validation')

    obj = rnaseq_data.mouse_nsc_validation_samples(annotate_by='Ensembl Gene ID')
    data = obj.data.loc[obj.data.index.str.contains('ENS')]
    meta = obj.meta

    # reorder for plotting niceness
    samples = ['eNSC%dmed' % i for i in (3, 5, 6)] \
              + ['eNSC%dmouse' % i for i in (3, 5, 6)] \
              + ['mDura%smouse' % i for i in ('3N1', '5N24A', '6N6')] \
              + ['mDura%shuman' % i for i in ('3N1', '5N24A', '6N6')]

    meta = meta.loc[samples]
    data = data.loc[:, meta.index]

    cpm = data.divide(meta.loc[:, 'read_count'].values, axis=1) * 1e6
    keep = (cpm > .5).sum(axis=1) > 5

    the_dat = np.log2(data.loc[keep] + 1)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax = corr.plot_correlation_coefficient_array(the_dat, vmin=0.4, ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'corr_coeff.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'corr_coeff.pdf'))

    cg = clustering.plot_correlation_clustermap(the_dat)
    cg.fig.savefig(os.path.join(outdir, 'corr_clustermap_all_genes.png'), dpi=200)
    cg.fig.savefig(os.path.join(outdir, 'corr_clustermap_all_genes.pdf'))

    for ng in [2500, 2000, 1500, 1000]:
        cg = clustering.plot_correlation_clustermap(the_dat, n_gene=ng)
        cg.fig.savefig(os.path.join(outdir, 'corr_clustermap_%d_genes.png' % ng), dpi=200)
        cg.fig.savefig(os.path.join(outdir, 'corr_clustermap_%d_genes.pdf' % ng))

    # compare within mice
    groups = [
        ['eNSC3med', 'eNSC3mouse', 'mDura3N1mouse', 'mDura3N1human'],
        ['eNSC5med', 'eNSC5mouse', 'mDura5N24Amouse', 'mDura5N24Ahuman'],
        ['eNSC6med', 'eNSC6mouse', 'mDura6N6mouse', 'mDura6N6human'],
    ]

    # create a dataframe for bar plotting
    within_cols = ["mouse", "correlation", "comparison"]
    the_corr = the_dat.corr()
    within_dat = []
    for i, g in zip([3, 5, 6], groups):
        within_dat.append([i, the_corr.loc[g[0], g[1]], "eNSCmed - eNSCmouse"])
        within_dat.append([i, the_corr.loc[g[0], g[2]], "eNSCmed - iNSCmouse"])
        within_dat.append([i, the_corr.loc[g[0], g[3]], "eNSCmed - iNSChuman"])
    within_mouse = pd.DataFrame(within_dat, columns=within_cols)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = sns.barplot(x='comparison', y='correlation', hue='mouse', data=within_mouse, ax=ax)
    ax.set_ylabel('Correlation')

    means = within_mouse.groupby('comparison').mean().loc[:, 'correlation']
    ax.plot([-.5, .5], [means.loc['eNSCmed - eNSCmouse']] * 2, 'k-')
    ax.plot([.5, 1.5], [means.loc['eNSCmed - iNSCmouse']] * 2, 'k-')
    ax.plot([1.5, 2.5], [means.loc['eNSCmed - iNSChuman']] * 2, 'k-')

    fig.savefig(os.path.join(outdir, "correlation_within_mouse.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "correlation_within_mouse.pdf"))

    def pwise_corr_boxplot(dat, n_gene=None):
        # all pairwise comparisons
        if n_gene is not None:
            mad = transformations.median_absolute_deviation(dat).sort_values(ascending=False)
            genes = mad.index[:n_gene]
            dat = dat.loc[genes]

        the_corr = dat.corr()
        pwise_corr = []
        for s in [r'eNSC.mouse', r'mDura.*mouse', r'mDura.*human']:
            ii = the_corr.index.str.contains(r'eNSC.med')
            jj = the_corr.columns.str.contains(s)
            pwise_corr.append(the_corr.loc[ii, jj].values.flatten())

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.boxplot(pwise_corr)
        ax.set_xticklabels(['eNSCmed - eNSCmouse', 'eNSCmed - iNSCmouse', 'eNSCmed - iNSChuman'])
        ax.set_ylabel('pairwise correlation')

        return ax

    ax = pwise_corr_boxplot(the_dat)
    ax.figure.savefig(os.path.join(outdir, "pairwise_correlation.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "pairwise_correlation.pdf"))