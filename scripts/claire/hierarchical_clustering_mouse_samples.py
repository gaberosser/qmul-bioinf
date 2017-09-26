from scipy import stats
from load_data import rnaseq_data
from plotting import corr, clustering
from stats import transformations
from matplotlib import pyplot as plt
from utils.output import unique_output_dir
import os
import numpy as np
import pandas as pd
import seaborn as sns
import references


if __name__ == "__main__":
    outdir = unique_output_dir('mouse_validation')

    # load mouse data

    obj = rnaseq_data.mouse_nsc_validation_samples(annotate_by='Ensembl Gene ID')

    # reorder for plotting niceness

    samples = ['eNSC%dmed' % i for i in (3, 5, 6)] \
              + ['eNSC%dmouse' % i for i in (3, 5, 6)] \
              + ['mDura%smouse' % i for i in ('3N1', '5N24A', '6N6')] \
              + ['mDura%shuman' % i for i in ('3N1', '5N24A', '6N6')]

    obj.meta = obj.meta.loc[samples]
    obj.data = obj.data.loc[:, obj.meta.index]

    # load reference samples

    obj52564 = rnaseq_data.gse52564(
        annotate_by='Ensembl Gene ID',
        samples=[
            'Astrocyte1',
            'Astrocyte2',
            'Neuron1',
            'Neuron2',
            'OPC1',
            'OPC2',
        ]
    )
    obj43916 = rnaseq_data.gse43916(annotate_by='Ensembl Gene ID', samples=['NSCs'])
    obj86248 = rnaseq_data.gse86248(annotate_by='Ensembl Gene ID')
    obj36114 = rnaseq_data.gse36114(annotate_by='Ensembl Gene ID')

    obj_all = rnaseq_data.MultipleBatchLoader([
        obj,
        obj52564,
        obj43916,
        obj86248,
        obj36114
    ])

    data_all = obj_all.data.loc[obj_all.data.index.str.contains('ENS')]
    meta_all = obj_all.meta

    cpm_all = data_all.divide(data_all.sum(axis=0), axis=1) * 1e6
    keep = (cpm_all > .5).sum(axis=1) > 5

    the_dat = np.log2(data_all.loc[keep] + 1)

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

    # for ng in [2500, 2000, 1500, 1000]:
    #     cg = clustering.plot_correlation_clustermap(the_dat, n_gene=ng)
    #     cg.fig.savefig(os.path.join(outdir, 'corr_clustermap_%d_genes.png' % ng), dpi=200)
    #     cg.fig.savefig(os.path.join(outdir, 'corr_clustermap_%d_genes.pdf' % ng))#

    # mean and stdev for mouse X vs reference NSC
    groups_by_protocol = [
        ['eNSCmed', samples[:3],],
        ['eNSCmouse', samples[3:6],],
        ['mDuramouse', samples[6:9],],
        ['mDurahuman', samples[9:12],],
    ]
    lbl = [t[0] for t in groups_by_protocol]
    aa = the_dat.corrwith(the_dat.loc[:, 'NSCs'])
    val_dat = [aa.loc[t[1]] for t in groups_by_protocol]
    val_mean = [aa.loc[t[1]].mean() for t in groups_by_protocol]
    val_std = [aa.loc[t[1]].std() for t in groups_by_protocol]

    fig = plt.figure(figsize=(3.4, 5.5))
    ax = fig.add_subplot(111)
    ax.bar(range(len(groups_by_protocol)), val_mean, 0.7)
    for i, m in enumerate(val_mean):
        ax.plot([i, i], [m - val_std[i], m + val_std[i]], 'k-')
    ax.set_ylabel("Correlation with reference NSC")
    ax.set_xticks(range(len(groups_by_protocol)))
    ax.set_xticklabels(lbl, rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'corr_with_reference_nsc.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'corr_with_reference_nsc.pdf'))

    # can run a regular T test as follows
    # from scipy import stats
    # stats.ttest_rel(aa.loc[groups_by_protocol[2][1]].values, aa.loc[groups_by_protocol[3][1]].values)

    # same plot but this time we compare within-mouse to eNSCmed

    # which genes are responsible for the observed difference between iNSC and eNSC?
    the_ensc = the_dat.loc[:, groups_by_protocol[0][1]]
    the_insc = the_dat.loc[:, groups_by_protocol[2][1] + groups_by_protocol[3][1]]
    t_values = pd.Series(index=the_ensc.index)
    p_values = pd.Series(index=the_ensc.index)
    log2_mfc = pd.Series(index=the_ensc.index)

    for g in the_ensc.index:
        # T statistic
        tmp = stats.ttest_ind(the_insc.loc[g].values, the_ensc.loc[g].values)
        the_stat = tmp.statistic
        if np.abs(the_stat) > 500:
            the_stat = np.sign(the_stat) * 500
        t_values.loc[g] = the_stat
        p_values.loc[g] = tmp.pvalue
        # median fold-change (iNSC rel to eNSC)
        mfc = np.log2(np.median(the_insc.loc[g].values) / (np.median(the_ensc.loc[g].values) + 1e-6))
        if np.abs(mfc) > 20:
            mfc = np.sign(mfc) * 20
        log2_mfc.loc[g] = mfc

    t_values.dropna(inplace=True)
    log2_mfc.dropna(inplace=True)
    idx = t_values.index.intersection(log2_mfc.index)
    t_values = t_values.loc[idx]
    log2_mfc = log2_mfc.loc[idx]
    p_values = p_values.loc[idx]

    from statsmodels.sandbox.stats import multicomp
    tmp = multicomp.multipletests(p_values.values, method='fdr_bh', alpha=0.001)
    # get the genes responsible for the observed changes
    references.ensembl_to_gene_symbol(the_insc.loc[tmp[0]].index, tax_id=10090)

    # compare within mice

    data = obj.data.loc[obj.data.index.str.contains('ENS')]
    meta = obj.meta

    # cpm = data.divide(meta.loc[:, 'read_count'].values, axis=1) * 1e6
    cpm = data.divide(data.sum(axis=0), axis=1) * 1e6
    keep = (cpm > .5).sum(axis=1) > 5

    the_dat_cv = np.log2(data.loc[keep] + 1)

    groups_by_mouse = [
        ['eNSC3med', 'eNSC3mouse', 'mDura3N1mouse', 'mDura3N1human'],
        ['eNSC5med', 'eNSC5mouse', 'mDura5N24Amouse', 'mDura5N24Ahuman'],
        ['eNSC6med', 'eNSC6mouse', 'mDura6N6mouse', 'mDura6N6human'],
    ]

    # create a dataframe for bar plotting
    within_cols = ["mouse", "correlation", "comparison"]
    the_corr = the_dat_cv.corr()
    within_dat = []
    for i, g in zip([3, 5, 6], groups_by_mouse):
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

    # different representation of same data
    lbl = ['eNSCmouse', 'iNSCmouse', 'iNSChuman']
    aa = pd.DataFrame(index=samples[3:], columns=samples[:3])
    for t in samples[:3]:
        aa.loc[:, t] = the_dat_cv.loc[:, samples[3:]].corrwith(the_dat.loc[:, t])

    val_mean = []
    val_std = []
    val_dat = []
    for t in groups_by_protocol[1:]:
        tmp = aa.loc[t[1], groups_by_protocol[0][1]].values.diagonal()
        val_dat.append(tmp)
        val_mean.append(tmp.mean())
        val_std.append(tmp.std())

    fig = plt.figure(figsize=(3.4, 5.5))
    ax = fig.add_subplot(111)
    ax.bar(range(len(val_mean)), val_mean, 0.7)
    for i, m in enumerate(val_mean):
        ax.plot([i, i], [m - val_std[i], m + val_std[i]], 'k-')
    ax.set_ylabel("Correlation with endogenous NSC")
    ax.set_xticks(range(len(val_mean)))
    ax.set_xticklabels(lbl, rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'corr_with_ensc.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'corr_with_ensc.pdf'))

    # all possible cross-mouse comparisons

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

    ax = pwise_corr_boxplot(the_dat_cv)
    ax.figure.savefig(os.path.join(outdir, "pairwise_correlation.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "pairwise_correlation.pdf"))