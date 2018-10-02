import pandas as pd
import os
import re
from settings import HGIC_LOCAL_DIR
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import clustering
from utils import output
from scipy import stats


if __name__ == '__main__':
    alpha = 0.05
    renorm = False
    # minimum number of samples with that cell type detectable
    min_n_samples = 2

    fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/characterisation/ffpe_cell_type_determination/xcell/xcell_results_salmon_tpm.xlsx'
    )
    outdir = output.unique_output_dir()

    res = pd.read_excel(fn, sheet_name=None)
    df = res['Proportions'].copy()
    df_pval = res['P values'].copy()

    # extract then drop the PID row
    pid = df.loc['Patient ID'].copy()
    df.drop('Patient ID', axis=0, inplace=True)
    df_pval.drop('Patient ID', axis=0, inplace=True)


    # zero out any non-significant results
    df[df_pval > alpha] = 0.

    # improve naming clarity
    df.columns = df.columns.str.replace(r'(?P<n>[0-9])DEF', '\g<n>_DEF')
    df.columns = [re.sub(r'NH1[56]_[0-9]*_', '%s_' % u, t) for t, u in zip(df.columns, pid.values)]
    df_pval.columns = df.columns

    # separate the 'scores' (last 3 rows)
    df_scores = df.iloc[-3:]
    df = df.iloc[:-3]

    # number of samples with detectable levels of cell types
    fig = plt.figure(figsize=(6, 9))
    ax = fig.add_subplot(111)
    (df != 0).sum(axis=1).sort_values(ascending=False).plot.barh(ax=ax, color='k')
    ax.set_xlabel("Number of samples with detectable levels (%d total)" % df.shape[1])
    ax.set_xlim([0, df.shape[1]])
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "number_detectable_by_cell_type.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "number_detectable_by_cell_type.tiff"), dpi=200)
    fig.savefig(os.path.join(outdir, "number_detectable_by_cell_type.pdf"), dpi=200)

    # filter based on this
    keep_ix = (df != 0).sum(axis=1) > min_n_samples
    print "The following %d cell types do not meet the requirement for detection in %d or more samples: \n%s" % (
        (~keep_ix).sum(),
        min_n_samples,
        ', '.join(keep_ix.index[keep_ix])
    )

    df = df.loc[keep_ix]

    print "%d cell types remain." % df.shape[0]

    # confused here: the column sums are not even close to 1.0
    # do we need to renormalise??
    if renorm:
        df = df.divide(df.sum(axis=0), axis=1)

    # heatmap: proportions for each patient
    # standardise across columns, because each cell type has different mean proportion
    cg = clustering.plot_clustermap(
        df.astype(float).transpose(),
        metric='euclidean',
        show_gene_labels=True,
        show_gene_clustering=True,
        cmap='RdBu_r',
        z_score=1
    )
    cg.gs.update(bottom=0.22, right=0.8)

    # box and whisker
    palette = sns.color_palette("deep", 16)
    fig = plt.figure(figsize=(17, 5.5))
    ax = fig.add_subplot(111)
    sns.boxplot(
        data=df.transpose(),
        palette=palette,
        fliersize=0,
        ax=ax
    )
    sns.stripplot(
        data=df.transpose(),
        palette=palette,
        jitter=0.3,
        ax=ax,
        edgecolor='k',
        linewidth=0.5,
    )
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    fig.subplots_adjust(bottom=0.5, left=0.05, right=0.99, top=0.98)
    fig.savefig(os.path.join(outdir, "box_and_swarmplot_all_types.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "box_and_swarmplot_all_types.tiff"), dpi=200)

    # stdev across samples
    stdev = df.std(axis=1).sort_values(ascending=False)

    fig = plt.figure(figsize=(6, 9))
    ax = fig.add_subplot(111)
    stdev.plot.barh(color='k', ax=ax)
    fig.subplots_adjust(bottom=0.07, left=0.4, right=0.98, top=0.98)
    ax.set_xlabel('Stdev. across samples')
    fig.savefig(os.path.join(outdir, "stdev_across_samples.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "stdev_across_samples.tiff"), dpi=200)

    # CV across samples (may be more relevant?)
    cv = (df.std(axis=1) / df.mean(axis=1)).sort_values(ascending=False)
    fig = plt.figure(figsize=(6, 9))
    ax = fig.add_subplot(111)
    cv.plot.barh(color='k', ax=ax)
    fig.subplots_adjust(bottom=0.07, left=0.4, right=0.98, top=0.98)
    ax.set_xlabel('CV across samples')
    fig.savefig(os.path.join(outdir, "cv_across_samples.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "cv_across_samples.tiff"), dpi=200)

    # correlation between inflammatory response pathway (?) and different cell types
    corr_metric = 'spearman'

    fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/core_pipeline/rnaseq/merged_s1_s2/ipa/pathways',
        'ipa_results_s2_de_relevant.xlsx'
    )
    ipa_res = pd.read_excel(fn)
    p = ipa_res.loc[:, ipa_res.columns.str.contains('_syngeneic_-logp')]
    p.columns = p.columns.str.replace('_syngeneic_-logp', '')
    z = ipa_res.loc[:, ipa_res.columns.str.contains('_syngeneic_z')]
    z.columns = z.columns.str.replace('_syngeneic_z', '')

    # pathways of interest
    pathways = p.index
    # pathways = [
    #     'Antigen Presentation Pathway',
    #     'LXR/RXR Activation',
    #     'EIF2 Signaling',
    # ]

    # rename the columns, to enable cross-indexing
    tt = df.columns.str.replace('GBM', '').str.replace(r'_.*', '')
    df.columns = tt

    co = pd.DataFrame(index=df.index, columns=pathways, dtype=float)
    co_p = pd.DataFrame(index=df.index, columns=pathways, dtype=float)

    for pw in pathways:
        this_p = p.loc[pw].sort_index()
        this_z = z.loc[pw].sort_index()
        this_df = df.loc[:, this_p.index].sort_index(axis=1)

        for ct in df.index:
            if corr_metric == 'spearman':
                co.loc[ct, pw], co_p.loc[ct, pw] = stats.spearmanr(this_p, this_df.loc[ct])
            elif corr_metric == 'pearson':
                co.loc[ct, pw], co_p.loc[ct, pw] = stats.pearsonr(this_p, this_df.loc[ct])

    cg = clustering.plot_clustermap(
        co.transpose(),
        metric='euclidean',
        show_gene_labels=True,
    )
    cg.gs.update(bottom=0.21, right=0.68)
    plt.setp(cg.ax_heatmap.yaxis.get_ticklabels(), fontsize=6)
    sns.heatmap(data=co, cmap='RdBu_r', ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)

    cg.savefig(os.path.join(outdir, "cell_proportion_pathway_pval_%s_clustering.png" % corr_metric), dpi=200)
    cg.savefig(os.path.join(outdir, "cell_proportion_pathway_pval_%s_clustering.tiff" % corr_metric), dpi=200)
    cg.savefig(os.path.join(outdir, "cell_proportion_pathway_pval_%s_clustering.pdf" % corr_metric), dpi=200)