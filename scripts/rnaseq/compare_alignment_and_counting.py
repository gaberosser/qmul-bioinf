from load_data import rnaseq_data
import os
from scipy import stats
from matplotlib import pyplot as plt
from plotting import bar
import seaborn as sns
from utils.output import unique_output_dir


if __name__ == '__main__':
    OUTPUT = unique_output_dir("compare_alignment", reuse_empty=True)
    star = rnaseq_data.gbm_paired_samples_loader(source='star', annotate_by='Ensembl Gene ID')
    fcou = rnaseq_data.gbm_paired_samples_loader(source='featurecounts', annotate_by='Ensembl Gene ID')

    # pct assigned
    total_counts = star.data.sum(axis=0)
    star_assigned = star.data.loc[star.data.index.str.contains('ENSG')].sum() / total_counts
    fcou_assigned = fcou.data.sum() / total_counts

    bar.grouped_bar_chart([star_assigned, fcou_assigned])
    ax = plt.gca()
    ax.set_ylim([0, 1])
    ax.legend(['STAR', 'hisat2 + featureCounts'], loc='upper right')
    ax.set_ylabel('Proportion of reads assigned to a gene')
    ax.set_position([0.125, 0.3, 0.85, 0.65])
    ax.figure.savefig(os.path.join(OUTPUT, 'pct_reads_assigned.png'), dpi=200)
    ax.figure.savefig(os.path.join(OUTPUT, 'pct_reads_assigned.pdf'))

    # correlation between gene counts STAR vs featureCounts
    corr_rsq = []
    corr_p = []
    for t in star.data.columns:
        r, pval = stats.pearsonr(star.data.loc[fcou.data.index, t], fcou.data.loc[:, t])
        corr_rsq.append(r ** 2)
        corr_p.append(pval)

    fig = plt.figure()
    h = plt.scatter(star_assigned, corr_rsq)
    h.axes.set_xlabel('Proportion assigned by STAR')
    h.axes.set_ylabel('R^2 correlation STAR / featureCounts')
    h.axes.set_ylim([0.9, 1.])
    fig.savefig(os.path.join(OUTPUT, 'assigned_vs_correlation.png'), dpi=200)
    fig.savefig(os.path.join(OUTPUT, 'assigned_vs_correlation.pdf'))

    # for comparison: find the correlation between samples for STAR
    rsq = star.data.loc[fcou.data.index].corr() ** 2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.heatmap(rsq, vmin=0, vmax=1., cmap='RdBu_r', ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.tight_layout()
    ax.figure.savefig(os.path.join(OUTPUT, 'star_correlation_between_samples.png'), dpi=200)
    ax.figure.savefig(os.path.join(OUTPUT, 'star_correlation_between_samples.pdf'))