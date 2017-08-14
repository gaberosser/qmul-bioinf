from load_data import rnaseq_data
from microarray import process
from stats import transformations
import pandas as pd
from plotting import clustering, heatmap
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import os
from utils import output


def get_twodata(data, pattern, filter_idx=None):
    if filter_idx is not None:
        twodata = data.loc[filter_idx, obj.data.columns.str.contains(pattern)]
    else:
        twodata = data.loc[:, obj.data.columns.str.contains(pattern)]
    if twodata.shape[1] != 2:
        raise ValueError("Expected two samples for matching entry %s, got %d" % (m, twodata.shape[1]))
    return twodata


if __name__ == "__main__":
    matches = (
        'GBM018',
        'GBM026',
        'GBM044',
        'DURA018',
        'DURA044',
    )

    outdir = output.unique_output_dir("replicate_consistency")

    obj = rnaseq_data.all_hgic_loader(annotate_by='Ensembl Gene ID')
    # reorder for plotting
    obj.data = obj.data.loc[:, obj.data.columns.sort_values()]
    obj.meta = obj.meta.loc[obj.data.columns]

    # filters:
    # 1) real genes only
    # 2) only genes that have a count > 5 in 10 or more samples
    gene_idx = obj.data.index.str.contains('ENSG')
    count_sum = obj.data.loc[gene_idx].sum(axis=0)

    # test the effect of the filter
    co = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 400, 500, 750, 1000]
    nrem = [(gene_idx & ((obj.data > t).sum(axis=1) > 10)).sum() for t in co]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(co, nrem)
    ax.set_xlabel("10 / 21 samples must have count higher than x")
    ax.set_ylabel("Number of genes remaining")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "filtering_effect.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "filtering_effect.pdf"))

    idx = gene_idx & ((obj.data > 10).sum(axis=1) > 10)
    data = obj.data.loc[idx]
    mad = transformations.median_absolute_deviation(data).sort_values(ascending=False)
    logdata = np.log(data + 1)

    # start with a dendrogram
    col_colours = clustering.generate_colour_map_dict(
        obj.meta,
        'sample',
        matches,
        label='Patient',
        non_matching='gray'
    )
    out = clustering.dendrogram_with_colours(logdata, col_colours=col_colours, vertical=False)
    dist = clustering.dendrogram_threshold_by_nclust(out['linkage'], 3)
    out['dendrogram_ax'].axvline(dist, ls='--', c='gray')
    out['fig'].savefig(os.path.join(outdir, "dendrogram_all_genes.png"), dpi=200)
    out['fig'].savefig(os.path.join(outdir, "dendrogram_all_genes.pdf"))

    # repeat but now only use N genes (by MAD)
    # tested and the result is unchanged for most values in the region [500, 5000]
    n_gene = 1500
    out = clustering.dendrogram_with_colours(logdata.loc[mad.index[:n_gene]], col_colours=col_colours, vertical=False)
    dist = clustering.dendrogram_threshold_by_nclust(out['linkage'], 3)
    out['dendrogram_ax'].axvline(dist, ls='--', c='gray')
    out['fig'].savefig(os.path.join(outdir, "dendrogram_top_%d_by_mad.png" % n_gene), dpi=200)
    out['fig'].savefig(os.path.join(outdir, "dendrogram_top_%d_by_mad.pdf" % n_gene))

    # scatterplot showing mean vs MAD before and after
    mad_all = transformations.median_absolute_deviation(np.log10(obj.data.loc[gene_idx] + 1), axis=1)
    mean_all = np.log10(obj.data.loc[gene_idx] + 1).mean(axis=1)
    mad_filt = transformations.median_absolute_deviation(np.log10(obj.data.loc[idx] + 1), axis=1)
    mean_filt = np.log10(obj.data.loc[idx]).mean(axis=1)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(mean_all, mad_all, label='Before filtering')
    ax.scatter(mean_filt, mad_filt, color='r', label='After filtering')
    ax.legend(loc='upper left')
    ax.set_xlabel('Mean log2 gene count')
    ax.set_ylabel('Median absolute deviation log2 gene count')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "mean_vs_mad_before_after_filtering.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "mean_vs_mad_before_after_filtering.pdf"))

    # table of correlation coefficients for matching samples
    pair_corrcoeff = pd.Series(index=matches)
    for m in matches:
        twodata = get_twodata(logdata, m, filter_idx=idx)
        pair_corrcoeff.loc[m] = twodata.corr().iloc[0, 1]

    # look at the identity of the most abundant genes
    idx_sort_by_mean = data.mean(axis=1).sort_values(ascending=False).index
    top_data_sorted = rnaseq_data.annotate(data.loc[idx_sort_by_mean[:30]], annotate_by='Approved Symbol')
    top_data_sorted_n = top_data_sorted.divide(count_sum, axis=1)
    ax, cbar = heatmap.single_heatmap(top_data_sorted_n, orientation='vertical', cmap='Reds', vmin=0, vmax=0.03)
    ax.figure.tight_layout()
    fig.savefig(os.path.join(outdir, "top_genes_by_abundance.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "top_genes_by_abundance.pdf"))

    # look at the identity of the most variable genes
    mad = transformations.median_absolute_deviation(data, axis=1)
    idx_sort_by_std = mad.sort_values(ascending=False).index
    top_data_sorted = rnaseq_data.annotate(data.loc[idx_sort_by_std[:30]], annotate_by='Approved Symbol')
    top_data_sorted_n = top_data_sorted.divide(count_sum, axis=1)
    ax, cbar = heatmap.single_heatmap(top_data_sorted_n, orientation='vertical', cmap='Reds', vmin=0, vmax=0.03)
    ax.figure.tight_layout()
    fig.savefig(os.path.join(outdir, "top_genes_by_deviance.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "top_genes_by_deviance.pdf"))

    # plots of log(absolute difference) vs proportion difference for each pair
    for i, m in enumerate(matches):
        twodata = get_twodata(obj.data, m, filter_idx=gene_idx)
        # require that EITHER both have a minimum count OR one is zero and the other high
        this_idx = (twodata > 10).all(axis=1) | (twodata > 100).any(axis=1)
        twodata = twodata.loc[this_idx]

        abs_diff = np.abs(twodata.iloc[:, 0] - twodata.iloc[:, 1])
        mean_val = twodata.sum(axis=1)
        rel_diff = abs_diff / mean_val
        out = sns.jointplot(np.log2(abs_diff + 1.), rel_diff, alpha=0.3, stat_func=None)
        out.ax_joint.set_ylim([0, 1])
        out.ax_marg_y.set_ylim([0, 1])
        ax = out.ax_joint

        ax.set_xlabel("Log2 absolute difference")
        ax.set_ylabel("Relative difference = absolute difference / mean")
        ax.set_title(m)
        plt.tight_layout()
        ax.figure.savefig(os.path.join(outdir, "abs_rel_difference_%s.png" % m), dpi=200)
        ax.figure.savefig(os.path.join(outdir, "abs_rel_difference_%s.pdf" % m))
        # list the genes for which one is zero
        the_list = rnaseq_data.annotate(twodata.loc[(twodata == 0).any(axis=1)], annotate_by='Approved Symbol')
        if the_list.size > 0:
            print "%s. Genes for which one is zero and one >100: %s" % (m, ', '.join(the_list.index))

    cg = clustering.plot_correlation_clustermap(logdata)
    cg.gs.update(bottom=0.25, right=0.7)
    cg.fig.savefig(os.path.join(outdir, "clustermap_correlation_coeff.png"), dpi=200)
    cg.fig.savefig(os.path.join(outdir, "clustermap_correlation_coeff.pdf"))

