import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from load_data import microarray_data, rnaseq_data
from microarray import process
from scipy.cluster import hierarchy
from scripts.comparison_rnaseq_microarray import consts
from scripts.mb_subgroup_classifier.load import load_xz_rnaseq
import collections
from scripts.output import unique_output_dir
import os


def plot_clustermap(
        dat,
        show_gene_labels=False,
        **kwargs
):

    cg = sns.clustermap(
        dat,
        **kwargs
    )
    # check whether x ticks were requested - if so, rotate and raise the bottom of the axes
    if kwargs.get('xticklabels', False):
        plt.setp(cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90)
        bottom = 0.1

    else:
        bottom = 0.02
    # remove useless row dendrogram
    cg.ax_row_dendrogram.set_visible(False)
    # and remove the space created for it
    wr = cg.gs.get_width_ratios()
    wr[0] = 0.035
    wr[1] = 0.02
    cg.gs.set_width_ratios(wr)
    # reduce whitespace
    cg.gs.update(bottom=bottom, top=0.98, left=0.02)

    cg.ax_heatmap.yaxis.label.set_visible(False)
    cg.ax_heatmap.xaxis.label.set_visible(False)
    if show_gene_labels:
        plt.setp(
            cg.ax_heatmap.yaxis.get_ticklabels(),
            rotation=0,
            fontsize=14
        )
    else:
        cg.ax_heatmap.yaxis.set_ticklabels([])

    return cg


def show_dendrogram_cut_by_distance(clustermap, dist, axis=1):
    """
    Add a dashed line to the specified dendrogram showing the cutoff point.
    The line is specified by distance
    :param clustermap:
    :param dist:
    :param axis: 0 (row) or 1 (column)
    :return: None
    """
    if axis not in (0, 1):
        raise ValueError("Valid choices for `axis` are (0, 1).")

    if axis == 0:
        dend = clustermap.dendrogram_row
        ax = clustermap.ax_row_dendrogram
    else:
        dend = clustermap.dendrogram_col
        ax = clustermap.ax_col_dendrogram

    xmin = np.array(dend.independent_coord).min()
    xmax = np.array(dend.independent_coord).max()
    rng = xmax - xmin

    plot_kwargs = dict(
        linestyle='--',
        color='k',
        alpha=0.3
    )

    if axis == 0:
        ax.plot(
            [dist, dist],
            [xmin - 0.1 * rng, xmax + 0.1 * rng],
            **plot_kwargs
        )
    else:
        ax.plot(
            [xmin - 0.1 * rng, xmax + 0.1 * rng],
            [dist, dist],
            **plot_kwargs
        )


def show_dendrogram_cut_by_nclust(clustermap, nclust, axis=1, sample_subset_idx=None):
    """
    Add a dashed line to the specified dendrogram showing the cutoff point.
    The line is specified by the desired final number of clusters
    :param clustermap:
    :param nclust: The final number of clusters required. The distance is computed using this.
    :param axis: 0 (row) or 1 (column)
    :param sample_subset_idx: If supplied, this is a list of indices pointing to the samples that we should use. In this
    case, the distance is computed at the point where _only those leaf nodes_ have formed nclust clusters.
    :return: None
    """
    if axis == 0:
        dend = clustermap.dendrogram_row
    else:
        dend = clustermap.dendrogram_col
    z = dend.linkage

    if sample_subset_idx is None:
        # take the mean distance between the nclust-1 and nclust levels
        dist = z[-nclust:(2 - nclust), 2].mean()
    else:
        # call the distance based on a fixed number of clusters DEFINED FOR A SUBSET OF NODES
        node_ids = set(sample_subset_idx)
        n = len(z) + 1
        cutoff_idx0 = None
        cutoff_idx1 = None
        # loop through linkage rows, keeping track of the specified leaf nodes
        for i, (l0, l1) in enumerate(z[:, :2]):
            l0 = int(l0)
            l1 = int(l1)
            if l0 in node_ids:
                node_ids.remove(l0)
                node_ids.add(n + i)
            if l1 in node_ids:
                node_ids.remove(l1)
                # fine to do this twice since we are working with a set
                node_ids.add(n + i)
            if len(node_ids) == nclust:
                cutoff_idx0 = i
            if len(node_ids) == (nclust - 1):
                cutoff_idx1 = i
                break
        if cutoff_idx0 is None or cutoff_idx1 is None:
            raise ValueError("Failed to compute the requested cluster distance")
        dist = z[[cutoff_idx0, cutoff_idx1], 2].mean()

    show_dendrogram_cut_by_distance(clustermap, dist, axis=axis)

    return dist


if __name__ == '__main__':
    SAVE_FIG = True
    geneset = consts.NORTHCOTT_GENES
    show_gene_labels = False
    n_genes = 1500
    AGGR_METHOD = 'max_std'
    # Z_SCORE = 0  # normalise gene expression levels by rows
    Z_SCORE = 1  # normalise gene expr levels by patient
    # INCLUDE_ALL_1299 = True  # Include Scramble, shBMI1, shCHD7, shBMI1CHD7
    INCLUDE_ALL_1299 = False  # Include Scramble only

    subgroup_colours = {
        'WNT': '#2D438E',
        'SHH': '#E5161C',
        'Group C': '#F2EA00',
        'Group D': '#2A8C43'
    }

    all_nstring = []
    [all_nstring.extend(t) for _, t in consts.NANOSTRING_GENES]
    all_ncott = []
    [all_ncott.extend(t) for _, t in consts.NORTHCOTT_GENES]

    if SAVE_FIG:
        outdir = unique_output_dir('hierarchical_clustering', reuse_empty=True)

    # Load data
    # Where redundant probes are present, use the one with the highest stdev

    # Thompson unannotated dataset
    # data, _ = microarray_data.load_annotated_thompson2006(aggr_field='SYMBOL', aggr_method='max_std')

    # Robinson annotated
    STUDY = 'Robinson'
    data, meta = microarray_data.load_annotated_microarray_gse37418(aggr_field='SYMBOL', aggr_method='max_std')
    meta = meta.loc[meta.subgroup.isin(['WNT', 'SHH', 'G3', 'G4'])]
    meta.subgroup = meta.subgroup.str.replace('G3', 'Group C')
    meta.subgroup = meta.subgroup.str.replace('G4', 'Group D')
    meta.loc[:, 'study'] = STUDY
    data = data.loc[:, meta.index]
    if 'EYS' in data.index:
        idx = data.index.str.replace('EYS', 'EGFL11')
        data.index = idx

    # Kool
    # STUDY = 'Kool'
    # data, meta = microarray_data.load_annotated_microarray_gse10327(aggr_field='SYMBOL', aggr_method='max_std')

    # Northcott
    # STUDY = 'Northcott'
    # data, meta = microarray_data.load_annotated_microarray_gse37382(aggr_field='SYMBOL', aggr_method='max_std')

    # find top genes by MAD - all genes included
    mad = process.median_absolute_deviation(data, axis=1)
    top_genes = mad.sort_values(ascending=False).index[:n_genes]
    print "Selecting top %d genes by MAD from %s study..." % (n_genes, STUDY)
    print "%d / %d genes (nanostring)" % (len(top_genes.intersection(all_nstring)), len(all_nstring))
    print "%d / %d genes (northcott)" % (len(top_genes.intersection(all_ncott)), len(all_ncott))

    # Zhao data
    zhao_sample_names = (
        'Pt1299',
        'Pt1487',
        'Pt1595',
        'ICb1299-III',
        'ICb1299-IV',
        'ICb1487-I',
        'ICb1487-III',
        'ICb1595-I',
        'ICb1595-III',
    )
    # data_zhao, meta_zhao = microarray_data.load_annotated_gse28192(
    #     sample_names=zhao_sample_names,
    #     log2=True,
    #     aggr_field='SYMBOL',
    #     aggr_method=AGGR_METHOD
    # )
    data_zhao, meta_zhao = microarray_data.load_annotated_gse28192(
        sample_names=zhao_sample_names,
        log2=False
    )
    data_zhao = process.variance_stabilizing_transform(data_zhao)
    data_zhao = data_zhao.loc[:, zhao_sample_names]
    meta_zhao = meta_zhao.loc[zhao_sample_names, :]
    data_zhao = microarray_data.annotate_and_aggregate_gse28192(data_zhao, aggr_field='SYMBOL', aggr_method=AGGR_METHOD)

    # prepare Zhao meta
    meta_zhao.loc[:, 'northcott classification'] = meta_zhao.loc[:, 'northcott classification'].str.replace('C', 'Group C')
    meta_zhao.loc[:, 'northcott classification'] = meta_zhao.loc[:, 'northcott classification'].str.replace('D', 'Group D')
    meta_zhao.loc[:, 'subgroup'] = meta_zhao.loc[:, 'northcott classification']
    meta_zhao.loc[:, 'study'] = 'Zhao'

    # RNA-Seq count data (FPKM): 1595 late passage
    data_sb, _ = rnaseq_data.mb_zhao_cultures(units='fpkm', annotate_by='Approved Symbol')
    data_sb = data_sb.loc[:, ['1595']]
    data_sb.columns = ['Zhang 1595']
    data_sb = np.log2(data_sb + 1.)
    meta_sb = pd.DataFrame(data=[['Zhang', 'Group C']], index=['Zhang 1595'], columns=['study', 'subgroup'])

    # RNA-Seq count data (FPKM): 1299 late passage
    data_xz = load_xz_rnaseq(kind='cuff', yugene=False).fillna(0.)
    if not INCLUDE_ALL_1299:
        # restrict to single control sample
        data_xz = data_xz.loc[:, ['Scramble.1', 'Scramble.2']]
        data_xz.columns = ['Zhang 1299', 'Zhang 1299 (rpt)']
    else:
        idx = data_xz.columns.str.replace(r'BMI1_CHD7', 'BC')
        idx = idx.str.replace(r'BMI1', 'B')
        idx = idx.str.replace(r'CHD7', 'C')
        idx = idx.str.replace(r'Scramble', 'Ctrl')
        data_xz.columns = idx
    data_xz = np.log2(data_xz + 1.)
    meta_xz = pd.DataFrame(index=data_xz.columns, columns=['study', 'subgroup'])
    meta_xz.loc[:, 'study'] = 'Zhang'
    meta_xz.loc[:, 'subgroup'] = 'Group D'

    # matching genes
    common_genes = (
        data.index.intersection(data_zhao.index)
            .intersection(data_sb.index)
            .intersection(data_xz.index)
    )

    # find top genes by MAD of classifying dataset
    mad = process.median_absolute_deviation(data.loc[common_genes, :], axis=1)
    top_genes = mad.sort_values(ascending=False).index[:n_genes]

    # plot: key gene expression plus clustering
    alternatives = {
        'EGFL11': 'EYS'
    }

    gene_class = {}

    g_nano = []
    gcount_nano = collections.Counter()

    for cls, gs in consts.NANOSTRING_GENES:
        for gg in gs:
            # if gg in data.index:
            if gg in common_genes:
                g_nano.append(gg)
                gcount_nano[cls] += 1
            elif gg in alternatives and alternatives[gg] in data.index:
                g_nano.append(alternatives[gg])
                gcount_nano[cls] += 1

    g_ncot = []

    for cls, gs in consts.NORTHCOTT_GENES:
        for gg in gs:
            gene_class[gg] = cls
            # if gg in data.index:
            if gg in common_genes:
                g_ncot.append(gg)
            elif gg in alternatives and alternatives[gg] in data.index:
                g_ncot.append(alternatives[gg])


    # Reference cohort expression data and linkage
    expr_topmad = data.loc[top_genes]
    z_topmad = hierarchy.linkage(expr_topmad.transpose(), method='average', metric='correlation')
    expr_nano = data.loc[g_nano]
    z_nano = hierarchy.linkage(expr_nano.transpose(), method='average', metric='correlation')
    expr_ncot = data.loc[g_ncot]
    z_ncot = hierarchy.linkage(expr_ncot.transpose(), method='average', metric='correlation')

    # integrate all data sources

    all_meta = [
        meta,
        meta_zhao,
        meta_xz,
        meta_sb
    ]
    meta_int = pd.concat(all_meta, axis=0).loc[:, ['study', 'subgroup']]

    all_data = [
        data,
        data_zhao,
        data_xz,
        data_sb
    ]
    data_int = pd.concat(all_data, axis=1).dropna()

    expr_int_topmad = data_int.loc[top_genes].dropna()
    z_int_topmad = hierarchy.linkage(expr_int_topmad.transpose(), method='average', metric='correlation')
    expr_int_nano = data_int.loc[g_nano].dropna()
    z_int_nano = hierarchy.linkage(expr_int_nano.transpose(), method='average', metric='correlation')
    expr_int_ncot = data_int.loc[g_ncot].dropna()
    z_int_ncot = hierarchy.linkage(expr_int_ncot.transpose(), method='average', metric='correlation')

    # maintain a list of original sample IDs for distance computation in clustering
    sample_idx = np.where(meta_int.index.str.startswith('GSM'))[0]

    row_colors_nano = pd.Series(dict([(k, gene_class[k]) for k in g_nano]), name='').apply(lambda x: subgroup_colours[x])
    row_colors_ncot = pd.Series(gene_class, name='').apply(lambda x: subgroup_colours[x])

    # now we assume that 4 clusters have been found and distinguish by greyscale shade
    cluster_colours = {
        1: '#000000',
        2: '#a5a5a5',
        3: '#6d6d6d',
        4: '#dbdbdb',
    }
    col_colors = pd.DataFrame(index=data.columns, columns=['inferred', 'subgroup'])
    col_colors.loc[:, 'subgroup'] = meta.subgroup.apply(subgroup_colours.get)
    col_colors.loc[:, 'inferred'] = pd.Series(
        hierarchy.fcluster(z_nano, 4, criterion='maxclust'),
        index=data.columns
    ).apply(cluster_colours.get)

    # by top MAD genes
    col_colors.loc[:, 'inferred'] = pd.Series(
        hierarchy.fcluster(z_topmad, 4, criterion='maxclust'),
        index=data.columns
    ).apply(cluster_colours.get)

    cg = plot_clustermap(
        expr_topmad,
        cmap='RdBu_r',
        col_colors=col_colors,
        col_linkage=z_topmad,
        z_score=0,
        xticklabels=False,
    )
    show_dendrogram_cut_by_nclust(cg, 4)
    ttl = "%s_top%d_heatmap" % (STUDY.lower(), n_genes)
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        # cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))  # slow due to large number of cells

    # by ncott
    col_colors.loc[:, 'inferred'] = pd.Series(
        hierarchy.fcluster(z_ncot, 4, criterion='maxclust'),
        index=data.columns
    ).apply(cluster_colours.get)

    cg = plot_clustermap(
        expr_ncot,
        cmap='RdBu_r',
        show_gene_labels=True,
        row_colors=row_colors_ncot,
        col_colors=col_colors,
        row_cluster=None,
        col_linkage=z_ncot,
        z_score=Z_SCORE,
        xticklabels=False,
    )
    plt.setp(
        cg.ax_heatmap.yaxis.get_ticklabels(),
        rotation=0,
        fontsize=8
    )
    show_dendrogram_cut_by_nclust(cg, 4)
    ttl = "%s_ncott_heatmap" % STUDY.lower()
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    # by nano
    col_colors.loc[:, 'inferred'] = pd.Series(
        hierarchy.fcluster(z_nano, 4, criterion='maxclust'),
        index=data.columns
    ).apply(cluster_colours.get)

    cg = plot_clustermap(
        expr_nano,
        show_gene_labels=True,
        cmap='RdBu_r',
        row_colors=row_colors_nano,
        col_colors=col_colors,
        row_cluster=None,
        col_linkage=z_nano,
        z_score=Z_SCORE,
        xticklabels=False,
    )
    cut_dist_nano = show_dendrogram_cut_by_nclust(cg, 4)
    ttl = "%s_nano_heatmap" % STUDY.lower()
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    study_colours = {
        STUDY: '#cccccc',
        'Zhao': '#000000',
        'Zhang': '#666666',
    }

    col_colors_int = pd.DataFrame(index=data_int.columns, columns=['study', 'subgroup'])
    col_colors_int.loc[:, 'subgroup'] = meta_int.subgroup.apply(subgroup_colours.get)
    col_colors_int.loc[:, 'study'] = meta_int.study.apply(study_colours.get)

    # by top MAD genes
    cg = plot_clustermap(
        expr_int_topmad,
        cmap='RdBu_r',
        col_colors=col_colors_int,
        col_linkage=z_int_topmad,
        z_score=Z_SCORE,
        xticklabels=True,
    )
    # add dashed line to show cutoff
    show_dendrogram_cut_by_nclust(cg, 4, sample_subset_idx=sample_idx)
    # show_dendrogram_cut_by_distance(cg, 0.55, axis=1)
    ttl = "%s_zhao_zhang_heatmap_cluster_by_top%d" % (STUDY.lower(), n_genes)
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        # cg.savefig(os.path.join(outdir, "%s.pdf" % ttl)) # slow due to large number of cells

    # by ncott
    cg = plot_clustermap(
        expr_int_ncot,
        cmap='RdBu_r',
        show_gene_labels=True,
        row_cluster=None,
        row_colors=row_colors_ncot,
        col_colors=col_colors_int,
        col_linkage=z_int_ncot,
        z_score=Z_SCORE,
        xticklabels=True,
    )
    plt.setp(
        cg.ax_heatmap.yaxis.get_ticklabels(),
        rotation=0,
        fontsize=8
    )
    # add dashed line to show cutoff
    # show_dendrogram_cut_by_nclust(cg, 4, sample_subset_idx=sample_idx)
    show_dendrogram_cut_by_nclust(cg, 5)
    # show_dendrogram_cut_by_distance(cg, 0.5, axis=1)
    ttl = "%s_zhao_zhang_heatmap_cluster_by_ncot" % STUDY.lower()
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    # by nano
    cg = plot_clustermap(
        expr_int_nano,
        show_gene_labels=True,
        cmap='RdBu_r',
        row_colors=row_colors_nano,
        row_cluster=None,
        col_colors=col_colors_int,
        col_linkage=z_int_nano,
        z_score=Z_SCORE,
        xticklabels=True,
    )
    # add dashed line to show cutoff
    show_dendrogram_cut_by_nclust(cg, 4, sample_subset_idx=sample_idx)
    # show_dendrogram_cut_by_distance(cg, 0.5, axis=1)
    ttl = "%s_zhao_zhang_heatmap_cluster_by_nano" % STUDY.lower()
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))
