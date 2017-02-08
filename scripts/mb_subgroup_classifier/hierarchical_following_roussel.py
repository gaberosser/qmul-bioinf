import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from load_data import microarray_data, rnaseq_data
from microarray import process
from scipy.cluster import hierarchy
from scripts.comparison_rnaseq_microarray import consts
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
        cutoff_idx = None
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
            if len(node_ids) == (nclust - 1):
                cutoff_idx = i
                break
        if cutoff_idx is None:
            raise ValueError("Failed to compute the requested cluster distance")
        dist = z[cutoff_idx - 1:cutoff_idx + 1, 2].mean()

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

    # pare down zhao meta
    meta_zhao.loc[:, 'northcott classification'] = meta_zhao.loc[:, 'northcott classification'].str.replace('C', 'Group C')
    meta_zhao.loc[:, 'northcott classification'] = meta_zhao.loc[:, 'northcott classification'].str.replace('D', 'Group D')
    mz = pd.DataFrame(index=meta_zhao.index, columns=['subgroup', 'study'])
    mz.loc[:, 'study'] = 'Zhao'
    mz.loc[:, 'subgroup'] = meta_zhao.loc[:, 'northcott classification']

    # matching genes
    common_genes = data.index.intersection(data_zhao.index)

    # find top genes by MAD of classifying dataset
    mad = process.median_absolute_deviation(data.loc[common_genes, :], axis=1)
    top_genes = mad.sort_values(ascending=False).index[:n_genes]

    X = data.loc[top_genes]

    # plot: key gene expression plus clustering
    alternatives = {
        'EGFL11': 'EYS'
    }

    g_nano = []
    gcount_nano = collections.Counter()

    for cls, gs in consts.NANOSTRING_GENES:
        for gg in gs:
            if gg in data.index:
                g_nano.append(gg)
                gcount_nano[cls] += 1
            elif gg in alternatives and alternatives[gg] in data.index:
                g_nano.append(alternatives[gg])
                gcount_nano[cls] += 1

    g_ncot = []
    gcount_ncot = collections.Counter()

    for cls, gs in consts.NORTHCOTT_GENES:
        for gg in gs:
            if gg in data.index:
                g_ncot.append(gg)
                gcount_ncot[cls] += 1
            elif gg in alternatives and alternatives[gg] in data.index:
                g_ncot.append(alternatives[gg])
                gcount_ncot[cls] += 1

    expr_nano = data.loc[g_nano]
    expr_ncot = data.loc[g_ncot]

    # build row colours Series
    colour_cycle = {
        'WNT': '#2D438E',
        'SHH': '#E5161C',
        'Group C': '#F2EA00',
        'Group D': '#2A8C43'
    }

    rc_nano = []
    rc_ncot = []

    # iterate over the genesets again to ensure classes are in the same order
    for cls, _ in consts.NANOSTRING_GENES:
        rc_nano.extend([colour_cycle[cls]] * gcount_nano[cls])
    row_colors_nano = pd.Series(rc_nano, index=g_nano, name='')

    for cls, _ in consts.NORTHCOTT_GENES:
        rc_ncot.extend([colour_cycle[cls]] * gcount_ncot[cls])
    row_colors_ncot = pd.Series(rc_ncot, index=g_ncot, name='')

    # now we assume that 4 clusters have been found and use those for colours
    col_colors = pd.DataFrame(index=X.columns, columns=['inferred', 'labelled'])
    for cls in colour_cycle:
        col_colors.loc[meta.subgroup == cls, 'labelled'] = colour_cycle[cls]

    # unsupervised hierarchical clustering

    # Nanostring heatmap
    z_nano = hierarchy.linkage(expr_nano.transpose(), method='average', metric='correlation')
    lbl_nano = hierarchy.fcluster(z_nano, 4, criterion='maxclust')
    col_colors.loc[lbl_nano == 1, 'inferred'] = '#000000'
    col_colors.loc[lbl_nano == 2, 'inferred'] = '#a5a5a5'
    col_colors.loc[lbl_nano == 3, 'inferred'] = '#6d6d6d'
    col_colors.loc[lbl_nano == 4, 'inferred'] = '#dbdbdb'

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

    # Northcott heatmap
    z_ncot = hierarchy.linkage(expr_ncot.transpose(), method='average', metric='correlation')
    lbl_ncot = hierarchy.fcluster(z_ncot, 4, criterion='maxclust')
    col_colors.loc[lbl_ncot == 1, 'inferred'] = '#000000'
    col_colors.loc[lbl_ncot == 2, 'inferred'] = '#a5a5a5'
    col_colors.loc[lbl_ncot == 3, 'inferred'] = '#6d6d6d'
    col_colors.loc[lbl_ncot == 4, 'inferred'] = '#dbdbdb'

    cg = plot_clustermap(
        expr_ncot,
        cmap='RdBu_r',
        row_colors=row_colors_ncot,
        col_colors=col_colors,
        row_cluster=None,
        col_linkage=z_ncot,
        z_score=Z_SCORE,
        xticklabels=False,
    )
    cut_dist_ncot = show_dendrogram_cut_by_nclust(cg, 4)
    ttl = "%s_ncott_heatmap" % STUDY.lower()
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    # Global (n_genes) heatmap
    z_all = hierarchy.linkage(X.transpose(), method='average', metric='correlation')
    lbl_all = hierarchy.fcluster(z_all, 4, criterion='maxclust')
    col_colors.loc[lbl_all == 1, 'inferred'] = '#000000'
    col_colors.loc[lbl_all == 2, 'inferred'] = '#a5a5a5'
    col_colors.loc[lbl_all == 3, 'inferred'] = '#6d6d6d'
    col_colors.loc[lbl_all == 4, 'inferred'] = '#dbdbdb'
    cg = plot_clustermap(
        X,
        cmap='RdBu_r',
        col_colors=col_colors,
        col_linkage=z_all,
        z_score=0,
        xticklabels=False,
    )
    cut_dist_all = show_dendrogram_cut_by_nclust(cg, 4)
    ttl = "%s_top%d_heatmap" % (STUDY.lower(), n_genes)
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        # cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))  # slow due to large number of cells

    # integrate the Zhao data
    # only keep common genes BUT use the same MAD genes as before

    meta_int = pd.concat((meta, mz), axis=0)
    data_int = pd.concat(
        (data.loc[common_genes, :], data_zhao.loc[common_genes, :]),
        axis=1
    )
    X_int = data_int.loc[top_genes]
    expr_int_nano = data_int.loc[g_nano].dropna()
    expr_int_ncot = data_int.loc[g_ncot].dropna()

    z_int = hierarchy.linkage(X_int.transpose(), method='average', metric='correlation')

    # maintain a list of original saple IDs for distance computation in clustering
    sample_idx = np.where(meta_int.index.str.startswith('GSM'))[0]

    col_colors_int = pd.DataFrame(index=data_int.columns, columns=['study', 'labelled'])

    col_colors_int.loc[meta_int.study == STUDY] = '#dbdbdb'
    col_colors_int.loc[meta_int.study == 'Zhao'] = '#000000'
    for cls in colour_cycle:
        col_colors_int.loc[meta_int.subgroup == cls, 'labelled'] = colour_cycle[cls]

    cg = plot_clustermap(
        X_int,
        cmap='RdBu_r',
        col_colors=col_colors_int,
        col_linkage=z_int,
        z_score=Z_SCORE,
        xticklabels=True,
    )
    # add dashed line to show cutoff
    show_dendrogram_cut_by_nclust(cg, 4, sample_subset_idx=sample_idx)
    # show_dendrogram_cut_by_distance(cg, cut_dist_all, axis=1)
    # show_dendrogram_cut_by_distance(cg, 0.55, axis=1)
    ttl = "%s_zhao_top%d_heatmap" % (STUDY.lower(), n_genes)
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        # cg.savefig(os.path.join(outdir, "%s.pdf" % ttl)) # slow due to large number of cells

    # try clustering with only Ncott or Nano
    z_int_nano = hierarchy.linkage(expr_int_nano.transpose(), method='average', metric='correlation')

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

    # try a new approach - call the distance based on a fixed number of clusters DEFINED FOR A SUBSET OF NODES
    show_dendrogram_cut_by_nclust(cg, 4, sample_subset_idx=sample_idx)
    # show_dendrogram_cut_by_distance(cg, cut_dist_nano, axis=1)
    # show_dendrogram_cut_by_distance(cg, 0.5, axis=1)

    ttl = "%s_zhao_ncott_heatmap_cluster_by_nano" % STUDY.lower()
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    z_int_ncot = hierarchy.linkage(expr_int_ncot.transpose(), method='average', metric='correlation')

    cg = plot_clustermap(
        expr_int_ncot,
        cmap='RdBu_r',
        row_cluster=None,
        col_colors=col_colors_int,
        col_linkage=z_int_ncot,
        z_score=Z_SCORE,
        xticklabels=True,
    )
    # add dashed line to show cutoff
    show_dendrogram_cut_by_nclust(cg, 4, sample_subset_idx=sample_idx)
    # show_dendrogram_cut_by_distance(cg, cut_dist_ncot, axis=1)
    # show_dendrogram_cut_by_distance(cg, 0.5, axis=1)
    ttl = "%s_zhao_ncott_heatmap_cluster_by_ncot" % STUDY.lower()
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))


    # load RNA-Seq count data
    from scripts.mb_subgroup_classifier.load import load_xz_rnaseq
    xz_expr = load_xz_rnaseq(kind='cuff', yugene=False).fillna(0.)

    # variance stabilising transform
    xz_vst = process.variance_stabilizing_transform(xz_expr)

    # integrate (dropping genes that have missing values)
    data_int = pd.concat(
        (data, data_zhao, xz_vst),
        axis=1
    ).dropna(axis=0)
    meta_xz = pd.DataFrame(index=xz_expr.columns, columns=['study', 'subgroup'])
    meta_xz.loc[:, 'study'] = 'Zhang'
    meta_xz.loc[:, 'subgroup'] = 'Group D'

    meta_int = pd.concat(
        (meta, mz, meta_xz),
        axis=0
    )

    col_colors_int = pd.DataFrame(index=data_int.columns, columns=['study', 'labelled'])

    col_colors_int.loc[meta_int.study == STUDY] = '#cccccc'
    col_colors_int.loc[meta_int.study == 'Zhao'] = '#000000'
    col_colors_int.loc[meta_int.study == 'Zhang'] = '#666666'
    for cls in colour_cycle:
        col_colors_int.loc[meta_int.subgroup == cls, 'labelled'] = colour_cycle[cls]

    X_int = data_int.loc[top_genes].dropna()
    expr_int_nano = data_int.loc[g_nano].dropna()
    expr_int_ncot = data_int.loc[g_ncot].dropna()

    z_int = hierarchy.linkage(X_int.transpose(), method='average', metric='correlation')

    cg = plot_clustermap(
        X_int,
        cmap='RdBu_r',
        col_colors=col_colors_int,
        col_linkage=z_int,
        z_score=Z_SCORE,
        xticklabels=True,
    )
    # add dashed line to show cutoff
    show_dendrogram_cut_by_distance(cg, cut_dist_all, axis=1)  # 0.55
    ttl = "%s_zhao_xz_top%d_heatmap" % (STUDY.lower(), n_genes)
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        # cg.savefig(os.path.join(outdir, "%s.pdf" % ttl)) # slow due to large number of cells

    # try clustering with only Ncott or Nano
    z_int_nano = hierarchy.linkage(expr_int_nano.transpose(), method='average', metric='correlation')

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
    show_dendrogram_cut_by_distance(cg, cut_dist_nano, axis=1)  # 0.5
    ttl = "%s_zhao_xz_ncott_heatmap_cluster_by_nano" % STUDY.lower()
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    z_int_ncot = hierarchy.linkage(expr_int_ncot.transpose(), method='average', metric='correlation')

    cg = plot_clustermap(
        expr_int_ncot,
        cmap='RdBu_r',
        row_cluster=None,
        col_colors=col_colors_int,
        col_linkage=z_int_ncot,
        z_score=Z_SCORE,
        xticklabels=True,
    )
    # add dashed line to show cutoff
    show_dendrogram_cut_by_distance(cg, cut_dist_ncot, axis=1)  # 0.5
    ttl = "%s_zhao_xz_ncott_heatmap_cluster_by_ncot" % STUDY.lower()
    if SAVE_FIG:
        cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
        cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))