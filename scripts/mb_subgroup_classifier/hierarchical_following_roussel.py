import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from load_data import microarray_data, rnaseq_data
from microarray import process
from scipy.cluster import hierarchy
from scripts.comparison_rnaseq_microarray import consts
import collections


if __name__ == '__main__':
    geneset = consts.NORTHCOTT_GENES
    show_gene_labels = False
    n_genes = 1500

    all_nstring = []
    [all_nstring.extend(t) for _, t in consts.NANOSTRING_GENES]
    all_ncott = []
    [all_ncott.extend(t) for _, t in consts.NORTHCOTT_GENES]

    # Load data
    # Where redundant probes are present, use the one with the highest stdev

    # Thompson unannotated dataset
    # data, _ = microarray_data.load_annotated_thompson2006(aggr_field='SYMBOL', aggr_method='max_std')

    # Robinson annotated
    # data, meta = microarray_data.load_annotated_microarray_gse37418(aggr_field='SYMBOL', aggr_method='max_std')
    # meta = meta.loc[meta.subgroup.isin(['WNT', 'SHH', 'G3', 'G4'])]
    # meta.subgroup = meta.subgroup.str.replace('G3', 'Group C')
    # meta.subgroup = meta.subgroup.str.replace('G4', 'Group D')
    # meta.loc[:, 'study'] = 'Robinson'
    # data = data.loc[:, meta.index]

    # Kool
    # data, meta = microarray_data.load_annotated_microarray_gse10327(aggr_field='SYMBOL', aggr_method='max_std')

    # Northcott
    data, meta = microarray_data.load_annotated_microarray_gse37382(aggr_field='SYMBOL', aggr_method='max_std')

    # find top genes by MAD - all genes included
    mad = process.median_absolute_deviation(data, axis=1)
    top_genes = mad.sort_values(ascending=False).index[:n_genes]
    print "%d / %d genes (nanostring)" % (len(top_genes.intersection(all_nstring)), len(all_nstring))
    print "%d / %d genes (northcott)" % (len(top_genes.intersection(all_ncott)), len(all_ncott))

    raise ValueError()

    # Zhao data
    zhao_sample_names = (
        # 'Pt1299',
        # 'Pt1487',
        # 'Pt1595',
        'ICb1299-III',
        'ICb1299-IV',
        'ICb1487-I',
        'ICb1487-III',
        'ICb1595-I',
        'ICb1595-III',
    )
    data_zhao, meta_zhao = microarray_data.load_annotated_gse28192(
        aggr_field='SYMBOL',
        aggr_method='max_std',
        log2=True,
        sample_names=zhao_sample_names
    )
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

    # unsupervised hierarchical clustering
    z = hierarchy.linkage(X.transpose(), method='average', metric='correlation')

    den = hierarchy.dendrogram(z)

    # plot: key gene expression plus clustering
    alternatives = {
        'EGFL11': 'EYS'
    }
    g = []
    gcount = collections.Counter()

    for cls, gs in geneset:
        for gg in gs:
            if gg in data.index:
                g.append(gg)
                gcount[cls] += 1
            elif gg in alternatives and alternatives[gg] in data.index:
                g.append(alternatives[gg])
                gcount[cls] += 1

    expr = data.loc[g]

    # build row colours Series
    colour_cycle = {
        'WNT': '#2D438E',
        'SHH': '#E5161C',
        'Group C': '#F2EA00',
        'Group D': '#2A8C43'
    }

    rc_data = []
    # iterate over the geneset again to ensure classes are in the same order
    for cls, _ in geneset:
        rc_data.extend([colour_cycle[cls]] * gcount[cls])
    row_colors = pd.Series(rc_data, index=g, name='')

    cg = sns.clustermap(
        expr,
        cmap='RdBu_r',
        row_colors=row_colors,
        row_cluster=None,
        col_linkage=z,
        z_score=0,
        xticklabels=False,
    )
    cg.ax_heatmap.yaxis.label.set_visible(False)
    if show_gene_labels:
        plt.setp(
            cg.ax_heatmap.yaxis.get_ticklabels(),
            rotation=0,
            fontsize=14
        )
    else:
        cg.ax_heatmap.yaxis.set_ticklabels([])


    # now we assume that 4 clusters have been found and use those for colours
    lbl = hierarchy.fcluster(z, 4, criterion='maxclust')

    col_colors = pd.DataFrame(index=expr.columns, columns=['inferred', 'labelled'])
    col_colors.loc[lbl == 1, 'inferred'] = '#000000'
    col_colors.loc[lbl == 2, 'inferred'] = '#a5a5a5'
    col_colors.loc[lbl == 3, 'inferred'] = '#6d6d6d'
    col_colors.loc[lbl == 4, 'inferred'] = '#dbdbdb'
    for cls in colour_cycle:
        col_colors.loc[meta.subgroup == cls, 'labelled'] = colour_cycle[cls]

    cg = sns.clustermap(
        expr,
        cmap='RdBu_r',
        row_colors=row_colors,
        col_colors=col_colors,
        row_cluster=None,
        col_linkage=z,
        z_score=0,
        xticklabels=False,
    )
    # reduce whitespace
    cg.gs.update(bottom=0.02, top=0.98, left=0.02)
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

    cg = sns.clustermap(
        X,
        cmap='RdBu_r',
        col_colors=col_colors,
        col_linkage=z,
        z_score=0,
        xticklabels=False,
    )
    # remove useless row dendrogram
    cg.ax_row_dendrogram.set_visible(False)
    # reduce whitespace
    cg.gs.update(bottom=0.02, top=0.98, left=0.02)
    cg.ax_heatmap.yaxis.label.set_visible(False)
    cg.ax_heatmap.xaxis.label.set_visible(False)
    cg.ax_heatmap.yaxis.set_ticklabels([])

    # integrate the Zhao data
    # only keep common genes BUT use the same MAD genes as before

    meta_int = pd.concat((meta, mz), axis=0)
    data_int = pd.concat(
        (data.loc[common_genes, :], data_zhao.loc[common_genes, :]),
        axis=1
    )
    X_int = data_int.loc[top_genes]
    expr_int = data_int.loc[g]

    z_int = hierarchy.linkage(X_int.transpose(), method='average', metric='correlation')

    col_colors_int = pd.DataFrame(index=data_int.columns, columns=['study', 'labelled'])

    col_colors_int.loc[col_colors_int.study == 'Robinson'] = '#dbdbdb'
    col_colors_int.loc[col_colors_int.study == 'Zhao'] = '#000000'
    for cls in colour_cycle:
        col_colors_int.loc[meta_int.subgroup == cls, 'labelled'] = colour_cycle[cls]

    cg = sns.clustermap(
        expr_int,
        cmap='RdBu_r',
        row_colors=row_colors,
        row_cluster=None,
        col_colors=col_colors_int,
        col_linkage=z_int,
        # standard_scale=0,
        z_score=0,
        xticklabels=False,
    )
    # remove useless row dendrogram
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_heatmap.yaxis.label.set_visible(False)
    if show_gene_labels:
        plt.setp(
            cg.ax_heatmap.yaxis.get_ticklabels(),
            rotation=0,
            fontsize=14
        )
    else:
        cg.ax_heatmap.yaxis.set_ticklabels([])

    cg = sns.clustermap(
        X_int,
        cmap='RdBu_r',
        col_colors=col_colors_int,
        col_linkage=z_int,
        z_score=0,
        xticklabels=False,
    )
    # remove useless row dendrogram
    cg.ax_row_dendrogram.set_visible(False)
    # reduce whitespace
    cg.gs.update(bottom=0.02, top=0.98, left=0.02)
    cg.ax_heatmap.yaxis.label.set_visible(False)
    cg.ax_heatmap.xaxis.label.set_visible(False)
    cg.ax_heatmap.yaxis.set_ticklabels([])
