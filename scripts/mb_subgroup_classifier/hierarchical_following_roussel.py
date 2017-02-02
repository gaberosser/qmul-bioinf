import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from load_data import microarray_data, rnaseq_data
from microarray import process
from scipy.cluster import hierarchy
from scripts.comparison_rnaseq_microarray import consts
import collections
from scripts.output import unique_output_dir
import os


def plot_clustermap(dat, show_gene_labels=False, **kwargs):
    cg = sns.clustermap(
        dat,
        **kwargs
    )
    # remove useless row dendrogram
    cg.ax_row_dendrogram.set_visible(False)
    # and remove the space created for it
    wr = cg.gs.get_width_ratios()
    wr[0] = 0.05
    wr[1] = 0.02
    cg.gs.set_width_ratios(wr)
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
    return cg


if __name__ == '__main__':
    SAVE_FIG = True
    geneset = consts.NORTHCOTT_GENES
    show_gene_labels = False
    n_genes = 1500
    AGGR_METHOD = 'max_std'

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
    data = data.loc[:, meta.index]

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
    data_zhao, meta_zhao = microarray_data.load_annotated_gse28192(
        sample_names=zhao_sample_names,
        log2=True,
        aggr_field='SYMBOL',
        aggr_method=AGGR_METHOD
    )
    # data_zhao, meta_zhao = microarray_data.load_annotated_gse28192(
    #     sample_names=zhao_sample_names,
    # )
    # data_zhao = process.variance_stabilizing_transform(data_zhao)
    # data_zhao = data_zhao.loc[:, zhao_sample_names]
    # meta_zhao = meta_zhao.loc[zhao_sample_names, :]
    # data_zhao = microarray_data.annotate_and_aggregate_gse28192(data_zhao, aggr_field='SYMBOL', aggr_method=AGGR_METHOD)


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
    lbl = hierarchy.fcluster(z, 4, criterion='maxclust')

    col_colors = pd.DataFrame(index=X.columns, columns=['inferred', 'labelled'])
    col_colors.loc[lbl == 1, 'inferred'] = '#000000'
    col_colors.loc[lbl == 2, 'inferred'] = '#a5a5a5'
    col_colors.loc[lbl == 3, 'inferred'] = '#6d6d6d'
    col_colors.loc[lbl == 4, 'inferred'] = '#dbdbdb'
    for cls in colour_cycle:
        col_colors.loc[meta.subgroup == cls, 'labelled'] = colour_cycle[cls]

    # Nanostring heatmap
    cg = plot_clustermap(
        expr_nano,
        show_gene_labels=True,
        cmap='RdBu_r',
        row_colors=row_colors_nano,
        col_colors=col_colors,
        row_cluster=None,
        col_linkage=z,
        z_score=0,
        xticklabels=False,
    )
    ttl = "%s_nano_heatmap" % STUDY.lower()
    cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    # Northcott heatmap
    cg = plot_clustermap(
        expr_ncot,
        cmap='RdBu_r',
        row_colors=row_colors_ncot,
        col_colors=col_colors,
        row_cluster=None,
        col_linkage=z,
        z_score=0,
        xticklabels=False,
    )
    ttl = "%s_ncott_heatmap" % STUDY.lower()
    cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    # Global (n_genes) heatmap
    cg = plot_clustermap(
        X,
        cmap='RdBu_r',
        col_colors=col_colors,
        col_linkage=z,
        z_score=0,
        xticklabels=False,
    )
    ttl = "%s_top%d_heatmap" % (STUDY.lower(), n_genes)
    cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    # integrate the Zhao data
    # only keep common genes BUT use the same MAD genes as before

    meta_int = pd.concat((meta, mz), axis=0)
    data_int = pd.concat(
        (data.loc[common_genes, :], data_zhao.loc[common_genes, :]),
        axis=1
    )
    X_int = data_int.loc[top_genes]
    expr_int_nano = data_int.loc[g_nano]
    expr_int_ncot = data_int.loc[g_ncot]

    z_int = hierarchy.linkage(X_int.transpose(), method='average', metric='correlation')

    col_colors_int = pd.DataFrame(index=data_int.columns, columns=['study', 'labelled'])

    col_colors_int.loc[meta_int.study == STUDY] = '#dbdbdb'
    col_colors_int.loc[meta_int.study == 'Zhao'] = '#000000'
    for cls in colour_cycle:
        col_colors_int.loc[meta_int.subgroup == cls, 'labelled'] = colour_cycle[cls]

    cg = plot_clustermap(
        expr_int_nano,
        show_gene_labels=True,
        cmap='RdBu_r',
        row_colors=row_colors_nano,
        row_cluster=None,
        col_colors=col_colors_int,
        col_linkage=z_int,
        z_score=0,
        xticklabels=False,
    )
    ttl = "%s_zhao_nano_heatmap" % STUDY.lower()
    cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    cg = plot_clustermap(
        expr_int_ncot,
        cmap='RdBu_r',
        row_colors=row_colors_ncot,
        row_cluster=None,
        col_colors=col_colors_int,
        col_linkage=z_int,
        z_score=0,
        xticklabels=False,
    )
    ttl = "%s_zhao_ncott_heatmap" % STUDY.lower()
    cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))

    cg = plot_clustermap(
        X_int,
        cmap='RdBu_r',
        col_colors=col_colors_int,
        col_linkage=z_int,
        z_score=0,
        xticklabels=False,
    )
    ttl = "%s_zhao_top%d_heatmap" % (STUDY.lower(), n_genes)
    cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))
