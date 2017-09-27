from load_data import microarray_data, allen_human_brain_atlas, rnaseq_data
from plotting.pca import pca_plot_by_group_2d, pca_plot_by_group_3d, cluster_ellipsoid
from analysis import process
from utils.output import unique_output_dir

from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os
import collections
from scipy.cluster import hierarchy

from scripts.mb_subgroup_classifier import load
## TODO: rename this:
from scripts.comparison_rnaseq_microarray import consts


def combine_expr(*args):
    common_genes = args[0].index
    sample_names = args[0].columns
    for t in args[1:]:
        common_genes = common_genes.intersection(t.index)
        sample_names = sample_names.append(t.columns)
    res = pd.DataFrame(index=common_genes, columns=sample_names)
    for t in args:
        res.loc[:, t.columns] = t.loc[common_genes]
    return res


def plot_2d(y, lbl, colour_map, marker_map, title=None, outdir=None, additional_data=None):
    # plots: PCA of classifier vs RNA-Seq
    fig, axs = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=False, figsize=(12, 5), num=title)
    for i, compo in enumerate([(0, 1), (1, 2), (0, 2)]):
        pca_plot_by_group_2d(
            y, lbl, components=compo, colour_map=colour_map,
            marker_map=marker_map,
            ax=axs[i],
            legend=False,
            additional_data_dict=additional_data
        )
    axs[-1].legend(loc='upper right')
    plt.tight_layout(pad=0.2, rect=[.02, .02, 1, 1])

    if title and outdir:
        plt.draw()
        fig.savefig(os.path.join(outdir, "%s.png" % title), dpi=300)
        fig.savefig(os.path.join(outdir, "%s.tiff" % title), dpi=200)
        fig.savefig(os.path.join(outdir, "%s.pdf" % title))

    return fig, axs


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


if __name__ == "__main__":

    N_PC = 3
    geneset = consts.NORTHCOTT_GENES
    outdir = unique_output_dir("pca_atcc_lines")

    # it's useful to maintain a list of known upregulated genes
    nano_genes = []
    for grp, arr in consts.NANOSTRING_GENES:
        if grp != 'WNT':
            nano_genes.extend(arr)
    nano_genes.remove('EGFL11')
    nano_genes.append('EYS')

    all_nstring = []
    [all_nstring.extend(t) for _, t in consts.NANOSTRING_GENES]
    all_ncott = []
    [all_ncott.extend(t) for _, t in consts.NORTHCOTT_GENES]


    # load Ncott data (285 non-WNT MB samples)
    # ncott, ncott_meta = microarray_data.load_annotated_microarray_gse37382(
    #     aggr_field='SYMBOL',
    #     aggr_method='max'
    # )
    # sort_idx = ncott_meta.subgroup.sort_values().index
    # ncott_meta = ncott_meta.loc[sort_idx]
    # ncott = ncott.loc[:, sort_idx]

    # load Kool dataset
    # kool, kool_meta = microarray_data.load_annotated_microarray_gse10327(
    #     aggr_field='SYMBOL',
    #     aggr_method='max',
    # )
    # sort_idx = kool_meta.subgroup.sort_values().index
    # kool_meta = kool_meta.loc[sort_idx]
    # kool = kool.loc[:, sort_idx]
    # kool_meta.loc[:, 'subgroup'] = (
    #     kool_meta.loc[:, 'subgroup'].str
    #         .replace('A', 'WNT')
    #         .replace('B', 'SHH')
    #         .replace('E', 'Group 3')
    #         .replace('C', 'Group 4')
    #         .replace('D', 'Group 4')
    # )

    # load Robinson dataset
    STUDY = 'Robinson'
    dat_ref, meta_ref = microarray_data.load_annotated_microarray_gse37418(aggr_field='SYMBOL', aggr_method='max')
    meta_ref = meta_ref.loc[~meta_ref.subgroup.isin(['U', 'SHH OUTLIER'])]
    sort_idx = meta_ref.subgroup.sort_values().index
    meta_ref = meta_ref.loc[sort_idx]
    dat_ref = dat_ref.loc[:, sort_idx]
    meta_ref.loc[:, 'subgroup'] = meta_ref.subgroup.str.replace('G3', 'Group 3').replace('G4', 'Group 4')
    meta_ref.loc[:, 'study'] = STUDY

    # load XZ RNA-Seq count data
    ## TODO: replace with a newer loader
    dat_xz = load.load_xz_rnaseq(kind='htseq', yugene=False)
    dat_xz = dat_xz.loc[:, ['XZ1', 'XZ2']]
    dat_xz.columns = ['Zhang 1299', 'Zhang 1299 (rpt)']
    meta_xz = pd.DataFrame([['Group 4']] * 2, index=dat_xz.columns, columns=['subgroup'])
    meta_xz.loc[:, 'study'] = 'Zhang'

    # load SB RNA-Seq count data
    # NB: have checked and using TPM rather than FPKM makes no difference, as expected
    obj_sb = rnaseq_data.zhao_mb_cultures(annotate_by='Approved Symbol')
    dat_sb = obj_sb.data.loc[:, ['ICb1595']]
    dat_sb.columns = ['Zhang 1595']
    meta_sb = pd.DataFrame([['Group 3']], index=dat_sb.columns, columns=['subgroup'])
    meta_sb.loc[:, 'study'] = 'Zhang'

    # load ATCC cell line data
    obj_atcc = rnaseq_data.atcc_cell_lines(annotate_by='Approved Symbol')
    dat_atcc = obj_atcc.data
    dat_atcc.columns = ['CRL3021']
    meta_atcc = pd.DataFrame([['Group 4']], index=dat_atcc.columns, columns=['subgroup'])
    meta_atcc.loc[:, 'study'] = 'Zhang'

    # load Xiao-Nan data
    xnan_sample_names = (
        'Pt1299',
        'ICb1299-I',
        'ICb1299-III',
        'ICb1299-IV',
        # 'ICb1487-I',
        # 'ICb1487-III',
        'Pt1595',
        'ICb1595-I',
        'ICb1595-III',
    )
    dat_xnan, meta_xnan = load.load_xiaonan_microarray(yugene=False, sample_names=xnan_sample_names)
    meta_xnan = pd.DataFrame(
        meta_xnan.loc[dat_xnan.columns, 'northcott classification'].replace('C', 'Group 3').replace('D', 'Group 4').values,
        index=dat_xnan.columns,
        columns=['subgroup']
    )
    meta_xnan.loc[:, 'study'] = 'Zhao'

    ###################################
    # extract only the desired samples
    # can switch classifier here.
    ###################################

    title = 'robi'
    all_data = combine_expr(dat_ref, dat_xz, dat_sb, dat_atcc, dat_xnan)

    xnan_sample_names1 = [
        'ICb1299-III',
        'ICb1299-IV',
        'ICb1595-I',
        'ICb1595-III',
    ]

    all_data = combine_expr(dat_ref, dat_xz, dat_sb, dat_atcc, dat_xnan.loc[:, xnan_sample_names1])

    all_data_yg = process.yugene_transform(all_data, resolve_ties=False)
    X = all_data_yg.loc[:, meta_ref.index].transpose()
    m = meta_ref.copy()

    # first fit the pca with lots of components to investigate the explained variance
    pca = PCA(n_components=10)
    pca.fit(X)
    y = pca.transform(X)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(range(1, 11), np.cumsum(pca.explained_variance_), '-o')
    plt.axis([0.8, 10.2, 0, 100.])
    ax.set_xticks(range(1, 11))
    plt.xlabel('Principal component')
    plt.ylabel('Cumulative % variance explained')
    fig.savefig(os.path.join(outdir, "pca_variance_explained.pdf"), dpi=200)
    fig.savefig(os.path.join(outdir, "pca_variance_explained.png"), dpi=200)

    pca = PCA(n_components=N_PC)
    pca.fit(X)
    y = pca.transform(X)

    y_xz = pca.transform(all_data_yg.loc[:, dat_xz.columns].transpose())  # only scramble
    y_sb = pca.transform(all_data_yg.loc[:, dat_sb.columns].transpose())  # 1595 only
    y_atcc = pca.transform(all_data_yg.loc[:, dat_atcc.columns].transpose())  # primary only
    y_xnan = pca.transform(all_data_yg.loc[:, xnan_sample_names1].transpose())

    y_1299 = y_xnan[0:2]
    y_1595 = y_xnan[2:]

    # get labels and plot by subgroup
    idx, labels = m.subgroup.factorize()

    # define colours and labels
    lbl_1299 = 'ICb1299 (Zhao et al.)'
    lbl_1595 = 'ICb1595 (Zhao et al.)'
    lbl_atcc = 'CHLA-01-Med'
    lbl_xz1299 = 'ICb1299'
    lbl_xz1595 = 'ICb1595'

    colour_map = {
        'Group 3': '#F2EA00',
        'Group 4': '#2A8C43',
        'WNT': '#2D438E',
        'SHH': '#E5161C',
        'control': 'gray',
        lbl_1299: 'c',
        lbl_xz1299: 'c',
        lbl_1595: 'k',
        lbl_xz1595: 'k',
        lbl_atcc: 'm'
    }

    marker_map = dict([(k, 'o') for k in colour_map])
    marker_map[lbl_1299] = 's'
    marker_map[lbl_xz1299] = 's'
    marker_map[lbl_1595] = '^'
    marker_map[lbl_xz1595] = '^'
    marker_map[lbl_atcc] = 'D'



    # plots: PCA of classifier vs RNA-Seq
    ttl = ("pca_%s-rnaseq_2d" % title)
    ad = collections.OrderedDict([
        (lbl_1299, y_1299),
        (lbl_1595, y_1595),
        (lbl_atcc, y_atcc),
    ])
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, outdir=outdir, additional_data=ad)

    # version 2: switch out Zhao microarray data for Zhang RNA-Seq data
    ttl = ("pca_%s-rnaseq-v2" % title)
    ad = collections.OrderedDict([
        (lbl_xz1299, y_xz),
        (lbl_xz1595, y_sb),
        (lbl_atcc, y_atcc),
    ])
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, outdir=outdir, additional_data=ad)

    # clusterplot

    # we assume that 4 clusters have been found and distinguish by greyscale shade
    cluster_colours = {
        1: '#000000',
        2: '#a5a5a5',
        3: '#6d6d6d',
        4: '#dbdbdb',
    }

    study_colours = {
        STUDY: '#cccccc',
        'Zhao': '#000000',
        'Zhang': '#666666',
    }

    ###################################
    # extract only the desired samples
    # can switch classifier here.
    ###################################

    title = 'robi'

    # prepare data

    xnan_sample_names2 = [
        'Pt1299',
        'ICb1299-I',
        'ICb1299-III',
        'ICb1299-IV',
        'Pt1595',
        'ICb1595-I',
        'ICb1595-III',
    ]

    xz_sample_names = ['Zhang 1299']

    data_for_clustering = combine_expr(dat_ref, dat_xz.loc[:, xz_sample_names], dat_sb, dat_atcc, dat_xnan.loc[:, xnan_sample_names2])
    meta_for_clustering = pd.concat((meta_ref, meta_xz.loc[xz_sample_names], meta_sb, meta_atcc, meta_xnan.loc[xnan_sample_names2]), axis=0)

    # rename CRL3021 -> CHLA-01-Med
    meta_for_clustering = meta_for_clustering.loc[data_for_clustering.columns]
    data_col = data_for_clustering.columns.tolist()
    data_col[data_col.index('CRL3021')] = 'Zhang CHLA-01-Med'
    data_for_clustering.columns = data_col
    meta_for_clustering.index = data_col

    if 'EYS' in data_for_clustering.index:
        idx = data_for_clustering.index.str.replace('EYS', 'EGFL11')
        data_for_clustering.index = idx


    gene_class = {}
    alternatives = {
        'EGFL11': 'EYS'
    }

    g_nano = []
    gcount_nano = collections.Counter()

    # matching genes
    for cls, gs in consts.NANOSTRING_GENES:
        for gg in gs:
            # if gg in data.index:
            if gg in data_for_clustering.index:
                g_nano.append(gg)
                gcount_nano[cls] += 1
            elif gg in alternatives and alternatives[gg] in data_for_clustering.index:
                g_nano.append(alternatives[gg])
                gcount_nano[cls] += 1

    g_ncot = []

    for cls, gs in consts.NORTHCOTT_GENES:
        for gg in gs:
            if cls == 'Group C':
                cls = 'Group 3'
            if cls == 'Group D':
                cls = 'Group 4'
            gene_class[gg] = cls
            # if gg in data.index:
            if gg in data_for_clustering.index:
                g_ncot.append(gg)
            elif gg in alternatives and alternatives[gg] in data_for_clustering.index:
                g_ncot.append(alternatives[gg])

    data_for_clustering_ncot = data_for_clustering.loc[g_ncot].dropna().astype(float)
    data_for_clustering_ncot = np.log2(data_for_clustering_ncot + 1)
    z_ncot = hierarchy.linkage(data_for_clustering_ncot.transpose(), method='average', metric='correlation')

    row_colors_nano = pd.Series(dict([(k, gene_class[k]) for k in g_nano]), name='').apply(colour_map.get)
    row_colors_ncot = pd.Series(gene_class, name='').apply(colour_map.get)

    col_colors = pd.DataFrame(index=data_for_clustering_ncot.columns, columns=['study', 'subgroup'])
    col_colors_subgroup = meta_for_clustering.subgroup.apply(colour_map.get)
    col_colors.loc[col_colors_subgroup.index, 'subgroup'] = col_colors_subgroup
    col_colors.loc[:, 'study'] = meta_for_clustering.study.apply(study_colours.get)
    # col_colors.loc[:, 'inferred'] = pd.Series(
    #     hierarchy.fcluster(z_ncot, 4, criterion='maxclust'),
    #     index=data_for_clustering_ncot.columns
    # ).apply(cluster_colours.get)

    cg = plot_clustermap(
        data_for_clustering_ncot,
        cmap='RdBu_r',
        show_gene_labels=True,
        row_colors=row_colors_ncot,
        col_colors=col_colors,
        row_cluster=None,
        col_linkage=z_ncot,
        z_score=1,
        xticklabels=True,
    )
    plt.setp(
        cg.ax_heatmap.yaxis.get_ticklabels(),
        rotation=0,
        fontsize=8
    )
    # clustering.show_dendrogram_cut_by_nclust(cg, 3)
    # clustering.show_dendrogram_cut_by_nclust(cg, 4)
    # clustering.show_dendrogram_cut_by_nclust(cg, 5)
    ttl = "%s_ncott_heatmap" % STUDY.lower()
    cg.savefig(os.path.join(outdir, "%s.png" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.tiff" % ttl), dpi=200)
    cg.savefig(os.path.join(outdir, "%s.pdf" % ttl))