from load_data import microarray_data, allen_human_brain_atlas, rnaseq_data
from plotting.pca import pca_plot_by_group_2d, pca_plot_by_group_3d, cluster_ellipsoid
from scripts.comparison_rnaseq_microarray import consts
from analysis import process
from utils.output import unique_output_dir

from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os

from scripts.mb_subgroup_classifier import load


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


def plot_2d(y, lbl, colour_map, marker_map, title=None, outdir=None, **additional_data):
    # plots: PCA of classifier vs RNA-Seq
    fig, axs = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=False, figsize=(12, 5), num=title)
    for i, compo in enumerate([(0, 1), (1, 2), (0, 2)]):
        pca_plot_by_group_2d(
            y, lbl, components=compo, colour_map=colour_map,
            marker_map=marker_map,
            ax=axs[i],
            legend=False,
            **additional_data
        )
    axs[-1].legend(loc='upper right')
    plt.tight_layout(pad=0.2, rect=[.02, .02, 1, 1])
    if title and outdir:
        fig.savefig(os.path.join(outdir, "%s.png" % title), dpi=300)
        fig.savefig(os.path.join(outdir, "%s.tiff" % title), dpi=200)
        fig.savefig(os.path.join(outdir, "%s.pdf" % title))

    return fig, axs


if __name__ == "__main__":

    N_PC = 3
    outdir = unique_output_dir("pca_atcc_lines")

    # it's useful to maintain a list of known upregulated genes
    nano_genes = []
    for grp, arr in consts.NANOSTRING_GENES:
        if grp != 'WNT':
            nano_genes.extend(arr)
    nano_genes.remove('EGFL11')
    nano_genes.append('EYS')

    # load Ncott data (285 non-WNT MB samples)
    ncott, ncott_meta = microarray_data.load_annotated_microarray_gse37382(
        aggr_field='SYMBOL',
        aggr_method='max'
    )
    sort_idx = ncott_meta.subgroup.sort_values().index
    ncott_meta = ncott_meta.loc[sort_idx]
    ncott = ncott.loc[:, sort_idx]
    ncott = process.yugene_transform(ncott)

    # load Allen (healthy cerebellum)

    he, he_meta = allen_human_brain_atlas.cerebellum_microarray_reference_data(agg_field='gene_symbol', agg_method='max')
    he_meta.loc[:, 'subgroup'] = 'control'

    # load Kool dataset
    kool, kool_meta = microarray_data.load_annotated_microarray_gse10327(
        aggr_field='SYMBOL',
        aggr_method='max',
    )
    sort_idx = kool_meta.subgroup.sort_values().index
    kool_meta = kool_meta.loc[sort_idx]
    kool = kool.loc[:, sort_idx]
    kool_meta.loc[:, 'subgroup'] = (
        kool_meta.loc[:, 'subgroup'].str
            .replace('A', 'WNT')
            .replace('B', 'SHH')
            .replace('E', 'Group 3')
            .replace('C', 'Group 4')
            .replace('D', 'Group 4')
    )
    kool = process.yugene_transform(kool)

    # load Robinson dataset
    robi, robi_meta = microarray_data.load_annotated_microarray_gse37418(aggr_field='SYMBOL', aggr_method='max')
    robi_meta = robi_meta.loc[~robi_meta.subgroup.isin(['U', 'SHH OUTLIER'])]
    sort_idx = robi_meta.subgroup.sort_values().index
    robi_meta = robi_meta.loc[sort_idx]
    robi = robi.loc[:, sort_idx]
    robi_meta.loc[:, 'subgroup'] = robi_meta.subgroup.str.replace('G3', 'Group 3').replace('G4', 'Group 4')
    robi = process.yugene_transform(robi)

    # combine them all - this ensures a common list of genes
    combs = []
    grps = []
    # combs.append(he); grps.extend(he_meta.subgroup.values)
    combs.append(ncott); grps.extend(ncott_meta.subgroup.values)
    combs.append(kool); grps.extend(kool_meta.subgroup.values)
    combs.append(robi); grps.extend(robi_meta.subgroup.values)

    all_expr = combine_expr(*combs)

    # load XZ RNA-Seq count data
    ## TODO: replace with a newer loader
    dat_xz = load.load_xz_rnaseq(kind='htseq', yugene=True)
    # X_xz = load_xz_rnaseq(kind='cuff', yugene=True, gene_symbols=X.columns).transpose()

    # load SB RNA-Seq count data
    # NB: have checked and using TPM rather than FPKM makes no difference, as expected
    obj_sb = rnaseq_data.zhao_mb_cultures(annotate_by='Approved Symbol')
    dat_sb = process.yugene_transform(obj_sb.data.loc[:, ['ICb1595']])

    # load ATCC cell line data
    obj_atcc = rnaseq_data.atcc_cell_lines(annotate_by='Approved Symbol')
    dat_atcc = process.yugene_transform(obj_atcc.data)

    all_data = combine_expr(*[all_expr, dat_xz, dat_sb, dat_atcc])


    ###################################
    # extract only the desired samples
    # can switch classifier here.
    ###################################
    title = 'robi'
    X = all_data.loc[:, robi_meta.index].transpose()
    m = robi_meta.copy()

    # title = 'kool'
    # X = all_expr.loc[:, kool_meta.index].transpose()
    # m = kool_meta.copy()

    # title = 'northcott'
    # X = all_expr.loc[:, ncott_meta.index].transpose()
    # m = ncott_meta.copy()

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

    y_xz = pca.transform(all_data.loc[:, ['XZ1', 'XZ2']].transpose())  # only scramble
    y_sb = pca.transform(all_data.loc[:, dat_sb.columns].transpose())  # 1595 only
    y_atcc = pca.transform(all_data.loc[:, dat_atcc.columns].transpose())  # primary only

    # load Xiao-Nan data
    xnan_sample_names = (
        'ICb1299-III',
        'ICb1299-IV',
        'ICb1487-I',
        'ICb1487-III',
        'ICb1595-I',
        'ICb1595-III',
    )
    xnan, xnan_meta = load.load_xiaonan_microarray(yugene=True, gene_symbols=X.columns, sample_names=xnan_sample_names)
    xnan = xnan.transpose()
    y_xnan = pca.transform(xnan)

    y_1299 = y_xnan[0:2]
    y_1487 = y_xnan[2:4]
    y_1595 = y_xnan[4:]

    # get labels and plot by subgroup
    idx, labels = m.subgroup.factorize()

    # define colours and labels
    lbl_1299_this = 'ICb1299 (XZ)'
    lbl_1595_this = 'ICb1595 (SB)'

    colour_map = {
        'Group 3': '#F2EA00',
        'Group 4': '#2A8C43',
        'WNT': '#2D438E',
        'SHH': '#E5161C',
        'control': 'gray',
        lbl_1299_this: 'b',
        lbl_1595_this: 'y',
        'ATCC': 'k'
    }

    marker_map = dict([(k, 'o') for k in colour_map])
    marker_map[lbl_1299_this] = 's'
    marker_map[lbl_1595_this] = 's'
    marker_map['ATCC'] = 's'


    # plots: PCA of classifier vs RNA-Seq
    ttl = ("pca_%s-rnaseq_2d" % title)
    ad = {
        lbl_1299_this: y_xz,
        lbl_1595_this: y_sb,
        'ATCC': y_atcc,
    }
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)
