from matplotlib import pyplot as plt
from plotting import heatmap
import pandas as pd
from references import known_genes
from load_data import microarray_data
from utils.output import unique_output_dir
from scripts.agdex_mouse_human_mb_microarray import generate_ortholog_table
from scripts.comparison_rnaseq_microarray.consts import NORTHCOTT_GENEID, NORTHCOTT_GENEID_MAP
from microarray import process
import os
import numpy as np
import seaborn as sns
from scipy.cluster import hierarchy


def mo_entrez_to_symbol(data, kg):
    id_to_sym = kg.loc[:, ['GeneID', 'Symbol']].set_index('GeneID').dropna()
    data_sym = data.loc[data.index.intersection(id_to_sym.index)].copy()
    data_sym.set_index(id_to_sym.loc[data_sym.index, 'Symbol'], inplace=True)
    return data_sym


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


if __name__ == '__main__':

    SAVE_PLOTS = True

    if SAVE_PLOTS:
        OUTDIR = unique_output_dir('mouse_ncott_ge')

    # ID <-> symbol translation
    df = known_genes(tax_id=10090)
    id_to_sym = df.loc[:, ['GeneID', 'Symbol']].set_index('GeneID').dropna()

    # load healthy mouse cerebellum data
    mo_he, meta_he = microarray_data.load_annotated_microarray_gse54650(aggr_field='ENTREZID', aggr_method='max_std')  # indexed by Entrez gene ID

    # load mouse MB data
    mo_mb, chd7 = microarray_data.load_annotated_microarray_sb_data(aggr_method='max_std')  # indexed by Entrez gene ID

    # reduce to common genes
    common_genes = mo_he.index.intersection(mo_mb.index)

    # mo_mb = mo_mb.loc[common_genes]
    # mo_he = mo_he.loc[common_genes]

    # combine
    mo_all = pd.concat(
        (mo_mb.loc[common_genes], mo_he.loc[common_genes]),
        axis=1
    )

    # translate to gene symbol
    mo_he_sym = mo_entrez_to_symbol(mo_he.loc[common_genes], df)
    mo_mb_sym = mo_entrez_to_symbol(mo_mb.loc[common_genes], df)
    mo_all_sym = pd.concat((mo_mb_sym, mo_he_sym), axis=1)


    # translate northcott genes to mouse
    homol1 = generate_ortholog_table.homologs(9606, 10090, field='gene_id')
    homol2 = generate_ortholog_table.homologs(9606, 10090, field='gene_symbol')
    homol = pd.merge(homol1, homol2, left_index=True, right_index=True)

    idx1 = homol.loc[:, 'gene_id_9606'].isin(NORTHCOTT_GENEID_MAP.values())  # 95/100 matches
    idx2 = homol.loc[:, 'gene_id_10090'].isin(mo_all.index)  # 93/100 matches
    idx = idx1 & idx2

    ncott_mo_map = homol.loc[idx].set_index('gene_id_9606')
    ncott_geneid_mo = []
    ncott_genesym_mo = []
    for grp, arr in NORTHCOTT_GENEID:
        mo_id = np.unique(ncott_mo_map.loc[arr, 'gene_id_10090'].dropna().astype(int).values)
        ncott_geneid_mo.append((grp, mo_id))
        mo_sym = np.unique(ncott_mo_map.loc[arr, 'gene_symbol_10090'].dropna().values)
        ncott_genesym_mo.append((grp, mo_sym))

    # create df containing translated mouse gene expression
    ncott_id_all_mo = []
    for grp, arr in ncott_geneid_mo:
        ncott_id_all_mo.extend(arr)

    # add housekeeping group for reference
    ncott_genesym_mo.append(
        ('Housekeeping', [
            # 'Rn18s', # missing
            'Actb',
            # 'Gapdh', # missing
            # 'Rpl13a', # missing
            'B2m',
            'Hmbs',
            'Pgk1',
            'Hsp90ab1',
            'Hprt',
        ])
    )

    # standardise
    mo_all_sym_n = mo_all_sym.subtract(mo_he_sym.mean(axis=1), axis=0).divide(mo_he_sym.std(axis=1), axis=0)
    # mo_all_sym_n = mo_all_sym.subtract(mo_all_sym.mean(axis=1), axis=0).divide(mo_all_sym.std(axis=1), axis=0)
    mo_all_sym_yg = process.yugene_transform(mo_all_sym)

    # distribution of all intensities
    if SAVE_PLOTS:
        # hist of intensities
        fig, axs = plt.subplots(2, 1, sharex=True)
        axs[0].hist(mo_mb.values.flatten(), 100)
        axs[0].set_title("Microarray intensities, SB")
        axs[1].hist(mo_he.values.flatten(), 100)
        axs[1].set_title("Microarray intensities, cerebellum")
        plt.tight_layout()
        fig.savefig(os.path.join(OUTDIR, "marr_sb-circad_hist_intensity.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "marr_sb-circad_hist_intensity.pdf"))

    if SAVE_PLOTS:
        # Plot: Mouse SB vs healthy circadian cerebellum, all Ncott genes standardised
        fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
            ncott_genesym_mo,
            mo_all_sym_n,
            vmax=10.,
            vmin=-10.,
            cbar=True,
            orientation='vertical',
            fig_kwargs={'figsize': [7, 11]},
            heatmap_kwargs={'square': False},
            gs_kwargs={'left': 0.25}
        )
        # reduce y label font size
        for ax in axs:
            plt.setp(ax.yaxis.get_ticklabels(), fontsize=8.5)
        # add dividing lines
        xbreaks = [3, 8]
        for ax in axs:
            for t in xbreaks:
                ax.axvline(t, color='w', linewidth=2.5)
                ax.axvline(t, color='0.4', linewidth=1.)
        fig.savefig(os.path.join(OUTDIR, "marr_sb-circad_ncott_homol_standardised.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "marr_sb-circad_ncott_homol_standardised.pdf"))

    if SAVE_PLOTS:
        # Plot: Mouse SB vs healthy circadian cerebellum, all Ncott genes
        fig, axs, cax, gs = heatmap.grouped_expression_heatmap(
            ncott_genesym_mo,
            mo_all_sym_yg,
            vmax=1.,
            vmin=0.,
            cbar=True,
            orientation='vertical',
            fig_kwargs={'figsize': [7, 11]},
            heatmap_kwargs={'square': False},
            gs_kwargs={'left': 0.25}
        )
        # reduce y label font size
        for ax in axs:
            plt.setp(ax.yaxis.get_ticklabels(), fontsize=8.5)
        # add dividing lines
        xbreaks = [3, 8]
        for ax in axs:
            for t in xbreaks:
                ax.axvline(t, color='w', linewidth=2.5)
                ax.axvline(t, color='0.4', linewidth=1.)
        fig.savefig(os.path.join(OUTDIR, "marr_sb-circad_ncott_homol_yg.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "marr_sb-circad_ncott_homol_yg.pdf"))

    # plot: clustered heatmap, Ncott + hkeeping only
    gg = []
    for _, arr in ncott_genesym_mo:
        gg.extend(arr)
    mo_ncott_sym = mo_all_sym.loc[gg]
    z = hierarchy.linkage(mo_ncott_sym.transpose(), method='average', metric='correlation')

    col_colors = pd.DataFrame(index=mo_ncott_sym.columns, columns=['chd7_status',])
    cluster_colours = {
        -1: '#000000',
        0: '#a5a5a5',
        1: '#dbdbdb',
    }

    a = np.zeros(len(chd7))
    a[chd7] = 1
    a = np.concatenate((a, -np.ones(meta_he.shape[0])))

    col_colors.loc[:, 'chd7_status'] = [cluster_colours.get(t) for t in a]

    cg = plot_clustermap(
        mo_ncott_sym,
        cmap='RdBu_r',
        col_colors=col_colors,
        col_linkage=z,
        z_score=0,
        xticklabels=True,
        show_gene_labels=True,
        row_cluster=None,
    )
    plt.setp(cg.ax_heatmap.yaxis.get_ticklabels(), fontsize=7)
    if SAVE_PLOTS:
        cg.savefig(os.path.join(OUTDIR, "clustermap_ncott.png"), dpi=200)
        cg.savefig(os.path.join(OUTDIR, "clustermap_ncott.pdf"))

    mo_ncott_sym_yg = mo_all_sym_yg.loc[gg]
    z_yg = hierarchy.linkage(mo_ncott_sym_yg.transpose(), method='average', metric='correlation')

    cg = plot_clustermap(
        mo_ncott_sym_yg,
        cmap='RdBu_r',
        col_colors=col_colors,
        col_linkage=z_yg,
        z_score=0,
        xticklabels=True,
        show_gene_labels=True
    )
    plt.setp(cg.ax_heatmap.yaxis.get_ticklabels(), fontsize=7)