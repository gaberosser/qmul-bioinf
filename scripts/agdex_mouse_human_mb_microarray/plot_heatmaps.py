from matplotlib import pyplot as plt
from plotting import heatmap
import pandas as pd
from references import known_genes
from scripts.agdex_mouse_human_mb_microarray import generate_ortholog_table, load_data
from scripts.comparison_rnaseq_microarray.consts import NORTHCOTT_GENEID, NORTHCOTT_GENEID_MAP
from microarray import process
import os
import numpy as np


def mo_entrez_to_symbol(data, kg):
    id_to_sym = kg.loc[:, ['GeneID', 'Symbol']].set_index('GeneID').dropna()
    data_sym = data.loc[data.index.intersection(id_to_sym.index)].copy()
    data_sym.set_index(id_to_sym.loc[data_sym.index, 'Symbol'], inplace=True)
    return data_sym


if __name__ == '__main__':

    SAVE_PLOTS = False

    if SAVE_PLOTS:
        OUTDIR = 'mouse_ncott_ge.0'
        i = 1
        while os.path.exists(OUTDIR):
            OUTDIR = 'mouse_ncott_ge.%d' % i
            i += 1
        print "Creating temp output dir %s" % OUTDIR
        os.makedirs(OUTDIR)

    # ID <-> symbol translation
    df = known_genes(tax_id=10090)
    id_to_sym = df.loc[:, ['GeneID', 'Symbol']].set_index('GeneID').dropna()

    # load healthy mouse cerebellum data
    mo_he = load_data.load_annotated_microarray_gse54650()  # indexed by Entrez gene ID

    # load mouse MB data
    mo_mb, chd7 = load_data.load_annotated_microarray_sb_data()  # indexed by Entrez gene ID

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