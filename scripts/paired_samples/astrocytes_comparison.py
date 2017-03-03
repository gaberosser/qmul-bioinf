import os
import re

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy

from load_data import rnaseq_data
from microarray import process
from utils.output import unique_output_dir
from plotting import clustering

from scripts.rnaseq import get_rrna_from_gtf

import references

if __name__ == "__main__":
    N_GENES = 1500
    OUTDIR = unique_output_dir("astrocytes", reuse_empty=True)

    # GSE73721 (reference astrocytes, oligos, ...)
    obj73721 = rnaseq_data.gse73721(source='star', annotate_by='Ensembl Gene ID')

    # GSE61794 (NSC x 2)
    obj61794 = rnaseq_data.gse61794(source='star', annotate_by='Ensembl Gene ID')

    # GBM paired samples
    objwtchg = rnaseq_data.gbm_astrocyte_nsc_samples_loader(source='star', annotate_by='Ensembl Gene ID')

    # rRNA gene IDs
    rrna_ensg = set(get_rrna_from_gtf.get_rrna())

    # combine the data
    data = pd.concat((obj73721.data, obj61794.data, objwtchg.data), axis=1)
    data = data.loc[data.index.str.contains('ENSG')]

    data_by_symbol = pd.concat((obj73721.data_by_symbol, obj61794.data_by_symbol, objwtchg.data_by_symbol), axis=1)

    mad = process.median_absolute_deviation(data, axis=1).sort_values(ascending=False)
    top_mad = mad.iloc[:N_GENES].index

    z = hierarchy.linkage(data.loc[top_mad].transpose(), method='average', metric='correlation')
    col_colors = pd.DataFrame(index=data.columns, columns=['group'])
    cmap = {
        'tumor': '#a6cee3',
        'Hippocampus': '#1f78b4',
        'Fetal': '#b2df8a',
        re.compile(r'ctx [0-9]* *astro', flags=re.IGNORECASE): '#33a02c',
        'ctx neuron': '#fb9a99',
        'oligo': '#e31a1c',
        'myeloid': '#fdbf6f',
        'ctx endo': '#ff7f00',
        'whole cortex': '#cab2d6',
        'DURA': 'black',
    }
    for t, col in cmap.items():
        col_colors.loc[col_colors.index.str.contains(t)] = col

    row_colors = pd.DataFrame(index=data.index, columns=['RNA type'])
    row_colors.loc[row_colors.index.isin(rrna_ensg)] = 'black'

    cg = clustering.plot_clustermap(
        data.loc[top_mad],
        col_linkage=z,
        z_score=0,
        show_gene_labels=False,
        col_colors=col_colors,
        row_colors=row_colors,
    )
    plt.setp(
        cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90
    )
    cg.gs.update(bottom=0.2)

    data_yg = process.yugene_transform(data)
    # reindex by gene name
    gg = references.ensembl_to_gene_symbol(data_yg.index)
    mad_yg = process.median_absolute_deviation(data_yg, axis=1).sort_values(ascending=False)
    top_mad_yg = mad_yg.iloc[:N_GENES].index

    z_yg = hierarchy.linkage(data_yg.loc[top_mad_yg].transpose(), method='average', metric='correlation')

    cg = clustering.plot_clustermap(
        data_yg.loc[top_mad_yg],
        cmap='RdBu_r',
        col_linkage=z_yg,
        col_colors=col_colors,
        row_colors=row_colors,
    )
    plt.setp(
        cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90
    )
    cg.gs.update(bottom=0.2)

    # investigate rRNA quantity
    aa = data.loc[rrna_ensg].sum() / data.sum()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = aa.plot.bar()
    ax.set_ylim([0, 1])
    ax.set_ylabel("Proportion rRNA")


    # remove rRNA
    data_rr = data.loc[~data.index.isin(rrna_ensg)]
    data_rr_yg = process.yugene_transform(data_rr)

    mad_rr = process.median_absolute_deviation(data_rr, axis=1).sort_values(ascending=False)
    top_mad_rr = mad_rr.iloc[:N_GENES].index
    mad_rr_yg = process.median_absolute_deviation(data_rr_yg, axis=1).sort_values(ascending=False)
    top_mad_rr_yg = mad_rr_yg.iloc[:N_GENES].index

    z_rr = hierarchy.linkage(data_rr.loc[top_mad_rr].transpose(), method='average', metric='correlation')
    z_rr_yg = hierarchy.linkage(data_rr_yg.loc[top_mad_rr_yg].transpose(), method='average', metric='correlation')

    cg = clustering.plot_clustermap(
        data_rr.loc[top_mad_rr],
        cmap='RdBu_r',
        z_score=0,
        col_linkage=z_yg,
        col_colors=col_colors,
    )
    plt.setp(
        cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90
    )
    cg.gs.update(bottom=0.2)

    cg = clustering.plot_clustermap(
        data_rr_yg.loc[top_mad_rr_yg],
        cmap='RdBu_r',
        col_linkage=z_yg,
        col_colors=col_colors,
    )
    plt.setp(
        cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90
    )
    cg.gs.update(bottom=0.2)


    # for reference, we can also load the original 'pre-processed' GSE73721 data
    # but these are indexed by mouse gene??
    fpkm, meta = rnaseq_data.brainrnaseq_preprocessed()
    fpkm_idx = np.array(fpkm.index.str.capitalize())
    fpkm_idx[fpkm_idx == 'Spata2l'] = 'Spata2L'
    fpkm.index = fpkm_idx


    # housekeeping?
    #
    # genes = [
    #     'Slc1a3',
    #     'Gja1',
    #     'Aldh1l1',
    #     'S100b',
    #     'Aqp4',
    #     'Gfap',
    #     'Cd44',
    #     'Aldoc',
    #     'Gfap',
    #     'Gapdh',
    #     'B2m'
    # ]
    # for g in genes:
    #     fig = plt.figure(figsize=(4.5, 7))
    #     ax = fig.add_subplot(111)
    #     ax.set_position([0.5, 0.08, 0.48, 0.9])
    #     data_by_symbol.loc[g].plot.barh(width=0.8, ax=ax)
    #     ax.set_xlabel(g)
    #     fig.savefig(os.path.join(OUTDIR, '%s.png' % g), dpi=200)
    #     fig.savefig(os.path.join(OUTDIR, '%s.pdf' % g))
        # plt.close(fig)

    # TODO: if required, make a single chart with all early-stage markers side by side
