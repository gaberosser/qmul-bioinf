import os
import re

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import ticker
import seaborn as sns
from scipy.cluster import hierarchy

from load_data import rnaseq_data
from microarray import process
from utils.output import unique_output_dir
from plotting import clustering, bar, heatmap

from scripts.rnaseq import gtf_reader

import references


N_GENES = 500


def plot_clustermap(data, yugene=False, n_genes=N_GENES, **kwargs):
    if yugene:
        data = process.yugene_transform(data)

    kwargs.setdefault('cmap', 'RdBu_r')

    mad = process.median_absolute_deviation(data, axis=1).sort_values(ascending=False)
    top_mad = mad.iloc[:n_genes].index
    z = hierarchy.linkage(data.loc[top_mad].transpose(), method='average', metric='correlation')
    cg = clustering.plot_clustermap(
        data.loc[top_mad],
        col_linkage=z,
        **kwargs
    )
    plt.setp(
        cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90
    )
    cg.gs.update(bottom=0.2)
    return cg


def plot_correlation_heatmap(data, yugene=False, n_genes=None, **kwargs):
    if yugene:
        data = process.yugene_transform(data)

    kwargs.setdefault('cmap', 'Reds')
    kwargs.setdefault('vmin', 0.)
    kwargs.setdefault('vmax', 1.)

    if n_genes is not None:
        mad = process.median_absolute_deviation(data, axis=1).sort_values(ascending=False)
        top_mad = mad.iloc[:n_genes].index
        data = data.loc[top_mad]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax = sns.heatmap(data.corr(), ax=ax, **kwargs)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.tight_layout()

    return ax


if __name__ == "__main__":

    SHOW_GENE_LABELS = False
    OUTDIR = unique_output_dir("astrocytes", reuse_empty=True)

    # GSE73721 (reference astrocytes, oligos, ...)
    obj73721 = rnaseq_data.gse73721(source='star', annotate_by='Ensembl Gene ID')

    # remove unneeded samples
    to_keep73721 = (
        ~obj73721.data.columns.str.contains('whole cortex')
        & ~obj73721.data.columns.str.contains('endo')
        & ~obj73721.data.columns.str.contains('tumor')
        & ~obj73721.data.columns.str.contains('myeloid')
    )


    # GSE61794 (H9-derived NSC x 2)
    obj61794 = rnaseq_data.gse61794(source='star', annotate_by='Ensembl Gene ID')

    # GBM paired samples
    objwtchg = rnaseq_data.gbm_astrocyte_nsc_samples_loader(source='star', annotate_by='Ensembl Gene ID')

    # Pollard (NSC x 2)
    objpollard = rnaseq_data.pollard_nsc(source='star', annotate_by='Ensembl Gene ID')

    # rRNA gene IDs
    rrna_ensg = set(gtf_reader.get_rrna())

    # MT gene_ids
    mt_ensg = set(gtf_reader.get_mitochondrial())

    # combine the data
    data = pd.concat((obj73721.data.loc[:, to_keep73721], obj61794.data, objwtchg.data, objpollard.data), axis=1)
    data = data.loc[data.index.str.contains('ENSG')]

    # compare with qPCR: comparing markers in the NSC samples and paired astrocytes
    astro_markers1 = [
        'S100B',
        'CD44',
        'ALDH1L1',
        'NFIA',
    ]
    data_markers = data.loc[
        references.gene_symbol_to_ensembl(astro_markers1),
        data.columns.str.contains('DURA')
    ]
    data_markers.index = astro_markers1
    series = [data_markers.iloc[i] for i in range(data_markers.shape[0])]
    colours = [
        '#0000FF',
        '#FF0000',
        '#00C000',
        '#AD07E3',
    ]

    # plot absolute counts

    fig = plt.figure(figsize=(8, 7.5))
    ax = fig.add_subplot(111)
    bar.grouped_bar_chart(series, ax=ax, colours=colours)
    ax.legend(data_markers.index)
    ax.set_position([0.125, .3, .85, .65])
    ax.set_ylabel('Raw counts')
    fig.savefig(os.path.join(OUTDIR, 'astro_markers_qpcr_abs_counts.pdf'))
    fig.savefig(os.path.join(OUTDIR, 'astro_markers_qpcr_abs_counts.png'), dpi=200)

    # plot relative fold change

    relfc18 = data_markers.loc[:, 'DURA018N2_ASTRO_DAY12'] / data_markers.loc[:, 'DURA018N2_NSC']
    relfc19 = data_markers.loc[:, 'DURA019N8C_ASTRO_DAY12'] / data_markers.loc[:, 'DURA019N8C_NSC']
    series = [
        pd.Series([1., relfc18.loc[g], relfc19.loc[g]], index=['NSC', 'ASTRO018', 'ASTRO019']) for g in relfc18.index
    ]
    fig = plt.figure(figsize=(8, 7.5))
    ax = fig.add_subplot(111)
    bar.grouped_bar_chart(series, ax=ax, colours=colours)
    ax.legend(data_markers.index, loc='upper left')
    ax.set_position([0.125, .15, .85, .8])
    ax.set_ylabel('Fold change')
    fig.savefig(os.path.join(OUTDIR, 'astro_markers_qpcr_fold_change.pdf'))
    fig.savefig(os.path.join(OUTDIR, 'astro_markers_qpcr_fold_change.png'), dpi=200)

    col_colors = pd.DataFrame(index=data.columns, columns=['group'])

    sample_groups = {
        'Hippocampal astrocytes': {'regex': 'Hippocampus', 'colour' : '#7fc97f'},
        'Fetal cortical astrocytes': {'regex': 'Fetal', 'colour' : '#beaed4'},
        'Cortical astrocytes': {'regex': re.compile(r'[0-9]*yo *ctx [0-9]* *astro', flags=re.IGNORECASE), 'colour' : '#fdc086'},
        'Cortical neurons': {'regex': 'ctx neuron', 'colour': '#ffff99'},
        'Oligodendrocytes': {'regex': 'oligo', 'colour': '#386cb0'},
        'Our iNSC': {'regex': '_NSC', 'colour': '#777777'},
        'H9 iNSC': {'regex': 'H9', 'colour': '#cccccc'},
        'Our induced astrocytes': {'regex': 'ASTRO', 'colour': 'black'},
    }

    cmap = {
        'Hippocampus': '#7fc97f',
        'Fetal': '#beaed4',
        re.compile(r'ctx [0-9]* *astro', flags=re.IGNORECASE): '#fdc086',
        'ctx neuron': '#ffff99',
        'oligo': '#386cb0',
        '_NSC': '#777777',
        'H9': '#cccccc',
        'ASTRO': 'black'
    }
    for t, col in cmap.items():
        col_colors.loc[col_colors.index.str.contains(t)] = col

    row_colors = pd.DataFrame(index=data.index, columns=['RNA type'])
    row_colors.loc[row_colors.index.isin(rrna_ensg)] = 'black'

    # plot clustermap with unmodified data (MT-rRNA still included)
    cg = plot_clustermap(
        data,
        z_score=0,
        col_colors=col_colors,
        show_gene_labels=SHOW_GENE_LABELS
    )
    cg.savefig(os.path.join(OUTDIR, 'clustermap_unmodified.png'), dpi=200)

    cg = plot_clustermap(
        data,
        yugene=True,
        col_colors=col_colors,
        show_gene_labels=SHOW_GENE_LABELS,
    )
    cg.savefig(os.path.join(OUTDIR, 'clustermap_unmodified_yg.png'), dpi=200)

    # plot rRNA quantity

    aa = data.loc[rrna_ensg].sum() / data.sum()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = aa.plot.bar(width=0.9)
    ax.set_ylim([0, 1])
    ax.set_ylabel("Proportion rRNA")
    plt.tight_layout()
    fig.savefig(os.path.join(OUTDIR, 'proportion_rrna.pdf'))
    fig.savefig(os.path.join(OUTDIR, 'proportion_rrna.png'), dpi=200)

    # remove rRNA

    data_rr = data.loc[~data.index.isin(rrna_ensg)]

    # plot MT RNA quantity (excluding rRNA)

    mt_no_rrna = mt_ensg.difference(rrna_ensg)
    bb = data_rr.loc[mt_no_rrna].sum() / data_rr.sum()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = bb.plot.bar(width=0.9)
    ax.set_ylim([0, .2])  # the maximum value observed is ~18%
    ax.set_ylabel("Proportion MT genes (excluding rRNA)")
    plt.tight_layout()
    fig.savefig(os.path.join(OUTDIR, 'proportion_mt.pdf'))
    fig.savefig(os.path.join(OUTDIR, 'proportion_mt.png'), dpi=200)

    # remove MT in addition to rRNA

    data_rr_mt = data_rr.loc[~data_rr.index.isin(mt_ensg)]

    # plot clustermap with rRNA removed
    cg = plot_clustermap(
        data_rr,
        z_score=0,
        col_colors=col_colors,
        show_gene_labels=SHOW_GENE_LABELS
    )
    cg.savefig(os.path.join(OUTDIR, 'clustermap_sub_rrna.png'), dpi=200)

    # correlation heatmap rRNA removed
    ax = plot_correlation_heatmap(data_rr.iloc[:, cg.dendrogram_col.reordered_ind])
    ax.figure.savefig(os.path.join(OUTDIR, 'correlation_sub_rrna.png'), dpi=200)
    ax.figure.savefig(os.path.join(OUTDIR, 'correlation_sub_rrna.pdf'))


    cg = plot_clustermap(
        data_rr,
        yugene=True,
        col_colors=col_colors,
        show_gene_labels=SHOW_GENE_LABELS
    )
    cg.savefig(os.path.join(OUTDIR, 'clustermap_sub_rrna_yg.png'), dpi=200)

    # plot clustermap with rRNA and MT RNA removed

    cg = plot_clustermap(
        data_rr_mt,
        z_score=0,
        col_colors=col_colors,
        show_gene_labels=SHOW_GENE_LABELS
    )
    cg.savefig(os.path.join(OUTDIR, 'clustermap_sub_rrna_mt.png'), dpi=200)

    # correlation heatmap rRNA and MT removed
    ax = plot_correlation_heatmap(data_rr_mt.iloc[:, cg.dendrogram_col.reordered_ind])
    ax.figure.savefig(os.path.join(OUTDIR, 'correlation_sub_rrna_mt.png'), dpi=200)
    ax.figure.savefig(os.path.join(OUTDIR, 'correlation_sub_rrna_mt.pdf'))

    cg = plot_clustermap(
        data_rr_mt,
        yugene=True,
        col_colors=col_colors,
        show_gene_labels=SHOW_GENE_LABELS
    )
    cg.savefig(os.path.join(OUTDIR, 'clustermap_sub_rrna_mt_yg.png'), dpi=200)

    # re-run clustering without Barres data (since this has different distribution)
    data_no_barres = pd.concat((obj61794.data, objwtchg.data, objpollard.data), axis=1)
    data_no_barres = data_no_barres.loc[data_no_barres.index.str.contains('ENSG')]

    cg = plot_clustermap(
        data_no_barres,
        z_score=0,
        col_colors=col_colors,
        show_gene_labels=SHOW_GENE_LABELS
    )
    cg.savefig(os.path.join(OUTDIR, 'clustermap_no_barres.png'), dpi=200)

    cg = plot_clustermap(
        data_no_barres,
        yugene=True,
        col_colors=col_colors,
        show_gene_labels=SHOW_GENE_LABELS
    )
    cg.savefig(os.path.join(OUTDIR, 'clustermap_no_barres_yg.png'), dpi=200)

    # for every sample, extract the top N by count and summarise

    topNs = [10, 50, 100]

    for topN in topNs:

        common_genes = set()
        top_dat = []
        for i in range(data_rr.shape[1]):
            t = data_rr_mt.iloc[:, i].sort_values(ascending=False)[:topN]
            common_genes.update(t.index)

        top_dat = data_rr_mt.loc[list(common_genes)].divide(data_rr.sum(), axis=1)
        top_dat.index = references.ensembl_to_gene_symbol(top_dat.index)

        cg = plot_clustermap(
            top_dat,
            n_genes=topN,  # plot all genes
            col_colors=col_colors,
            show_gene_labels=True,
            z_score=0
        )
        cg.savefig(os.path.join(OUTDIR, 'clustermap_sub_rrna_mt_top%dgenes.png' % topN), dpi=200)

    # bar charts of successive markers, used to characterise based on timeline
    # for this, only astrocytes and NSCs useful, so remove oligo and neuron
    astro_markers2 = [
        'NFIA',
        'SLC1A3',
        'ALDH1L1',
        'GJA1',
        'S100B',
        'CD44',
        'ALDOC',
        'GFAP',
        'AQP4',
        'SLC1A2'
    ]
    data_timeline = data.loc[
        references.gene_symbol_to_ensembl(astro_markers2),
        ~data.columns.str.contains('neuron')
        # & ~data.columns.str.contains('oligo')
    ]

    # normalise by dividing by the sum of non-rRNA, non-MT genes
    n = data_rr_mt.loc[
        :,
        ~data.columns.str.contains('neuron')
        # & ~data.columns.str.contains('oligo')
    ].sum(axis=0)
    data_timeline = data_timeline.divide(n, axis=1) * 1e6  # arbitrary multiplier to improve small number visualisation

    fig = plt.figure(figsize=(10.5, 4.8))
    gs = plt.GridSpec(1, len(astro_markers2), left=0.2, bottom=0.02, right=0.98, top=0.95, wspace=0.01)
    for i, m in enumerate(astro_markers2):
        ax = fig.add_subplot(gs[i])
        ax.barh(np.arange(data_timeline.shape[1]), data_timeline.iloc[i], height=0.9, align='center')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1.))
        if i == 0:
            ax.set_yticks(np.arange(data_timeline.shape[1]))
            ax.set_yticklabels(data_timeline.columns)
        else:
            ax.set_yticklabels([])
        ax.set_title(m)
        ax.set_xticks([])
        ax.set_ylim([-1, data_timeline.shape[1]])

    fig.savefig(os.path.join(OUTDIR, 'astro_markers.pdf'))
    fig.savefig(os.path.join(OUTDIR, 'astro_markers.png'), dpi=200)

    # more general neuronal lineage markers
    neuronal_lineage_markers = {
        'NSC': [
            'VIM',  # VIM
            'BMI1',
            'NES',
            'NEUROD1',
            'SOX2',
            'FABP7'
        ],
        'oligodendrocyte': [
            'GALC',
            'SOX10',
            'MOG',
        ],
        'OPC': [
            'PDGFRA',
            'NKX2-2',
            'OLIG2',
        ],
        'astrocyte': [
            'FGFR3',
            'GFAP',
            'S100B',
            'MSI1',
        ]
    }
    all_neuronal_markers = []
    for grp, arr in neuronal_lineage_markers.items():
        all_neuronal_markers.extend(arr)
    all_neuronal_markers_ens = references.gene_symbol_to_ensembl(all_neuronal_markers)

    n = data_rr_mt.loc[
        :,
        ~data.columns.str.contains('neuron')
        # & ~data.columns.str.contains('oligo')
    ].sum(axis=0)

    data_rr_mt_markers = data_rr_mt.loc[
        all_neuronal_markers_ens,
        ~data.columns.str.contains('neuron')
     ] / n
    data_rr_mt_markers.index = all_neuronal_markers
    data_rr_mt_markers = data_rr_mt_markers.divide(
        data_rr_mt_markers.sum(axis=1), axis=0
    )

    fig, axs, _, gs = heatmap.grouped_expression_heatmap(
        neuronal_lineage_markers.items(),
        data_rr_mt_markers,
        heatmap_kwargs={'cmap': 'Reds', 'square': False},
        cbar=False,
        vmin=0., vmax=0.1,
        fig_kwargs={'figsize': [5.5, 6.3]})
    gs.update(left=0.4, wspace=.1, bottom=0.2, top=0.98, right=.98)
    fig.savefig(os.path.join(OUTDIR, 'neuronal_lineage_marker_norm_expr.png'), dpi=200)
    fig.savefig(os.path.join(OUTDIR, 'neuronal_lineage_marker_norm_expr.pdf'))

    # playing around with dynamic range, etc.

    def dynamic_range_plot(data, ax):
        for lbl, d in sample_groups.items():
            this_data = data.loc[:, data.columns.str.contains(d['regex'])]
            first = True
            for t in this_data.columns:
                col = this_data.loc[:, t].sort_values()
                col /= col.max()
                x = np.linspace(0, 1, col.size)
                plt_lbl = None
                if first:
                    plt_lbl = lbl
                    first = False
                ax.plot(x, col.values, color=d['colour'], label=plt_lbl)


    data_rr_yg = process.yugene_transform(data_rr)
    data_rr_log = np.log(data_rr + 1)
    data_rr_log_yg = process.yugene_transform(data_rr_log)

    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)

    dynamic_range_plot(data_rr, axs[0, 0])
    dynamic_range_plot(data_rr_log, axs[0, 1])
    dynamic_range_plot(data_rr_yg, axs[1, 0])
    dynamic_range_plot(data_rr_log_yg, axs[1, 1])
    axs[0, 0].set_xlim(0.4, 1.)
    axs[0, 0].set_ylabel('Expression level (normalised)')
    axs[1, 0].set_ylabel('Expression level (normalised)')
    axs[1, 0].set_xlabel('Percentile')
    axs[1, 1].set_xlabel('Percentile')

    axs[0, 0].legend(loc='upper left', frameon=True, facecolor='w', framealpha=0.9)
    axs[0, 0].set_title('Raw counts')
    axs[0, 1].set_title('Log counts')
    axs[1, 0].set_title('YuGene counts')
    axs[1, 1].set_title('YuGene log counts')
    plt.tight_layout()

    fig.savefig(os.path.join(OUTDIR, "dynamic_range_transformations.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "dynamic_range_transformations.pdf"))

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