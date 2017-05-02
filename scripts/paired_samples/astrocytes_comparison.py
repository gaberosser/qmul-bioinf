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


def plot_clustermap(data, yugene=False, n_genes=N_GENES, yugene_resolve_ties=False, **kwargs):
    if yugene:
        data = process.yugene_transform(data, resolve_ties=yugene_resolve_ties)

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

    # it is helpful to have access to the row index so we'll add it here
    # I *think* certain kwargs might cause this to fail (if no row dend has been computed?) so add a generic try-exc
    try:
        cg.row_index = top_mad[cg.dendrogram_row.reordered_ind]
    except Exception:
        pass

    return cg


def plot_all_clustermaps(data, filestem, **kwargs):
    """
    Run the clustermap process and save to disk for each of the possible data transformation combinations:
    Counts
    Log (counts + 1)
    YuGene counts
    YuGene log (counts + 1)
    :param data:
    :param filestem:
    :param yugene:
    :param n_genes:
    :param kwargs:
    :return:
    """
    # store column ordering
    col_order = {}

    # 1) counts
    cg = plot_clustermap(
        data,
        z_score=0,
        **kwargs
    )
    cg.savefig('%s.png' % filestem, dpi=200)
    col_order['counts'] = cg.dendrogram_col.reordered_ind

    # 2) yugene counts
    cg = plot_clustermap(
        data,
        yugene=True,
        **kwargs
    )
    cg.savefig('%s_yg.png' % filestem, dpi=200)
    col_order['yg_counts'] = cg.dendrogram_col.reordered_ind

    # 3) log (counts + 1)
    cg = plot_clustermap(
        np.log(data + 1),
        z_score=0,
        **kwargs
    )
    cg.savefig('%s_log.png' % filestem, dpi=200)
    col_order['log_counts'] = cg.dendrogram_col.reordered_ind

    # 4) yugene log (counts + 1)
    cg = plot_clustermap(
        np.log(data + 1),
        yugene=True,
        **kwargs
    )
    cg.savefig('%s_log_yg.png' % filestem, dpi=200)
    col_order['yg_log_counts'] = cg.dendrogram_col.reordered_ind

    return col_order


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

def plot_all_correlation_heatmaps(data, filestem, col_orders, **kwargs):
    ax = plot_correlation_heatmap(data.iloc[:, col_orders['counts']], **kwargs)
    ax.figure.savefig("%s.png" % filestem, dpi=200)
    ax.figure.savefig("%s.pdf" % filestem)

    t = process.yugene_transform(data)
    ax = plot_correlation_heatmap(t.iloc[:, col_orders['yg_counts']], **kwargs)
    ax.figure.savefig("%s_yg.png" % filestem, dpi=200)
    ax.figure.savefig("%s_yg.pdf" % filestem)

    t = np.log(data + 1)
    ax = plot_correlation_heatmap(t.iloc[:, col_orders['log_counts']], **kwargs)
    ax.figure.savefig("%s_log.png" % filestem, dpi=200)
    ax.figure.savefig("%s_log.pdf" % filestem)

    t = process.yugene_transform(np.log(data + 1))
    ax = plot_correlation_heatmap(t.iloc[:, col_orders['yg_log_counts']], **kwargs)
    ax.figure.savefig("%s_log_yg.png" % filestem, dpi=200)
    ax.figure.savefig("%s_log_yg.pdf" % filestem)


if __name__ == "__main__":

    SHOW_GENE_LABELS = False
    OUTDIR = unique_output_dir("astrocytes", reuse_empty=True)
    INCLUDE_ALL_NSC = True
    COMBINE_REPLICATES = True  # merges the H9 replicates
    SIDETRACK_1 = False

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
    if COMBINE_REPLICATES:
        rc = obj61794.meta.read_count.sum()
        obj61794.meta = pd.DataFrame(
            data={
                'cell_type': 'NSC',
                'srr': 'SRR1586371-2',
                'read_count': rc,
                'sample': 'H9 NSC',
            },
            index=['SRR1586371-2']
        )
        obj61794.data = pd.DataFrame(obj61794.data.sum(axis=1), columns=['H9 NSC'])


    # GBM paired samples
    objwtchg_paired = rnaseq_data.gbm_astrocyte_nsc_samples_loader(source='star', annotate_by='Ensembl Gene ID')

    # WTCHG ALL samples
    objwtchg_all = rnaseq_data.all_wtchg_loader(source='star', annotate_by='Ensembl Gene ID')
    to_keepwtchg = (
        objwtchg_all.data.columns.str.contains('NSC')
        | objwtchg_all.data.columns.str.contains('ASTRO')
    )
    objwtchg_all.data = objwtchg_all.data.loc[:, to_keepwtchg]
    objwtchg_all.meta = objwtchg_all.meta.loc[to_keepwtchg, :]

    # Pollard (NSC x 2)
    objpollard = rnaseq_data.pollard_nsc(source='star', annotate_by='Ensembl Gene ID')

    # rRNA gene IDs
    rrna_ensg = set(gtf_reader.get_rrna())

    # MT gene_ids
    mt_ensg = set(gtf_reader.get_mitochondrial())

    # combine the data
    if INCLUDE_ALL_NSC:
        data = pd.concat((obj73721.data.loc[:, to_keep73721], obj61794.data, objwtchg_all.data, objpollard.data), axis=1)

        # combine the metadata
        meta = pd.concat((
            obj73721.meta.loc[to_keep73721],
            obj61794.meta,
            objwtchg_all.meta,
            objpollard.meta
        ), axis=0)

    else:
        data = pd.concat((obj73721.data.loc[:, to_keep73721], obj61794.data, objwtchg_paired.data, objpollard.data), axis=1)

        # combine the metadata
        meta = pd.concat((
            obj73721.meta.loc[to_keep73721],
            obj61794.meta,
            objwtchg_paired.meta,
            objpollard.meta
        ), axis=0)

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
        'Pollard NSC': {'regex': 'Pollard', 'colour': '#cccccc'},
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
    for d in sample_groups.values():
        col_colors.loc[col_colors.index.str.contains(d['regex'])] = d['colour']

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

    filestem = os.path.join(OUTDIR, 'clustermap_sub_rrna_mt')
    col_order = plot_all_clustermaps(data_rr_mt, filestem, col_colors=col_colors)

    filestem = os.path.join(OUTDIR, 'correlation_sub_rrna_mt')
    plot_all_correlation_heatmaps(data_rr_mt, filestem, col_order, vmin=0.5, vmax=1.)

    # re-run clustering without Barres data (since this has different distribution)
    data_no_barres = pd.concat((obj61794.data, objwtchg_paired.data, objpollard.data), axis=1)
    data_no_barres = data_no_barres.loc[data_no_barres.index.str.contains('ENSG')]

    filestem = os.path.join(OUTDIR, 'clustermap_no_barres')
    col_order = plot_all_clustermaps(data_no_barres, filestem, col_colors=col_colors)

    filestem = os.path.join(OUTDIR, 'correlation_no_barres')
    plot_all_correlation_heatmaps(data_no_barres, filestem, col_order, vmin=0.5, vmax=1.)

    # Requested by LG: all data except our astrocytes
    data_no_iastro = data_rr_mt.loc[:, ~data_rr_mt.columns.str.contains('_ASTRO')]
    filestem = os.path.join(OUTDIR, 'clustermap_no_iastro')
    col_order = plot_all_clustermaps(data_no_iastro, filestem, col_colors=col_colors)

    # all NSCs
    data_nsc = data_rr_mt.loc[:, data_rr_mt.columns.str.contains('NSC')]
    filestem = os.path.join(OUTDIR, 'clustermap_nsc_only')
    col_order = plot_all_clustermaps(data_nsc, filestem, col_colors=col_colors)

    # for every sample, extract the top N by count and summarise

    topNs = [10, 50, 100]

    for topN in topNs:

        common_genes = set()
        top_dat = []
        for i in range(data_rr.shape[1]):
            t = data_rr_mt.iloc[:, i].sort_values(ascending=False)[:topN]
            common_genes.update(t.index)

        top_dat = data_rr_mt.loc[list(common_genes)].divide(data_rr.sum(), axis=1)
        symb = references.ensembl_to_gene_symbol(top_dat.index)
        tidx = np.array(top_dat.index)
        tidx[~symb.isnull().values] = symb.loc[~symb.isnull()].values
        top_dat.index = tidx

        filestem = os.path.join(OUTDIR, 'clustermap_sub_rrna_mt_top_%d' % topN)
        col_order = plot_all_clustermaps(top_dat, filestem, col_colors=col_colors)

        filestem = os.path.join(OUTDIR, 'correlation_sub_rrna_mt_top_%d' % topN)
        plot_all_correlation_heatmaps(top_dat, filestem, col_order, vmin=0.5, vmax=1.)


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

    # for each marker, plot the MAX normalised expression for comparison
    marker_level = data_rr.loc[all_neuronal_markers_ens] / data_rr_mt.sum(axis=0)
    marker_level.index = all_neuronal_markers

    fig = plt.figure(figsize=(6.8, 3.2))
    ax = fig.add_subplot(111)
    sns.boxplot(marker_level.transpose(), ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    # max_marker_level.plot.bar(ax=ax)
    ax.set_ylabel("Proportion of reads")
    plt.tight_layout()
    fig.savefig(os.path.join(OUTDIR, 'neuronal_lineage_marker_level_boxplot.png'), dpi=200)
    fig.savefig(os.path.join(OUTDIR, 'neuronal_lineage_marker_level_boxplot.pdf'))


    # playing around with dynamic range, etc.

    def dynamic_range_plot(data, ax=None, renorm=True):
        """
        Show the dynamic range of the supplied data.
        :param data:
        :param ax:
        :param renorm: If True, rescale the x axis so that zero represents the first non-zero y value.
        :return:
        """
        for lbl, d in sample_groups.items():
            this_data = data.loc[:, data.columns.str.contains(d['regex'])]
            first = True
            for t in this_data.columns:
                col = this_data.loc[:, t].sort_values()
                col /= col.max()
                if renorm:
                    col = col.loc[col > 0]
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
    axs[0, 0].set_xlim(0., 1.)
    axs[0, 0].set_ylabel('Expression level (normalised)')
    axs[1, 0].set_ylabel('Expression level (normalised)')
    axs[1, 0].set_xlabel('Non-zero percentile')
    axs[1, 1].set_xlabel('Non-zero percentile')

    axs[0, 0].legend(loc='upper left', frameon=True, facecolor='w', framealpha=0.9)
    axs[0, 0].set_title('Raw counts')
    axs[0, 1].set_title('Log counts')
    axs[1, 0].set_title('YuGene counts')
    axs[1, 1].set_title('YuGene log counts')
    plt.tight_layout()

    fig.savefig(os.path.join(OUTDIR, "dynamic_range_transformations.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "dynamic_range_transformations.pdf"))

    # little sidetrack: resample reads from a high count Barres datum.
    # This is achieved by weighted random sampling of the genes based on the empirical distribution function
    # Key Q: Does this support the low count YG distribution?
    # Spoiler alert: Yes, it seems so.

    # takes a while, so allow disabling it
    if SIDETRACK_1:

        total_read_counts = [200000, 500000, 1000000, 2000000, 5000000]
        n_rpt = 100

        # set the reference (~7.5mi reads)
        ref_data = data_rr.loc[:, '13yo ctx astro'].sort_values()
        ref_data /= ref_data.sum()
        ref_data = ref_data.loc[ref_data != 0]

        sim_data = {}
        for trc in total_read_counts:
            print "Read count %d" % trc
            simd = np.random.choice(ref_data.index, p=ref_data.values, size=(trc, n_rpt))
            # count for each column
            counted = [pd.value_counts(simd[:, i]) for i in range(n_rpt)]
            sim_data[trc] = pd.DataFrame(counted, columns=ref_data.index).transpose().fillna(0)

        # yugene transform the ref_data
        ref_data_yg = process.yugene_transform(pd.DataFrame(ref_data))

        # yugene transform each of the expt runs
        sim_data_yg = dict([
            (k, process.yugene_transform(v)) for k, v in sim_data.iteritems()
        ])

        # plot them all together using a modified plot code and look for a difference in the distribution
        fig = plt.figure()
        ax = fig.add_subplot(111)

        sim_colours = {
            200000: 'y',
            500000: 'b',
            1000000: 'c',
            2000000: 'r',
            5000000: 'g'
        }
        ref_col = ref_data_yg.iloc[:, 0].sort_values()
        ref_col = ref_col.loc[ref_col > 0]
        ref_x = x = np.linspace(0, 1, ref_col.size)
        ax.plot(ref_x, ref_col, color='k')

        for n, arr in sim_data_yg.items():
            first = True
            for t in arr.columns:
                col = arr.loc[:, t].sort_values()
                col = col.loc[col > 0]
                x = np.linspace(0, 1, col.size)
                if first:
                    lbl = str(n)
                    first = False
                else:
                    lbl = None
                ax.plot(x, col.values, color=sim_colours[n], label=lbl)

        ax.set_xlabel('Rescaled percentile')
        ax.set_ylabel('Normalised expression value')
        ax.set_ylim([0, 1])
        ax.legend(loc='upper left', frameon=True, facecolor='w', framealpha=0.9)
        fig.savefig(os.path.join(OUTDIR, "dynamic_range_simulations.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "dynamic_range_simulations.pdf"))


    # for reference, we can also load the original 'pre-processed' GSE73721 data
    # but these are indexed by mouse gene??
    # fpkm, meta = rnaseq_data.brainrnaseq_preprocessed()
    # fpkm_idx = np.array(fpkm.index.str.capitalize())
    # fpkm_idx[fpkm_idx == 'Spata2l'] = 'Spata2L'
    # fpkm.index = fpkm_idx


    # Seaborn pairplot of NSC (scatter matrix) with two populations of points:
    # In both cases, require that more than half of the samples have a read count > MIN_COUNT
    # UNLESS one sample has a reading > HIGH_COUNT
    # top 500 by MAD
    # randomly drawn 500 from non-MAD
    MIN_COUNT = 4
    MAX_BELOW = int(np.ceil(data_nsc.shape[1] * 0.5))
    HIGH_COUNT = 100

    nz_idx = (
        ((data_nsc < MIN_COUNT).sum(axis=1) < MAX_BELOW) | ((data_nsc > HIGH_COUNT).any(axis=1))
    )
    data_nsc_nz = data_nsc.loc[nz_idx, :]

    # yugene
    data_nsc_nz_yg = process.yugene_transform(data_nsc_nz)

    # add one, norm, take log
    data_nsc_nz += 1
    data_nsc_nz = data_nsc_nz.divide(data_nsc_nz.sum(axis=0), axis=1)
    data_nsc_nz = np.log(data_nsc_nz + 1)

    # MAD - compute on normalised values
    mad_nsc_nz = process.median_absolute_deviation(data_nsc_nz).sort_values(ascending=False)
    top_idx = mad_nsc_nz.index[:N_GENES]
    rem_idx = mad_nsc_nz.index[N_GENES:]

    # reduce number of remainder for plotting purposes
    to_discard = rem_idx[np.random.permutation(rem_idx.size)[N_GENES:]]
    data_nsc_nz = data_nsc_nz.drop(to_discard)

    # add 'hue' column
    data_nsc_nz.loc[:, 'hue'] = 'Remainder'
    data_nsc_nz.loc[top_idx, 'hue'] = 'Top %d by MAD' % N_GENES

    # generate the plot
    # pg = sns.pairplot(data_nsc_nz, hue='hue')

    # repeat on YuGene
    mad_nsc_nz_yg = process.median_absolute_deviation(data_nsc_nz_yg).sort_values(ascending=False)
    top_idx = mad_nsc_nz_yg.index[:N_GENES]
    rem_idx = mad_nsc_nz_yg.index[N_GENES:]

    # reduce number of remainder for plotting purposes
    to_discard = rem_idx[np.random.permutation(rem_idx.size)[N_GENES:]]
    data_nsc_nz_yg = data_nsc_nz_yg.drop(to_discard)

    data_nsc_nz_yg.loc[:, 'hue'] = 'Remainder'
    data_nsc_nz_yg.loc[top_idx, 'hue'] = 'Top %d by MAD' % N_GENES

    pg = sns.pairplot(data_nsc_nz_yg, hue='hue')
    delta = 0.04
    pg.set(ylim=[-delta, 1 + delta], xlim=[-delta, 1 + delta])
    pg._legend.set_visible(False)
    pg.fig.subplots_adjust(bottom=0.05)
    pg.savefig(os.path.join(OUTDIR, 'scatter_matrix_yg.png'), dpi=200)
    pg.savefig(os.path.join(OUTDIR, 'scatter_matrix_yg.pdf'))

