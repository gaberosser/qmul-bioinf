import csv
import multiprocessing as mp
import os

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy as hc

from plotting import clustering
from rnaseq import gsea
from scripts.hgic_final import analyse_xcell_results
from scripts.hgic_final import consts
from settings import HGIC_LOCAL_DIR, GIT_LFS_DATA_DIR
from utils import log, output, reference_genomes

logger = log.get_console_logger()


XCELL_SIGNATURE_FN = os.path.join(GIT_LFS_DATA_DIR, 'xcell', 'ESM3_signatures.xlsx')



def load_ipa_signatures(fn):
    res = {}
    with open(fn, 'rb') as f:
        c = csv.reader(f)
        for row in c:
            res[row[0]] = row[1:]
    return res


def simplify_tcga_names(data):
    """
    Simplify the TCGA naming to enable comparison with the GlioVis data.
    e.g. TCGA-76-4931-01A-01 becomes TCGA.76.4931
    In the process, duplicates can arise. We just pick the first one in these cases.
    :param data: pd.DataFrame containing the names in the columns
    :return: pd.DataFrame
    """
    out = data.copy()
    cols = out.columns.str.replace('-', '.')
    cols = cols.str.replace(r'\.[0-9A-Z]{3}\.[0-9]{2}$', '')
    out.columns = cols

    # following this renaming, we have duplicates in the columns
    # I've checked these and they appear to have been sequenced twice. Pick one arbitrarily
    dupes = out.columns[out.columns.duplicated()]
    if len(dupes) > 0:
        logger.warn(
            "After relabelling, there are %d duplicate sample names. We'll keep the first instance in each case.",
            len(dupes)
        )
        out = out.loc[:, ~out.columns.duplicated()]
    return out


def run_ssgsea(data, pathways, **ssgsea_kwds):
    pool = mp.Pool()
    res = {}
    jobs = {}
    for sn in data.columns:
        for pw, genes in pathways.items():
            jobs[(sn, pw)] = pool.apply_async(
                gsea.run_one_ssgsea,
                args=(data[sn], genes),
                kwds=ssgsea_kwds
            )
    pool.close()
    pool.join()
    for (sn, pw), j in jobs.items():
        res.setdefault(sn, {})[pw] = j.get()

    return res


if __name__ == '__main__':

    pids = consts.PIDS

    # reference comparators
    comparisons = ['h9', 'gibco']

    # cutoff for significant correlation
    corr_alpha = 0.05

    # cutoff used to determine whether a cell type has been detected
    xcell_alpha = 0.05

    # set this True to erase non-significant xCell results
    zero_nonsign_xcell = False

    # minimum number of samples with that cell type detectable
    min_n_samples = 4

    # cutoff for pct of genes shared between pathway and cell type signature
    pct_shared_max = 10.

    # correlation metric
    corr_metric = 'spearman'

    xcell_tcga_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/tcga/xcell/xcell_results_rnaseq_fpkm.xlsx'
    )
    outdir = output.unique_output_dir()

    xcell_res = pd.read_excel(xcell_tcga_fn, sheet_name=None)
    xcell_tcga = xcell_res['Proportions'].copy()
    xcell_tcga_pval = xcell_res['P values'].copy()

    # convert sample name format and resolve duplicates
    xcell_tcga = simplify_tcga_names(xcell_tcga)
    xcell_tcga_pval = simplify_tcga_names(xcell_tcga_pval)

    # drop the 'scores' rows
    xcell_tcga = xcell_tcga.loc[xcell_tcga_pval.index]

    # erase non-significant results if requested
    if zero_nonsign_xcell:
        xcell_tcga.loc[xcell_tcga_pval > xcell_alpha] = 0.

    # number of samples with detectable levels of cell types
    detected = (xcell_tcga_pval <= xcell_alpha)

    # filter based on this
    keep_ix = detected.sum(axis=1) >= min_n_samples
    logger.info(
        "The following %d cell types do not meet the requirement for detection in %d or more samples: \n%s",
        (~keep_ix).sum(),
        min_n_samples,
        ', '.join(keep_ix.index[keep_ix])
    )

    xcell_tcga = xcell_tcga.loc[keep_ix]
    xcell_tcga_pval = xcell_tcga_pval.loc[keep_ix]
    detected = detected.loc[keep_ix]

    logger.info("%d cell types remain." % keep_ix.sum())

    # rnaseq GBM metadata (GlioVis)
    tcga_gbm_meta_fn = os.path.join(HGIC_LOCAL_DIR, 'current', 'input_data', 'tcga', 'GlioVis_TCGA_GBMLGG.meta.xlsx')
    tcga_gbm_meta_full = pd.read_excel(tcga_gbm_meta_fn, header=0, index_col=0)
    # filter to include only GBM IDH1wt
    tcga_gbm_meta = tcga_gbm_meta_full.loc[
        (tcga_gbm_meta_full.Histology == 'GBM') & (tcga_gbm_meta_full['IDH.status'] == 'WT')
    ]

    # metadata (Wang classification)
    wang_fn = os.path.join(HGIC_LOCAL_DIR, 'current', 'input_data', 'tcga', 'wang_table_s4_tcga_classification_idh1_status.xlsx')
    wang_meta = pd.read_excel(wang_fn, skiprows=2, sheet_name='GBM', header=0, index_col=0)
    wang_meta.index = wang_meta.index.str.replace(r'\.[0-9]{2}$', '')

    # load TCGA data in FPKM units
    tcga_rnaseq_data_fn = os.path.join(HGIC_LOCAL_DIR, 'current', 'input_data', 'tcga', 'rnaseq.xlsx')
    tcga_rnaseq_data = pd.read_excel(tcga_rnaseq_data_fn, header=0, index_col=0, sheet_name='fpkm')

    # convert sample name format and resolve duplicates
    tcga_rnaseq_data = simplify_tcga_names(tcga_rnaseq_data)

    # reduce to intersection of samples
    # we require the sample to be present in three sources:
    # - xCell cell type composition matrix
    # - TCGA data (RPKM)
    # - GlioVis metadata
    ix = tcga_gbm_meta.index.intersection(tcga_rnaseq_data.columns).intersection(xcell_tcga.columns)
    logger.info(
        "Of the %d samples in the GlioVis metadata, we retain %d that are also in the xCell and TCGA matrices.",
        tcga_gbm_meta.shape[0],
        ix.size
    )
    tcga_gbm_meta = tcga_gbm_meta.loc[ix]
    tcga_rnaseq_data = tcga_rnaseq_data.loc[:, ix]
    xcell_tcga = xcell_tcga.loc[:, ix]
    xcell_tcga_pval = xcell_tcga_pval.loc[:, ix]

    fig = plt.figure(figsize=(6, 9))
    ax = fig.add_subplot(111)
    detected.sum(axis=1).sort_values(ascending=False).plot.barh(ax=ax, color='k')
    ax.set_xlabel("Number of samples with detectable levels (%d total)" % xcell_tcga.shape[1])
    ax.set_xlim([0, xcell_tcga.shape[1]])
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "number_detectable_by_cell_type.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "number_detectable_by_cell_type.tiff"), dpi=200)
    fig.savefig(os.path.join(outdir, "number_detectable_by_cell_type.pdf"), dpi=200)

    # box and whisker
    palette = sns.color_palette("deep", 16)
    fig = plt.figure(figsize=(11, 5.5))
    ax = fig.add_subplot(111)
    sns.boxplot(
        data=xcell_tcga.transpose(),
        fliersize=0,
        ax=ax,
        zorder=2
    )
    # remove face fill
    plt.setp(ax.artists, edgecolor='k', facecolor='none', alpha=1., zorder=2)
    sns.stripplot(
        data=xcell_tcga.transpose(),
        palette=palette,
        jitter=0.3,
        ax=ax,
        edgecolor='none',
        linewidth=0.5,
        # alpha=0.5,
        size=2,
        zorder=1
    )
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    fig.subplots_adjust(bottom=0.5, left=0.05, right=0.99, top=0.98)
    fig.savefig(os.path.join(outdir, "box_and_swarmplot_all_types.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "box_and_swarmplot_all_types.tiff"), dpi=200)

    # stdev across samples
    stdev = xcell_tcga.std(axis=1).sort_values(ascending=False)

    fig = plt.figure(figsize=(6, 9))
    ax = fig.add_subplot(111)
    stdev.plot.barh(color='k', ax=ax)
    fig.subplots_adjust(bottom=0.07, left=0.4, right=0.98, top=0.98)
    ax.set_xlabel('Stdev. across samples')
    fig.savefig(os.path.join(outdir, "stdev_across_samples.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "stdev_across_samples.tiff"), dpi=200)

    # CV across samples (may be more relevant?)
    cv = (xcell_tcga.std(axis=1) / xcell_tcga.mean(axis=1)).sort_values(ascending=False)

    fig = plt.figure(figsize=(6, 9))
    ax = fig.add_subplot(111)
    cv.plot.barh(color='k', ax=ax)
    fig.subplots_adjust(bottom=0.07, left=0.4, right=0.98, top=0.98)
    ax.set_xlabel('CV across samples')
    fig.savefig(os.path.join(outdir, "cv_across_samples.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "cv_across_samples.tiff"), dpi=200)

    # heatmap: proportions for each patient
    # standardise across columns, because each cell type has different mean proportion
    rl = hc.linkage(xcell_tcga.astype(float).transpose(), method='average', metric='euclidean')

    # row colours: identify transcription and methylation classification
    wang_class = pd.Series('Unknown', index=xcell_tcga.columns)
    ix = wang_meta.index.intersection(xcell_tcga.columns)
    wang_class.loc[ix] = wang_meta.loc[ix, 'Group']

    sturm_class = pd.Series('Unknown', index=xcell_tcga.columns)
    ix = tcga_gbm_meta['Random.Forest.Sturm.Cluster'].dropna().index
    sturm_class.loc[ix] = tcga_gbm_meta.loc[ix, 'Random.Forest.Sturm.Cluster']

    row_colours = pd.DataFrame(index=xcell_tcga.columns)
    sturm_cmap = {
        "RTK I 'PDGFRA'": consts.SUBGROUP_SET_COLOURS['RTK I full'],
        "RTK II 'Classic'": consts.SUBGROUP_SET_COLOURS['RTK II full'],
        "Mesenchymal": consts.SUBGROUP_SET_COLOURS['MES full'],
        "Unknown": 'gray'
    }
    wang_cmap = {
        'PN': consts.SUBGROUP_SET_COLOURS['RTK I partial'],
        'CL': consts.SUBGROUP_SET_COLOURS['RTK II partial'],
        'MS': consts.SUBGROUP_SET_COLOURS['MES partial'],
        'Unknown': 'grey'
    }
    row_colours.insert(0, 'Sturm', sturm_class.map(sturm_cmap))
    row_colours.insert(0, 'Verhaak', wang_class.map(wang_cmap))

    cg = clustering.plot_clustermap(
        xcell_tcga.astype(float).transpose(),
        metric='euclidean',
        show_gene_labels=False,
        show_gene_clustering=True,
        cmap='YlOrRd',
        row_linkage=rl,
        z_score=1,
        vmin=-1.5,
        vmax=6.,
        row_colors=row_colours
    )
    cg.gs.update(left=0.03, bottom=0.22, right=0.9)
    c_labels = [''] * len(cg.cax.get_yticks())
    c_labels[0] = 'Low'
    c_labels[-1] = 'High'
    cg.cax.set_yticklabels(c_labels)
    cg.cax.set_ylabel('Normalised proportion', labelpad=-70)  # bit hacky, but this places the label correctly
    cg.savefig(os.path.join(outdir, "cell_proportion_cluster_by_patient.png"), dpi=200)
    cg.savefig(os.path.join(outdir, "cell_proportion_cluster_by_patient.tiff"), dpi=200)
    cg.savefig(os.path.join(outdir, "cell_proportion_cluster_by_patient.pdf"), dpi=200)

    ## TODO: add legend to plot, explaining row_colors

    # this demonstrates that TCGA.06.5410 is an outlier
    # remove it? (TODO?)

    # load IPA signatures
    ipa_sign_fn = os.path.join(HGIC_LOCAL_DIR, 'current', 'input_data', 'ipa_pathways', 'ipa_exported_pathways.csv')
    ipa_signatures = load_ipa_signatures(ipa_sign_fn)

    # run ssGSEA separately for each sample
    ssgsea_rnaseq_data = run_ssgsea(tcga_rnaseq_data, ipa_signatures)
    ssgsea_rnaseq_data = pd.DataFrame(ssgsea_rnaseq_data)

    # carry out correlation analysis
    co, co_p = analyse_xcell_results.pathway_cell_type_composition_correlation_analysis(
        ssgsea_rnaseq_data,
        xcell_tcga,
        corr_metric=corr_metric
    )

    # heatmap showing correlation between pathways and cell types

    # precursor: check for cases where there is a substantial overlap in genes in pathways and cell type signatures
    # load xCell signatures
    xcell_s = pd.read_excel(XCELL_SIGNATURE_FN, header=0, index_row=0)
    xcell_signatures = {}
    for i, row in xcell_s.iterrows():
        xcell_signatures[row.Celltype_Source_ID] = set(row.iloc[2:].dropna().values)

    # convert IPA pathway Ensembl IDs to symbols for compatibility
    ipa_signatures_symb = {}
    for k, v in ipa_signatures.items():
        ipa_signatures_symb[k] = reference_genomes.ensembl_to_gene_symbol(v).dropna()

    # compute overlap between cell type signatures and IPA signatures
    pct_shared = analyse_xcell_results.compute_cell_type_pathway_overlap(
        ipa_signatures_symb,
        xcell_signatures,
    )

    # aggregate taking max over pathways
    cc = pct_shared.columns.str.replace(r'(?P<ct>[^_]*)_.*', r'\g<ct>')
    pct_shared_aggr = pct_shared.groupby(cc, axis=1).max()

    # set of pathways with any significance
    logger.info("%d pathways enriched in at least one patient and retained after correlation analysis" % co.shape[1])

    # run clustering to order the rows/cols nicely
    rl = hc.linkage(co.fillna(0.).transpose(), method='average', metric='euclidean')
    row_ix = hc.leaves_list(rl)
    cl = hc.linkage(co.fillna(0.), method='average', metric='euclidean')
    col_ix = hc.leaves_list(cl)

    # reorder the data based on the clustering
    co = co.iloc[col_ix, row_ix]
    co_p = co_p.iloc[col_ix, row_ix]

    # for plotting, we only need an indicator of which values are significant
    plot_dict = analyse_xcell_results.plot_heatmap_with_quantification(
        co,
        co_p,
        alpha=corr_alpha,
        hatch_df=pct_shared_aggr.loc[co.columns, co.index] > pct_shared_max,
        figsize=(9., 7.6)
    )
    plt.setp(plot_dict['main_ax'].yaxis.get_ticklabels(), fontsize=8)
    gs = plot_dict['gs']
    gs.update(left=0.45, bottom=0.25, top=0.99, right=0.93, wspace=0.03)
    fig = plot_dict['fig']
    fig.savefig(os.path.join(outdir, "cell_proportion_pathway_ssgsea_%s_clustering_tcga.png" % corr_metric), dpi=200)
    fig.savefig(os.path.join(outdir, "cell_proportion_pathway_ssgsea_%s_clustering_tcga.tiff" % corr_metric), dpi=200)
    fig.savefig(os.path.join(outdir, "cell_proportion_pathway_ssgsea_%s_clustering_tcga.pdf" % corr_metric), dpi=200)

    # show one scatterplot to exemplify this
    pw = 'Chondroitin Sulfate Biosynthesis'
    ct = 'Tregs'
    this_ssgsea = ssgsea_rnaseq_data.loc[pw]
    this_xcell = xcell_tcga.loc[ct, this_ssgsea.index]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(this_xcell, this_ssgsea)
    ax.set_xlabel('Proportion %s' % ct)
    ax.set_ylabel('ssGSEA score for %s' % pw)
    fig.tight_layout()
    fig.savefig(
        os.path.join(
            outdir,
            "example_scatter_%s_%s.png" % (ct.lower().replace(' ', '_'), pw.lower().replace(' ', '_'))
        ),
        dpi=200
    )


