"""
Here we follow exactly the same route taken in the Sulman validation, but with our own GIC RNA-Seq samples.
The idea behind this is to work out how different the ssGSEA approach is (run on GIC only) to the IPA approach
(run on DE gene list from a GBM-iNSC comparison)
"""

import pandas as pd
import os
import re
import csv
from settings import HGIC_LOCAL_DIR, GIT_LFS_DATA_DIR
from matplotlib import pyplot as plt
import seaborn as sns
import multiprocessing as mp
from scipy.cluster import hierarchy as hc
from plotting import clustering
from utils import log, output
from rnaseq import gsea
from scripts.hgic_final import consts
from scripts.hgic_final import analyse_xcell_results
import references
logger = log.get_console_logger()


XCELL_SIGNATURE_FN = os.path.join(GIT_LFS_DATA_DIR, 'xcell', 'ESM3_signatures.xlsx')


def load_ipa_signatures(fn):
    res = {}
    with open(fn, 'rb') as f:
        c = csv.reader(f)
        for row in c:
            ix = row[0].decode('utf-8')
            res[ix] = row[1:]
    return res


if __name__ == "__main__":
    """
    Here, we repeat the IPA - xCell analysis previously run on our own data, this time applying it to a dataset
    of matched GBM and GIC lines 'published' by Verhaak. The dataset is not actually public, and is held by Erik Sulman.
    His PhD student Jie performed the xCell and ssGSEA analysis (ssGSEA is a proxy for IPA results here).

    Outline:
    1) Establish overlap between xCell signatures and the pathways being considered.
    2) Compute correlation between the two
    3) Plot
    4) Submit directly to Nature; acceptance with no revisions required.

    """
    indir = os.path.join(GIT_LFS_DATA_DIR, 'gic_ssgsea')

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

    xcell_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/characterisation/ffpe_cell_type_determination/xcell/xcell_results_salmon_tpm.xlsx'
    )

    ssgsea_norm_fn = os.path.join(
        indir,
        'norm.score.txt'
    )
    ssgsea_raw_fn = os.path.join(
        indir,
        'raw.score.txt'
    )
    ssgsea_pathways_fn = os.path.join(
        indir,
        'ipa_pathway_numbering.csv'
    )

    # relevant FFPE samples
    rna_ffpe_samples = [
        'NH15_1661DEF2C',
        'NH15_1877_SP1C',
        'NH15_2101_DEF1A',
        'NH16_270_DEF1Ereplacement',
        'NH16_616DEF1B',
        'NH16_677_SP1A',
        'NH16_2063_DEF1Areplacement',
        'NH16_2214DEF1A',
        'NH16_2255DEF1B2',
        'NH16_2806DEF3A1'
    ]

    outdir = output.unique_output_dir()

    xcell_res = pd.read_excel(xcell_fn, sheet_name=None)
    xcell_prop = xcell_res['Proportions'].copy()
    xcell_pval = xcell_res['P values'].copy()

    # drop the 'scores' rows
    xcell_prop = xcell_prop.loc[xcell_pval.index]

    # erase non-significant results if requested
    if zero_nonsign_xcell:
        xcell_prop.loc[xcell_pval > xcell_alpha] = 0.

    # reduce to samples in study
    xcell_prop = xcell_prop[rna_ffpe_samples]
    xcell_pval = xcell_pval[rna_ffpe_samples]

    # number of samples with detectable levels of cell types
    detected = (xcell_pval <= xcell_alpha)

    # filter based on this
    keep_ix = detected.sum(axis=1) >= min_n_samples
    if not keep_ix.all():
        logger.info(
            "The following %d cell types do not meet the requirement for detection in %d or more samples: \n%s",
            (~keep_ix).sum(),
            min_n_samples,
            ', '.join(keep_ix.index[~keep_ix])
        )

        xcell_prop = xcell_prop.loc[keep_ix]
        xcell_pval = xcell_pval.loc[keep_ix]
        detected = detected.loc[keep_ix]

        logger.info(
            "%d cell types remain: %s.",
            keep_ix.sum(),
            ', '.join(keep_ix.index[keep_ix])
        )

    fig = plt.figure(figsize=(6, 9))
    ax = fig.add_subplot(111)
    detected.sum(axis=1).sort_values(ascending=False).plot.barh(ax=ax, color='k')
    ax.set_xlabel("Number of samples with detectable levels (%d total)" % xcell_prop.shape[1])
    ax.set_xlim([0, xcell_prop.shape[1]])
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "number_detectable_by_cell_type.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "number_detectable_by_cell_type.tiff"), dpi=200)
    fig.savefig(os.path.join(outdir, "number_detectable_by_cell_type.pdf"), dpi=200)

    # box and whisker
    palette = sns.color_palette("deep", 16)
    fig = plt.figure(figsize=(11, 5.5))
    ax = fig.add_subplot(111)
    sns.boxplot(
        data=xcell_prop.transpose(),
        fliersize=0,
        ax=ax,
        zorder=2
    )
    # remove face fill
    plt.setp(ax.artists, edgecolor='k', facecolor='none', alpha=1., zorder=2)
    sns.stripplot(
        data=xcell_prop.transpose(),
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
    stdev = xcell_prop.std(axis=1).sort_values(ascending=False)

    fig = plt.figure(figsize=(6, 9))
    ax = fig.add_subplot(111)
    stdev.plot.barh(color='k', ax=ax)
    fig.subplots_adjust(bottom=0.07, left=0.4, right=0.98, top=0.98)
    ax.set_xlabel('Stdev. across samples')
    fig.savefig(os.path.join(outdir, "stdev_across_samples.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "stdev_across_samples.tiff"), dpi=200)

    # CV across samples (may be more relevant?)
    cv = (xcell_prop.std(axis=1) / xcell_prop.mean(axis=1)).sort_values(ascending=False)

    fig = plt.figure(figsize=(6, 9))
    ax = fig.add_subplot(111)
    cv.plot.barh(color='k', ax=ax)
    fig.subplots_adjust(bottom=0.07, left=0.4, right=0.98, top=0.98)
    ax.set_xlabel('CV across samples')
    fig.savefig(os.path.join(outdir, "cv_across_samples.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "cv_across_samples.tiff"), dpi=200)

    # load ssGSEA normalised scores
    # NB normalised scores are just raw scores rescaled to [0, 1] (separately by pathway)
    ssgsea = pd.read_csv(ssgsea_norm_fn, skiprows=2, header=0, index_col=0, sep='\t').drop('Description', axis=1)

    # load true names for IPA pathways
    ssgsea_pathway_names = pd.read_csv(ssgsea_pathways_fn, index_col=0, header=None).squeeze().str.decode('utf-8')
    ssgsea.index = ssgsea_pathway_names.reindex(ssgsea.index.str.replace('_', ' '))

    # heatmap: proportions for each patient
    # standardise across columns, because each cell type has different mean proportion
    rl = hc.linkage(xcell_prop.astype(float).transpose(), method='average', metric='euclidean')

    cg = clustering.plot_clustermap(
        xcell_prop.astype(float).transpose(),
        metric='euclidean',
        show_gene_labels=False,
        show_gene_clustering=True,
        cmap='YlOrRd',
        row_linkage=rl,
        z_score=1,
        vmin=-1.5,
        vmax=6.,
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

    # load IPA signatures
    ipa_sign_fn = os.path.join(HGIC_LOCAL_DIR, 'current', 'input_data', 'ipa_pathways', 'ipa_exported_pathways_ensembl_ids.csv')
    ipa_signatures = load_ipa_signatures(ipa_sign_fn)

    # carry out correlation analysis

    # 1. rename FFPE samples to patient ID for easier linking
    # metafile makes this easier
    m_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/our_data/rnaseq/bulk_ffpe/salmon_pipeline',
        'ffpe_salmon_meta.xlsx'
    )
    m = pd.read_excel(m_fn, header=0, index_col=0)

    new_cols = m.loc[xcell_prop.columns, 'reference_id']
    xcell_prop.columns = new_cols

    # 2. aggregate over ssGSEA replicates

    ssgsea_matched = ssgsea.copy()
    cg = clustering.plot_clustermap(ssgsea_matched, cmap='RdBu_r')
    cg.gs.update(bottom=0.15)

    # this identifies that GBM026_P8 is an outlier - remove?
    ssgsea_matched = ssgsea_matched.drop('GBM026_P8', axis=1)

    # take mean
    pp = ssgsea_matched.columns.str.replace(r'_P.*', '').str.replace('GBM', '')
    ssgsea_matched = ssgsea_matched.groupby(pp, axis=1).mean()

    co, co_p = analyse_xcell_results.pathway_cell_type_composition_correlation_analysis(
        ssgsea_matched,
        xcell_prop,
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
        ipa_signatures_symb[k] = references.ensembl_to_gene_symbol(v).dropna()

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
    if pct_shared_aggr is not None:
        hatch_df = pct_shared_aggr.loc[co.columns, co.index] > pct_shared_max
    else:
        hatch_df = None

    plot_dict = analyse_xcell_results.plot_heatmap_with_quantification(
        co,
        co_p,
        alpha=corr_alpha,
        hatch_df=hatch_df,
        figsize=(8., 9.)
    )
    plt.setp(plot_dict['main_ax'].yaxis.get_ticklabels(), fontsize=8)
    gs = plot_dict['gs']
    gs.update(left=0.45, bottom=0.2, top=0.99, right=0.93, wspace=0.03)
    fig = plot_dict['fig']
    fig.savefig(os.path.join(outdir, "cell_proportion_pathway_ssgsea_%s_clustering.png" % corr_metric), dpi=200)
    fig.savefig(os.path.join(outdir, "cell_proportion_pathway_ssgsea_%s_clustering.tiff" % corr_metric), dpi=200)
    fig.savefig(os.path.join(outdir, "cell_proportion_pathway_ssgsea_%s_clustering.pdf" % corr_metric), dpi=200)


