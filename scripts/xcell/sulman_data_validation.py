import csv
import os

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy as hc

from plotting import clustering
from scripts.hgic_final import analyse_xcell_results
from settings import HGIC_LOCAL_DIR, GIT_LFS_DATA_DIR
from utils import log, output

logger = log.get_console_logger()


XCELL_SIGNATURE_FN = os.path.join(GIT_LFS_DATA_DIR, 'xcell', 'ESM3_signatures.xlsx')
SYN_INDIR = os.path.join(
    HGIC_LOCAL_DIR,
    'current/core_pipeline/rnaseq/merged_s1_s2/ipa/xcell'
)


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
        'current/input_data/verhaak_sulman/xcell/xcell_results_rnaseq_fpkm.xlsx'
    )
    ssgsea_norm_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/verhaak_sulman/ssgsea/norm.score.txt'
    )
    ssgsea_raw_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/verhaak_sulman/ssgsea/raw.score.txt'
    )
    ssgsea_pathways_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/verhaak_sulman/ssgsea/ipa_pathway_numbering.csv'
    )
    syngeneic_correlation_fn = os.path.join(SYN_INDIR, 'correlation_%s_syngeneic.xlsx' % corr_metric)

    outdir = output.unique_output_dir()

    xcell_res = pd.read_excel(xcell_fn, sheet_name=None)
    xcell_prop = xcell_res['Proportions'].copy()
    xcell_pval = xcell_res['P values'].copy()

    # drop the 'scores' rows
    xcell_prop = xcell_prop.loc[xcell_pval.index]

    # erase non-significant results if requested
    if zero_nonsign_xcell:
        xcell_prop.loc[xcell_pval > xcell_alpha] = 0.

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
    ipa_sign_fn = os.path.join(HGIC_LOCAL_DIR, 'current', 'input_data', 'ipa_pathways', 'ipa_exported_pathways_symbols.csv')
    ipa_signatures_symb = load_ipa_signatures(ipa_sign_fn)

    # carry out correlation analysis
    # reduce ssGSEA results down to GIC only, then match column names (even though this means they're labelled as bulk)
    ssgsea_matched = ssgsea.drop(xcell_prop.columns, axis=1)
    ssgsea_matched.columns = ssgsea_matched.columns.str.replace('PC', 'TU')

    # reduce ssGSEA results to include only pathways identified in the syngeneic comparisons
    # load syngeneic IPA results
    res_syn = pd.read_excel(syngeneic_correlation_fn, sheet_name=None)
    ssgsea_matched = ssgsea_matched.reindex(res_syn[corr_metric].columns).dropna()

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

    # compute overlap between cell type signatures and IPA signatures
    pct_shared = analyse_xcell_results.compute_cell_type_pathway_overlap(
        ipa_signatures_symb,
        xcell_signatures,
    )

    # aggregate taking max over pathways
    cc = pct_shared.columns.str.replace(r'(?P<ct>[^_]*)_.*', r'\g<ct>')
    pct_shared_aggr = pct_shared.groupby(cc, axis=1).max()

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
        figsize=(8., 11.)
    )
    plt.setp(plot_dict['main_ax'].yaxis.get_ticklabels(), fontsize=8)
    gs = plot_dict['gs']
    gs.update(left=0.45, bottom=0.17, top=0.99, right=0.93, wspace=0.03)
    fig = plot_dict['fig']
    fig.savefig(os.path.join(outdir, "cell_proportion_pathway_ssgsea_%s_clustering.png" % corr_metric), dpi=200)
    fig.savefig(os.path.join(outdir, "cell_proportion_pathway_ssgsea_%s_clustering.tiff" % corr_metric), dpi=200)
    fig.savefig(os.path.join(outdir, "cell_proportion_pathway_ssgsea_%s_clustering.pdf" % corr_metric), dpi=200)


