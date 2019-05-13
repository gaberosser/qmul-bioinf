import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import itertools
import collections
import os

from rnaseq import loader, filter, general
from methylation import loader as methylation_loader, process
from scripts.hgic_final import consts
from plotting import common, pca, _plotly, scatter, adjuster
from utils import output, log, setops, excel, dictionary
from stats import decomposition, transformations
import references
from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR

import plotly.plotly as py
from plotly import graph_objs as go


logger = log.get_console_logger()

def generate_plotly_plot(
        res,
        filename,
        feature_size_scaling,
        feature_text,
        sample_text,
        sample_colours,
        sample_markers,
        components=(1, 2),
        sample_marker_size=12.,
        auto_open=False,
        feature_text_mask=None
):
    """

    :param res:
    :param filename:
    :param feature_size_scaling:
    :param feature_text:
    :param sample_text:
    :param sample_colours:
    :param sample_markers:
    :param sample_marker_size:
    :param auto_open:
    :param feature_text_mask: If supplied, this should be a boolean array the same length as feature_text. Where it is
     True, the corresponding text will _not_ be drawn (but is still available on hover).
     If absent, all text is hidden and only available on hover.
    :return:
    """
    x, y = [res['feat_dat'][i] for i in components]

    msize_in = (x ** 2 + y ** 2) ** .5

    tmp = (msize_in - feature_size_scaling[0][0]) / (feature_size_scaling[0][1] - feature_size_scaling[0][0])
    tmp[tmp < 0] = 0
    tmp[tmp > 1] = 1
    msize = tmp * (feature_size_scaling[1][1] - feature_size_scaling[1][0]) + feature_size_scaling[1][0]

    if feature_text_mask is None:
        feature_text_mask = np.ones(len(x), dtype=bool)

    # plot into two groups: text shown permanently and text hidden
    feat_trace = go.Scatter(
        x=x[feature_text_mask],
        y=y[feature_text_mask],
        mode='markers',
        text=feature_text[feature_text_mask],
        marker={
            'size': msize[feature_text_mask],
            'color': 'rgb(155, 155, 155, .7)'
        }
    )

    feat_trace_text_shown = go.Scatter(
        x=x[~feature_text_mask],
        y=y[~feature_text_mask],
        mode='markers+text',
        text=feature_text[~feature_text_mask],
        marker={
            'size': msize[~feature_text_mask],
            'color': 'rgb(155, 155, 155, .7)'
        }
    )

    plotly_markers = np.array([_plotly.MATPLOTLIB_TO_PLOTLY_MARKERS[sample_markers[t]] for t in sample_text])

    sx, sy = [res['sample_dat'][i] for i in components]

    sample_trace = go.Scatter(
        x=sx,
        y=sy,
        mode='markers+text',
        text=sample_text,
        textposition='bottom center',
        textfont={
            'color': [sample_colours[t] for t in sample_text]
        },
        marker={
            'size': sample_marker_size,
            'color': [sample_colours[t] for t in sample_text],
            'symbol': plotly_markers,
            'line': {'color': 'black', 'width': 2.}
        }
    )

    ev_pct = res['explained_variance_ratio'] * 100.

    layout = go.Layout(
        showlegend=False,
        hovermode='closest',
        width=800,
        height=600,
        xaxis=go.layout.XAxis(title='PC %d (%.2f %%)' % (components[0] - 1, ev_pct[components[0] - 1])),
        yaxis=go.layout.YAxis(title='PC %d (%.2f %%)' % (components[1] - 1, ev_pct[components[1] - 1])),
        dragmode='pan'
    )
    fig = go.Figure(data=[feat_trace, feat_trace_text_shown, sample_trace], layout=layout)
    p = py.plot(fig, filename=filename, auto_open=auto_open)

    return p


def plot_biplot(
        dat,
        meta,
        dims,
        scatter_colours,
        scatter_markers,
        annotate_features=None,
        adjust_annotation=True,
        adjust_annotation_kwargs=None,
        **kwargs
):
    """
    :param dat:
    :param meta: pd.DataFrame, must have columns entitled `type` and `patient_id`
    :param dims:
    :param scatter_colours:
    :param scatter_markers:
    :param annotate_features_radius: If supplied, this is the biplot radius outside of which we annotate genes (by
    symbol).
    :param **kwargs: Passed to pca.biplot()
    :return:
    """

    if adjust_annotation_kwargs is None:
        adjust_annotation_kwargs = {}

    sample_colours = meta.patient_id.map(scatter_colours.get).to_dict()
    sample_markers = meta.type.map(scatter_markers.get).to_dict()

    res = pca.biplot(
        dat,
        plot_dims=dims,
        sample_colours=sample_colours,
        sample_markers=sample_markers,
        **kwargs
    )

    sample_x, sample_y = res['sample_data']
    feat_x, feat_y = res['feature_data']
    ax = res['ax']
    fig = res['fig']

    typ_ix, typ = meta.type.factorize()

    # connect patients
    for pid in meta.patient_id.unique():
        ix = meta.patient_id == pid
        for t0, t1 in itertools.combinations(typ, 2):
            # draw all possible connections between these two cell types (in one direction only)
            ix0 = meta.index[ix & (meta.type == t0)]
            ix1 = meta.index[ix & (meta.type == t1)]

            for a, b in itertools.product(ix0, ix1):
                ax.plot(
                    [sample_x[meta.index == a][0], sample_x[meta.index == b][0]],
                    [sample_y[meta.index == a][0], sample_y[meta.index == b][0]],
                    lw=1.5,
                    color=scatter_colours[pid],
                    zorder=9
                )

    # custom legend outside of plot
    line_kwargs = {
        'class': 'line',
        'markerfacecolor': 'none',
        'markeredgecolor': 'k',
        'markeredgewidth': 1.0,
        'linestyle': 'none'
    }
    patch_kwargs = {
        'class': 'patch',
        'edgecolor': 'k',
        'linewidth': 1.
    }
    legend_dict = {'Patient': collections.OrderedDict(), 'Cell type': collections.OrderedDict()}
    for pid in consts.PIDS:
        ll = dict(patch_kwargs)
        ll['facecolor'] = scatter_colours[pid]
        legend_dict['Patient'][pid] = ll
    for t in typ:
        pp = dict(line_kwargs)
        pp['marker'] = scatter_markers[t]
        legend_dict['Cell type'][t] = pp

    res['legend_dict'] = legend_dict

    common.add_custom_legend(ax, legend_dict, loc_outside=True)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)

    if annotate_features is not None:
        symbols_selected = references.ensembl_to_gene_symbol(annotate_features)

        # add gene symbol annotations
        text_handles = []
        for ix, gs in zip(annotate_features, symbols_selected):
            if not pd.isnull(gs):
                text_handles.append(ax.text(feat_x[ix], feat_y[ix], gs, zorder=10))
        # rearrange them to avoid overlaps
        if adjust_annotation:
            adjuster.adjust_text_radial_plus_repulsion(text_handles, **adjust_annotation_kwargs)

    return fig, ax, res


def inverse_covariance_projection(feat_dat):
    """
    Project data into the covariance representation (square rooting weights) then return the Euclidean norm.
    Key idea: quantify the distance from the origin of data points in 2D where there is some clear covariance structure.
    :param feat_dat: 2D matrix, either N x 1 (in which case the transformation is trivial) or N x 2.
    :return:
    """
    # distance based on projection into covariance matrix
    x = np.array(feat_dat)
    # centre
    xm = x.mean(axis=0)
    x0 = x - xm

    n = x.shape[1]
    cov = np.cov(x0.transpose())
    if n == 1:
        y = cov ** -0.5 * x0
        dist = y.squeeze()
    elif n == 2:
        A = np.linalg.cholesky(cov)
        y = np.linalg.solve(A, x0.T)
        dist = (y.T ** 2).sum(axis=1) ** .5
    else:
        raise ValueError("Input data must be an N x 1 or N x 2 matrix (or similar).")

    if isinstance(feat_dat, pd.DataFrame):
        dist = pd.Series(dist, index=feat_dat.index)
    return dist


def get_topmost_quantile_by_loading(feat_dat, quantile):
    feat_dat = pd.DataFrame(feat_dat)

    # extract gene lists for pathway analysis
    rad = (feat_dat ** 2).sum(axis=1) ** .5
    val = np.quantile(rad, quantile)
    return feat_dat.index[rad >= val]


def extract_gois(loading_scale, de_res, quantile, mean_logfc=True, alpha=None):
    """
    Extract genes of interest based on an upper quantile of highest loadings.
    These are then combined with DE data.
    The selection is based on the Euclidean distance from the origin across all included component loadings.
    :param loading_scale: pandas Series (or castable to one), with entries corresponding to features (genes) and
    representing a summary (e.g. L2 norm) of the magnitude.
    :param de_res:
    :param quantile:
    :param mean_logfc:
    :param alpha:
    :return:
    """
    val = np.quantile(loading_scale, quantile)
    ens = pd.Series(loading_scale).index[loading_scale >= val]
    # ens = get_topmost_quantile_by_loading(feat_dat, quantile)
    dif = ens.difference(de_res.index)
    if len(dif):
        logger.warning(
            "%d of the genes selected on the biplot are NOT in the DE results: %s",
            len(dif),
            ', '.join(dif)
        )
    ens = ens.intersection(de_res.index)
    this_de_res = de_res.loc[ens]

    this_logfc = this_de_res[["%s_logFC" % p for p in pids]]

    if alpha is not None:
        this_fdr = this_de_res[["%s_FDR" % p for p in pids]].fillna(1.)
        ens_keep = (this_fdr <= alpha).sum(axis=1) > 0
        this_logfc = this_logfc.loc[ens_keep]
        this_fdr = this_fdr.loc[ens_keep]
        ens = ens[ens_keep]
        logger.info("Removing %d genes that were not significantly DE in any comparison", (~ens_keep).sum())


    if mean_logfc:
        this_mean_logfc = np.nanmean(this_logfc, axis=1)
        # are there any non-concordant genes? if so, warn and remove
        non_conc = np.array(
            [np.sign(v[v.abs() > min_logfc]).dropna().unique().size > 1 for _, v in this_logfc.iterrows()])
        if non_conc.any():
            logger.warning(
                "%d selected genes were not concordant across DE results: %s. We'll set the logFC to blank for these.",
                non_conc.sum(),
                ens[non_conc]
            )
            this_mean_logfc[non_conc] = np.nan
        return pd.Series(this_mean_logfc, index=ens)
    else:
        res = this_logfc
        if alpha is not None:
            to_drop = (this_fdr > alpha)
            this_fdr.values[to_drop] = np.nan
            this_logfc.values[to_drop] = np.nan
            res = pd.concat((this_logfc, this_fdr), axis=1)
        return res


if __name__ == '__main__':
    """
    Idea here: recreate the analysis Sven carried out, generating biplots for the RNA-Seq data.
    We can then extend this idea to methylation (?)

    Absolutely amazing resource here:
    https://www.fbbva.es/microsite/multivariate-statistics/biplots.html

    Code snippet inspiration here:
    https://stackoverflow.com/questions/39216897/plot-pca-loadings-and-loading-in-biplot-in-sklearn-like-rs-autoplot
    """
    outdir = output.unique_output_dir()
    outdir_rna = os.path.join(outdir, "rnaseq")
    outdir_meth = os.path.join(outdir, "methylation")
    outdir_joint = os.path.join(outdir, "joint")
    for d in [outdir_rna, outdir_meth, outdir_joint]:
        if not os.path.exists(d):
            os.makedirs(d)

    publish_plotly = True
    pids = consts.PIDS

    alpha = 0.05
    eps = 1.  # offset applied during log transform

    min_logfc = consts.DE_PARAMS['lfc']

    # quantile to use for defining genes of interest
    # these lists can then be sent to IPA or similar
    quantile = 0.99

    # dimensions (components) to investigate
    dims = [0, 1, 2]

    # number of genes / probes to select in biplots
    top_n_features_rna = 40
    top_n_features_meth = 40

    selection_radii_for_plotting = {
        0: 0.6,
        1: 0.30,
        2: 0.25
    }

    # path to syngeneic DE results
    fn_de_res = os.path.join(
        HGIC_LOCAL_DIR,
        'current/core_pipeline/rnaseq/',
        'full_de_syngeneic_only.xlsx'
    )

    # load DE results
    de_res = pd.read_excel(fn_de_res, index_col=0)

    # load data for iNSC and GBM (Salmon TPM)
    rna_obj = loader.load_by_patient(consts.PIDS, source='salmon')
    rna_obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES_GIC + consts.S1_RNASEQ_SAMPLES_INSC)
    rna_obj.meta.insert(0, 'patient_id', rna_obj.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    # load methylation data
    meth_obj = methylation_loader.load_by_patient(consts.PIDS, include_control=False)
    meth_obj.filter_by_sample_name(consts.S1_METHYL_SAMPLES_GIC + consts.S1_METHYL_SAMPLES_INSC)
    meth_obj.meta.insert(0, 'patient_id', meth_obj.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    cmap = common.get_best_cmap(len(consts.PIDS))
    scatter_colours = dict(zip(consts.PIDS, cmap))

    scatter_markers = {
        'GBM': 's',
        'iNSC': 'o'
    }

    # scaling parameter applied during SVD
    scale_preserved = 0.05

    ######## RNA-Seq ########

    sample_colours = rna_obj.meta.patient_id.map(scatter_colours.get).to_dict()
    sample_markers = rna_obj.meta.type.map(scatter_markers.get).to_dict()

    dat_rna = filter.filter_by_cpm(rna_obj.data, min_n_samples=2)

    dat_rna = np.log(dat_rna + eps)
    # dat_rna = transformations.variance_stabilizing_transform(dat_rna)
    # copy of dat with gene symbols
    dat_with_gs_rna = dat_rna.copy()
    general.add_gene_symbols_to_ensembl_data(dat_with_gs_rna)
    # fill back in with ENS where no gene symbol is available
    dat_with_gs_rna.loc[dat_with_gs_rna['Gene Symbol'].isnull(), 'Gene Symbol'] = dat_with_gs_rna.index[dat_with_gs_rna['Gene Symbol'].isnull()]

    # recreate Sven's original (unweighted) plot
    # we're not going to use this, it's just for validation / reproducibility
    # fig, ax, _ = plot_biplot(
    #     dat,
    #     rna_obj.meta,
    #     (0, 1),
    #     scatter_colours,
    #     scatter_markers,
    #     annotate_features_radius=0.4,
    #     include_weighting=False,
    #     scale=10.
    # )
    # fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated_unweighted.png" % (0, 1)), dpi=200)

    # dict to store 'selected' genes for later use
    rna_selected_by_quantile_mean_logfc = {}
    rna_selected_by_quantile_separate_logfc = {}

    svd_res = decomposition.svd_for_biplot(dat_rna, feat_axis=0, scale_preserved=scale_preserved)

    for first_dim in dims:
        dims_pair = (first_dim, first_dim + 1)
        selection_radius = selection_radii_for_plotting[first_dim]

        # extract gene lists for pathway analysis
        for dim in [(first_dim,), dims_pair]:
            this_feat = svd_res['feat_dat'][[t + 1 for t in dim]]
            this_dist = inverse_covariance_projection(this_feat)

            # mean logFC
            rna_selected_by_quantile_mean_logfc[dim] = extract_gois(
                this_dist,
                de_res,
                quantile,
                alpha=alpha
            )

            # separate logFC
            # remove non-significant results
            rna_selected_by_quantile_separate_logfc[dim] = extract_gois(
                this_dist,
                de_res,
                quantile,
                mean_logfc=False,
                alpha=alpha
            )

        # selection_radius = selection_radii_for_plotting[first_dim]
        ix = this_dist.sort_values(ascending=False).index[:top_n_features_rna]
        fig, ax, res = plot_biplot(
            dat_rna,
            rna_obj.meta,
            dims_pair,
            scatter_colours,
            scatter_markers,
            annotate_features=this_dist.sort_values(ascending=False).index[:top_n_features_rna],
            scale=0.05
        )
        fig.savefig(os.path.join(outdir_rna, "pca_biplot_dims_%d-%d_annotated.png" % dims_pair), dpi=200)

        fig, ax, res = plot_biplot(
            dat_rna,
            rna_obj.meta,
            dims_pair,
            scatter_colours,
            scatter_markers,
            scale=0.05
        )
        fig.savefig(os.path.join(outdir_rna, "pca_biplot_dims_%d-%d_unannotated.png" % dims_pair), dpi=200)

        if publish_plotly:
            size_scaling = [
                [0.1, selection_radius],
                [2., 10.]
            ]
            feature_text = dat_with_gs_rna['Gene Symbol']
            sample_text = dat_rna.columns
            rad = (np.array(zip(*[svd_res['feat_dat'][i + 1] for i in dims_pair])) ** 2).sum(axis=1) ** .5
            to_annotate = rad > selection_radius
            p1 = generate_plotly_plot(
                svd_res,
                filename="pca_biplot_dims_%d-%d" % tuple(t + 1 for t in dims_pair),
                feature_size_scaling=size_scaling,
                feature_text=feature_text,
                sample_text=sample_text,
                sample_colours=sample_colours,
                sample_markers=sample_markers,
                feature_text_mask=~to_annotate,
                components=tuple(i + 1 for i in dims_pair),
            )

    # export lists for IPA

    ix_all = sorted(setops.reduce_union(*[t.index for t in rna_selected_by_quantile_mean_logfc.values()]))

    ipa_mean_logfc = pd.DataFrame(index=ix_all)
    for k, v in rna_selected_by_quantile_mean_logfc.items():
        ipa_mean_logfc.insert(0, "pc_%s_logFC" % '-'.join([str(t+1) for t in k]), v)
    ipa_mean_logfc.to_excel(os.path.join(outdir, "for_ipa_mean_logfc.xlsx"))

    # for separated data, combine single and paired PC for maximum efficiency
    for first_dim in dims:
        dims_pair = (first_dim, first_dim + 1)
        ix_all = setops.reduce_union(*[rna_selected_by_quantile_separate_logfc[k].index for k in [(first_dim,), dims_pair]])
        this_df = pd.DataFrame(index=ix_all)
        for k in [(first_dim,), dims_pair]:
            tt = rna_selected_by_quantile_separate_logfc[k].copy()
            tt = tt.loc[:, tt.columns.str.contains('logFC')]
            tt.columns = tt.columns.str.replace('_logFC', '_%s_logFC' % '-'.join([str(t+1) for t in k]))
            this_df = pd.concat((this_df, tt), axis=1, sort=True)
        this_df.to_excel(os.path.join(outdir_rna, "for_ipa_separate_logfc_pc%d.xlsx" % (first_dim + 1)))

    # combine with DE results and export to table
    for_export = {}
    for first_dim in dims:
        dims_pair = (first_dim, first_dim + 1)
        for dim in [(first_dim,), dims_pair]:
            the_key = "PC_%s" % '-'.join([str(t + 1) for t in dim])
            this_feat = svd_res['feat_dat'][[i + 1 for i in dim]]
            this_ens = get_topmost_quantile_by_loading(this_feat, quantile).intersection(de_res.index)
            for_export[the_key] = de_res.loc[this_ens]
    excel.pandas_to_excel(for_export, os.path.join(outdir_rna, "full_de_syngeneic_only_filtered_by_biplot.xlsx"))

    ######## DNA methylation ########

    sample_colours = meth_obj.meta.patient_id.map(scatter_colours.get).to_dict()
    sample_markers = meth_obj.meta.type.map(scatter_markers.get).to_dict()

    dat = process.m_from_beta(meth_obj.data)

    probe_prop = 0.01

    # dict to store 'selected' probes for later use
    meth_selected_by_quantile_mean_logfc = {}
    meth_selected_by_quantile_separate_logfc = {}
    meth_selected_by_quantile_gene_only_all = {}

    svd_res = decomposition.svd_for_biplot(dat, feat_axis=0, scale_preserved=scale_preserved)

    # downsample features (probes) at random for plotting purposes
    probe_ix = dat.index.to_native_types()
    np.random.shuffle(probe_ix)
    probe_ix = probe_ix[:int(len(probe_ix) * probe_prop)]

    for first_dim in dims:
        dims_pair = (first_dim, first_dim + 1)

        fig, ax, res = plot_biplot(
            dat,
            meth_obj.meta,
            dims_pair,
            scatter_colours,
            scatter_markers,
            scale=0.05,
            loading_idx=probe_ix
        )
        fig.savefig(os.path.join(outdir_meth, "pca_biplot_dims_%d-%d.png" % dims_pair), dpi=200)

    # hypothesis: PC2 identifies subgroups
    subgroups = consts.SUBGROUPS
    subgroups_inv = dictionary.complement_dictionary_of_iterables(subgroups, squeeze=True)
    groups = meth_obj.meta.patient_id.apply(subgroups_inv.get)

    colours = dict([(k, consts.SUBGROUP_SET_COLOURS["%s partial" % k]) for k in groups.unique()])

    ix = meth_obj.meta.type == 'GBM'

    for dim in range(1, 4):

        fig, ax = plt.subplots(figsize=(6.5, 3.2))
        sns.boxplot(svd_res['sample_dat'][dim].loc[ix], y=groups[ix], ax=ax, palette=colours)
        sns.swarmplot(
            svd_res['sample_dat'][dim].loc[ix],
            y=groups[ix],
            ax=ax,
            color='0.3',
            edgecolor='k',
            linewidth=1.0,
            alpha=0.5
        )
        ax.set_ylabel('Subgroup')
        ax.set_xlabel('PC%d' % dim)
        fig.tight_layout()

        x = svd_res['sample_dat'][dim].loc[(meth_obj.meta.type == 'GBM') & (groups == 'RTK I')]
        y = svd_res['sample_dat'][dim].loc[(meth_obj.meta.type == 'GBM') & (groups == 'RTK II')]

        u_test = stats.mannwhitneyu(x, y)
        print "MWU test of PC%d value (RTK I vs RTK II): p=%.3e" % (dim, u_test.pvalue)

        fig.savefig(os.path.join(outdir_meth, "PC%d_vs_subgroup.png" % dim), dpi=200)

    # similar, but now focus on the features (probes)
    # for each PC, get the top N probes by absolute loading and test distribution
    top_n_probe = 10

    # fig, axs = plt.subplots(nrows=top_n_probe, ncols=3, sharex=True, figsize=(11, 9))

    for dim in range(1, 4):
        this_feat = svd_res['feat_dat'][dim]
        ix2 = this_feat.abs().sort_values(ascending=False).index[:top_n_probe]
        ix2 = this_feat.loc[ix2].sort_values().index
        x = dat.loc[ix2, ix]

        fig, axs = plt.subplots(nrows=top_n_probe, figsize=(3.5, 6.5), sharex=True)
        for i, k in enumerate(ix2):
            ax = axs[i]
            sns.boxplot(x.loc[k], y=groups[ix], ax=ax, palette=colours)
            sns.swarmplot(
                x.loc[k],
                y=groups[ix],
                ax=ax,
                color='0.3',
                edgecolor='k',
                linewidth=1.0,
                alpha=0.5
            )
            ax.set_ylabel(k, fontsize=8)
            ax.xaxis.label.set_visible(False)
        ax.set_xlabel('PC%d' % dim)
        fig.subplots_adjust(bottom=0.06, top=0.98, left=0.2, right=0.98, hspace=0.1)
