import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from adjustText import adjust_text
import pandas as pd
import itertools
import collections
import os

from rnaseq import loader, filter, general
from scripts.hgic_final import consts
from plotting import common, pca, _plotly, scatter, adjuster
from utils import output, log, setops, excel
from stats import decomposition
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
        annotate_features_radius=None,
        annotate_features_quantile=None,
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
    if annotate_features_radius is not None and annotate_features_quantile is not None:
        raise AttributeError("Supply EITHER annotate_features_radius OR annotate_features_quantile.")

    if annotate_features_quantile is not None:
        assert 0 < annotate_features_quantile < 1, "annotate_features_quantile must be between 0 and 1 (not inclusive)."

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

    selected = None

    if annotate_features_radius is not None:
        # annotate most influential genes
        selected = pca.highlight_biplot_features(feat_x, feat_y, annotate_features_radius, ax)

    if annotate_features_quantile is not None:
        rad = (feat_x ** 2 + feat_y ** 2) ** .5
        cut = sorted(rad)[int(len(rad) * annotate_features_quantile)]
        selected = rad >= cut

    if selected is not None:
        genes_selected = dat.index[selected]
        symbols_selected = references.ensembl_to_gene_symbol(genes_selected)

        # add gene symbol annotations
        text_handles = []
        for ix, gs in zip(np.where(selected)[0], symbols_selected):
            if not pd.isnull(gs):
                text_handles.append(ax.text(feat_x[ix], feat_y[ix], gs, zorder=10))
        # rearrange them to avoid overlaps
        if adjust_annotation:
            adjuster.adjust_text_radial_plus_repulsion(text_handles, **adjust_annotation_kwargs)

    return fig, ax, res


def get_topmost_quantile_by_loading(feat_dat, quantile):
    feat_dat = pd.DataFrame(feat_dat)

    # extract gene lists for pathway analysis
    rad = (feat_dat ** 2).sum(axis=1) ** .5
    val = np.quantile(rad, quantile)
    return dat.index[rad >= val]


def extract_gois(feat_dat, de_res, quantile, mean_logfc=True, alpha=None):
    """
    Extract genes of interest from the SVD representation, based on an upper quantile of highest loadings.
    These are then combined with DE data.
    The selection is based on the Euclidean distance from the origin across all included component loadings.
    :param feat_dat: pandas DataFrame (or castable to one), with genes on the rows and however many components required
    on the columns.
    :param de_res:
    :param quantile:
    :param mean_logfc:
    :param alpha:
    :return:
    """
    # feat_dat = pd.DataFrame(feat_dat)

    # extract gene lists for pathway analysis
    # rad = (feat_dat ** 2).sum(axis=1) ** .5
    # val = np.quantile(rad, quantile)
    # ens = dat.index[rad >= val]
    ens = get_topmost_quantile_by_loading(feat_dat, quantile)
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
    obj = loader.load_by_patient(consts.PIDS, source='salmon')
    obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES_GIC + consts.S1_RNASEQ_SAMPLES_INSC)
    obj.meta.insert(0, 'patient_id', obj.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    cmap = common.get_best_cmap(len(consts.PIDS))
    scatter_colours = dict(zip(consts.PIDS, cmap))

    scatter_markers = {
        'GBM': 's',
        'iNSC': 'o'
    }

    # scaling parameter applied during SVD
    scale_preserved = 0.05

    sample_colours = obj.meta.patient_id.map(scatter_colours.get).to_dict()
    sample_markers = obj.meta.type.map(scatter_markers.get).to_dict()

    dat = filter.filter_by_cpm(obj.data, min_n_samples=2)
    # TODO: include VST or similar here
    dat = np.log(dat + eps)
    # copy of dat with gene symbols
    dat_with_gs = dat.copy()
    general.add_gene_symbols_to_ensembl_data(dat_with_gs)
    # fill back in with ENS where no gene symbol is available
    dat_with_gs.loc[dat_with_gs['Gene Symbol'].isnull(), 'Gene Symbol'] = dat_with_gs.index[dat_with_gs['Gene Symbol'].isnull()]

    # recreate Sven's original (unweighted) plot
    # we're not going to use this, it's just for validation / reproducibility
    fig, ax, _ = plot_biplot(
        dat,
        obj.meta,
        (0, 1),
        scatter_colours,
        scatter_markers,
        annotate_features_radius=0.4,
        include_weighting=False,
        scale=10.
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated_unweighted.png" % (0, 1)), dpi=200)

    # dict to store 'selected' genes for later use
    selected_by_quantile_mean_logfc = {}
    selected_by_quantile_separate_logfc = {}
    selected_by_quantile_gene_only_all = {}

    svd_res = decomposition.svd_for_biplot(dat, feat_axis=0, scale_preserved=scale_preserved)

    for first_dim in dims:
        dims_pair = (first_dim, first_dim + 1)
        selection_radius = selection_radii_for_plotting[first_dim]

        # extract gene lists for pathway analysis
        for dim in [(first_dim,), dims_pair]:
            this_feat = svd_res['feat_dat'][[i + 1 for i in dim]]

            # mean logFC
            selected_by_quantile_mean_logfc[dim] = extract_gois(
                this_feat,
                de_res,
                quantile,
                alpha=alpha
            )

            # separate logFC
            # remove non-significant results
            selected_by_quantile_separate_logfc[dim] = extract_gois(
                this_feat,
                de_res,
                quantile,
                mean_logfc=False,
                alpha=alpha
            )

        selection_radius = selection_radii_for_plotting[first_dim]
        fig, ax, res = plot_biplot(
            dat,
            obj.meta,
            dims_pair,
            scatter_colours,
            scatter_markers,
            annotate_features_radius=selection_radius,
            scale=0.05
        )
        fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated.png" % dims_pair), dpi=200)

        if publish_plotly:
            size_scaling = [
                [0.1, selection_radius],
                [2., 10.]
            ]
            feature_text = dat_with_gs['Gene Symbol']
            sample_text = dat.columns
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

    ix_all = sorted(setops.reduce_union(*[t.index for t in selected_by_quantile_mean_logfc.values()]))

    ipa_mean_logfc = pd.DataFrame(index=ix_all)
    for k, v in selected_by_quantile_mean_logfc.items():
        ipa_mean_logfc.insert(0, "pc_%s_logFC" % '-'.join([str(t+1) for t in k]), v)
    ipa_mean_logfc.to_excel(os.path.join(outdir, "for_ipa_mean_logfc.xlsx"))

    # for separated data, combine single and paired PC for maximum efficiency
    for first_dim in dims:
        dims_pair = (first_dim, first_dim + 1)
        ix_all = setops.reduce_union(*[selected_by_quantile_separate_logfc[k].index for k in [(first_dim,), dims_pair]])
        this_df = pd.DataFrame(index=ix_all)
        for k in [(first_dim,), dims_pair]:
            tt = selected_by_quantile_separate_logfc[k].copy()
            tt = tt.loc[:, tt.columns.str.contains('logFC')]
            tt.columns = tt.columns.str.replace('_logFC', '_%s_logFC' % '-'.join([str(t+1) for t in k]))
            this_df = pd.concat((this_df, tt), axis=1, sort=True)
        this_df.to_excel(os.path.join(outdir, "for_ipa_separate_logfc_pc%d.xlsx" % (first_dim + 1)))

    # combine with DE results and export to table
    for_export = {}
    for first_dim in dims:
        dims_pair = (first_dim, first_dim + 1)
        for dim in [(first_dim,), dims_pair]:
            the_key = "PC_%s" % '-'.join([str(t + 1) for t in dim])
            this_feat = svd_res['feat_dat'][[i + 1 for i in dim]]
            this_ens = get_topmost_quantile_by_loading(this_feat, quantile).intersection(de_res.index)
            for_export[the_key] = de_res.loc[this_ens]
    excel.pandas_to_excel(for_export, os.path.join(outdir, "full_de_syngeneic_only_filtered_by_biplot.xlsx"))