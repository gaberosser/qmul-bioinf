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
from utils import output, log, setops
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
        text=dat.columns,
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

    typ_ix, typ = obj.meta.type.factorize()

    # connect patients
    for pid in meta.patient_id.unique():
        ix = obj.meta.patient_id == pid
        for t0, t1 in itertools.combinations(typ, 2):
            # draw all possible connections between these two cell types (in one direction only)
            ix0 = obj.meta.index[ix & (obj.meta.type == t0)]
            ix1 = obj.meta.index[ix & (obj.meta.type == t1)]

            for a, b in itertools.product(ix0, ix1):
                ax.plot(
                    [sample_x[obj.meta.index == a][0], sample_x[obj.meta.index == b][0]],
                    [sample_y[obj.meta.index == a][0], sample_y[obj.meta.index == b][0]],
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


def extract_gois(feat_dat, de_res, quantile, mean_logfc=True, alpha=None):
    # extract gene lists for pathway analysis
    rad = (np.array(zip(*feat_dat)) ** 2).sum(axis=1) ** .5
    val = np.quantile(rad, quantile)
    ens = dat.index[rad >= val]
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

    if alpha is not None:
        ens_keep = (this_fdr <= alpha).sum(axis=1) > 0
        this_logfc = this_logfc.loc[ens_keep]
        this_fdr = this_fdr.loc[ens_keep]
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

    # quantiles to use for defining genes of interest
    # these lists can then be sent to IPA or similar
    quantiles = [0.99, 0.995]

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
        for q in quantiles:
            for dim in [(first_dim,), (first_dim, first_dim + 1)]:
                this_feat = tuple(svd_res['feat_dat'][i + 1] for i in dim)

                # mean logFC
                selected_by_quantile_mean_logfc.setdefault(dim, {})[q] = extract_gois(
                    this_feat,
                    de_res,
                    q
                )
                # separate logFC
                selected_by_quantile_separate_logfc.setdefault(dim, {})[q] = extract_gois(
                    this_feat,
                    de_res,
                    q,
                    mean_logfc=False
                )

                # joint gene only
                selected_by_quantile_gene_only_all.setdefault(dim, {})[q] = extract_gois(
                    this_feat,
                    de_res,
                    q,
                ).dropna().index

            fig, ax, res = plot_biplot(
                dat,
                obj.meta,
                dims_pair,
                scatter_colours,
                scatter_markers,
                annotate_features_quantile=q,
                scale=0.05
            )
            fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_q%.3f.png" % (first_dim, first_dim + 1, q)), dpi=200)

        selection_radius = selection_radii_for_plotting[first_dim]
        fig, ax, res = plot_biplot(
            dat,
            obj.meta,
            (first_dim, first_dim + 1),
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
                feature_text_mask=~to_annotate
            )

    # export lists for IPA
    ix_all = sorted(setops.reduce_union(*[t[0.99].index for t in selected_by_quantile_mean_logfc.values()]))
    ipa_mean_logfc = pd.DataFrame(index=ix_all)
    for k, v in selected_by_quantile_mean_logfc.items():
        ipa_mean_logfc.insert(0, "pc_%s_q99_logFC" % '-'.join([str(t+1) for t in k]), v[0.99])
        ipa_mean_logfc.insert(0, "pc_%s_q995_logFC" % '-'.join([str(t+1) for t in k]), v[0.995])
    ipa_mean_logfc.to_excel(os.path.join(outdir, "for_ipa_mean_logfc.xlsx"))

    for k, v in selected_by_quantile_separate_logfc.items():
        ix_all = setops.reduce_union(*[v[q].index for q in quantiles])
        this = []
        for q in quantiles:
            tt = v[q]
            tt.columns = ["%s_%d_logFC" % (p, int(q * 1000)) for p in pids]
            this.append(tt)
        pd.concat(this, axis=1, sort=True).to_excel(os.path.join(outdir, "for_ipa_separate_logfc_pc%s.xlsx" % '-'.join([str(t+1) for t in k])))

    ipa_gene_only = pd.DataFrame(index=ix_all)
    for k, v in selected_by_quantile_gene_only_all.items():
        for q in quantiles:
            tt = pd.Series(0.001, index=v[q])  # this is a generically significant pvalue
            ipa_gene_only.insert(
                0,
                "pc_%s_q%d_FDR" % (
                    '-'.join([str(t + 1) for t in k]),
                    int(q * 1000)
                ),
                tt)
    ipa_gene_only.fillna(1., inplace=True)
    ipa_gene_only.to_excel(os.path.join(outdir, "for_ipa_gene_only_all.xlsx"))


    # bar chart showing the explained variance
    fig = plt.figure(figsize=(5.6, 3.))
    ax = fig.add_subplot(111)
    n_plot = 6
    ax.bar(range(1, n_plot + 1), res['explained_variance_ratio'][:n_plot] * 100.)
    ax.set_xlabel('Principal component')
    ax.set_ylabel('% variance explained')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "explained_variance_bar_chart.png"), dpi=200)

    # test out a few distinguishing genes
    goi = ['VGF', 'CDKN2A', 'STMN2', 'TMEFF2']
    aa = dat_with_gs.loc[dat_with_gs['Gene Symbol'].isin(goi)].set_index('Gene Symbol')


    def compare_two_gene_levels(
        dat_two_cols,
        meta,
        legend_dict,
        colour_map=scatter_colours,
        marker_map=scatter_markers,
    ):

        ax = scatter.scatter_with_colour_and_markers(
            dat_two_cols,
            colour_subgroups=meta.patient_id,
            colour_map=colour_map,
            marker_subgroups=meta.type,
            marker_map=marker_map
        )
        common.add_custom_legend(ax, legend_dict, loc_outside=True)
        ax.set_xlabel('%s (logTPM)' % dat_two_cols.columns[0])
        ax.set_ylabel('%s (logTPM)' % dat_two_cols.columns[1])
        fig = ax.figure
        fig.tight_layout()
        fig.subplots_adjust(right=0.8)

        return fig, ax

    fig, ax = compare_two_gene_levels(
        aa.transpose()[['VGF', 'CDKN2A']],
        obj.meta,
        res['legend_dict']
    )
    fig.savefig(os.path.join(outdir, "VGF_CDKN2A.png"), dpi=200)

    fig, ax = compare_two_gene_levels(
        aa.transpose()[['VGF', 'TMEFF2']],
        obj.meta,
        res['legend_dict']
    )
    fig.savefig(os.path.join(outdir, "VGF_TMEFF2.png"), dpi=200)

    ## TODO: consider annotating the biplot gene markers based on the DE results
    # For example, highlight patient-specific genes with different colours
    # For example, highlight the most deregulated genes

    dims = (0, 1)  # for copy paste convenience
    fig, ax, res = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        scale=0.05
    )

    feat_dat = pd.DataFrame(np.array(res['feature_data']).transpose(), index=dat.index)
    feat_dat.columns = ['x', 'y']

    # get mean logFC to order genes
    mean_logfc = pd.Series(np.nanmean(de_res[["%s_logFC" % p for p in pids]], axis=1), index=de_res.index)
    mean_logfc.dropna(inplace=True)

    ix = feat_dat.index.intersection(mean_logfc.index)
    mean_logfc = mean_logfc.loc[ix]

    mean_logfc = mean_logfc.loc[mean_logfc.abs().sort_values(ascending=False).index]

    ax.scatter(
        feat_dat.loc[mean_logfc.index[:50], 'x'],
        feat_dat.loc[mean_logfc.index[:50], 'y'],
        c='k',
        facecolor='k',
        marker='^',
    )
    gg = references.ensembl_to_gene_symbol(mean_logfc.index[:50]).dropna()
    for k, v in feat_dat.loc[mean_logfc.index[:50]].iterrows():
        g = gg[k] if k in gg else k
        ax.text(v['x'], v['y'], g)

    dims = (2, 3)  # for copy paste convenience
    fig, ax, res = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        scale=0.05
    )

    feat_dat = pd.DataFrame(np.array(res['feature_data']).transpose(), index=dat.index)
    feat_dat.columns = ['x', 'y']

    # get mean logFC to order genes
    mean_logfc = pd.Series(np.nanmean(de_res[["%s_logFC" % p for p in pids]], axis=1), index=de_res.index)
    mean_logfc.dropna(inplace=True)

    ix = feat_dat.index.intersection(mean_logfc.index)
    mean_logfc = mean_logfc.loc[ix]

    mean_logfc = mean_logfc.loc[mean_logfc.abs().sort_values(ascending=False).index]

    ax.scatter(
        feat_dat.loc[mean_logfc.index[:50], 'x'],
        feat_dat.loc[mean_logfc.index[:50], 'y'],
        c='k',
        facecolor='k',
        marker='^',
    )
    gg = references.ensembl_to_gene_symbol(mean_logfc.index[:50]).dropna()

    # repeat but with DE genes present in a few patients

    dims = (2, 3)  # for copy paste convenience
    fig, ax, res = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        scale=0.05
    )

    feat_dat = pd.DataFrame(np.array(res['feature_data']).transpose(), index=dat.index)
    feat_dat.columns = ['x', 'y']

    # get mean logFC to order genes
    ix = dat.loc[:, dat.columns.str.contains('GBM')].std(axis=1).sort_values(ascending=False).index[:50]

    ax.scatter(
        feat_dat.loc[ix, 'x'],
        feat_dat.loc[ix, 'y'],
        c='k',
        facecolor='k',
        marker='^',
    )
    gg = references.ensembl_to_gene_symbol(ix).dropna()
    for k, v in feat_dat.loc[ix].iterrows():
        g = gg[k] if k in gg else k
        ax.text(v['x'], v['y'], g)

    # multi-axis plot showing the distribution of sample values in each component
    res = decomposition.svd_for_biplot(dat)
    ix = obj.meta.type == 'GBM'

    n_plot = 6
    fig, axs = plt.subplots(n_plot, 1, sharex=True)
    big_ax = fig.add_subplot(111, frameon=False)
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.grid(False)

    for i in range(n_plot):
        ax = axs[i]

        sns.kdeplot(
            res['sample_dat'][i + 1][ix],
            color='b',
            ax=ax,
            shade=True,
            label='GBM' if i == 0 else None
        )
        ax.scatter(res['sample_dat'][i + 1][ix], np.zeros(ix.sum()), c='b', marker='|')
        sns.kdeplot(
            res['sample_dat'][i + 1][~ix],
            color='r',
            ax=ax,
            shade=True,
            label = 'iNSC' if i == 0 else None
        )
        ax.scatter(res['sample_dat'][i + 1][~ix], np.zeros((~ix).sum()), c='r', marker='|')

        if i == 0:
            ax.legend(loc='upper right')
    axs[-1].set_xlabel('Sample value in PC')
    big_ax.set_ylabel('Density (a.u.)')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "kde_pc_values.png"), dpi=200)

    ## Run GO analysis using GOAtools
    # TODO: if this works well, move to a module
    from goatools import base
    import wget

    obo_fn = os.path.join(LOCAL_DATA_DIR, 'gene_ontology', 'current', 'go-basic.obo')
    genetogo_fn = os.path.join(LOCAL_DATA_DIR, 'gene_ontology', 'current', 'gene2go')
    genetoens_fn = os.path.join(LOCAL_DATA_DIR, 'gene_ontology', 'current', 'gene2ensembl.gz')
    genetoens_url = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz"

    obo_fn = base.download_go_basic_obo(obo_fn)
    genetogo_fn = base.download_ncbi_associations(genetogo_fn)
    if not os.path.isfile(genetoens_fn):
        logger.info("Downloading RefGene-Ensembl converter from %s, saving to %s.", genetoens_url, genetoens_fn)
        wget.download(genetoens_url, out=genetoens_fn)

    def ens_to_entrez(ens, genetoens_fn):
        gene2ens = pd.read_csv(genetoens_fn, header=0, sep='\t')
        gene2ens = gene2ens.loc[gene2ens['#tax_id'] == 9606]
        conv_df = gene2ens.loc[gene2ens.Ensembl_gene_identifier.isin(ens), ['GeneID', 'Ensembl_gene_identifier']]
        # reduce to unique (Entrez ID, Ensembl ID) pairs
        conv = collections.defaultdict(list)
        for _, row in conv_df.iterrows():
            conv[row['Ensembl_gene_identifier']].append(conv['GeneID'])

        res = []
        for e in ens:
            res.extend(conv.get(e, []))
        # alternative: but this takes the LAST definition for any given Entrez ID, ignoring any previous duplicates
        # conv = dict(set([tuple(t) for t in conv_df.values.tolist()]))
        return res

