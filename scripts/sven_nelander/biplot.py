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
from utils import output
import references

import plotly.plotly as py
from plotly import graph_objs as go


def generate_plotly_plot(
        res,
        filename,
        feature_size_scaling,
        feature_text,
        sample_text,
        sample_colours,
        sample_markers,
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
    x, y = res['feature_data']

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

    sample_trace = go.Scatter(
        x=res['sample_data'][0],
        y=res['sample_data'][1],
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

    pcs = res['components']
    ev_pct = res['explained_variance_ratio'] * 100.

    layout = go.Layout(
        showlegend=False,
        hovermode='closest',
        width=800,
        height=600,
        xaxis=go.layout.XAxis(title='PC %d (%.2f %%)' % (pcs[0] + 1, ev_pct[pcs[0]])),
        yaxis=go.layout.YAxis(title='PC %d (%.2f %%)' % (pcs[1] + 1, ev_pct[pcs[1]])),
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

    # arrowprops = {
    #     'color': 'black',
    #     'linewidth': 1.,
    #     'arrowstyle': '-',
    # }

    if annotate_features_radius is not None:
        # annotate most influential genes
        selected = pca.highlight_biplot_features(feat_x, feat_y, annotate_features_radius, ax)
        genes_selected = dat.index[selected]
        symbols_selected = references.ensembl_to_gene_symbol(genes_selected)

        # add gene symbol annotations
        text_handles = []
        for ix, gs in zip(np.where(selected)[0], symbols_selected):
            if not pd.isnull(gs):
                text_handles.append(ax.text(feat_x[ix], feat_y[ix], gs, zorder=10))
        # rearrange them to avoid overlaps
        if adjust_annotation:
            # adjust_text(
            #     text_handles,
            #     arrowprops=arrowprops,
            #     text_from_points=False,
            #     draggable=False,
            #     ax=ax
            # )
            adjuster.adjust_text_radial_plus_repulsion(text_handles, **adjust_annotation_kwargs)

    return fig, ax, res




if __name__ == "__main__":
    """
    Idea here: recreate the analysis Sven carried out, generating biplots for the RNA-Seq data.
    We can then extend this idea to methylation (?)

    Absolutely amazing resource here:
    https://www.fbbva.es/microsite/multivariate-statistics/biplots.html

    Code snippet inspiration here:
    https://stackoverflow.com/questions/39216897/plot-pca-loadings-and-loading-in-biplot-in-sklearn-like-rs-autoplot
    """
    outdir = output.unique_output_dir()
    publish_plotly = False

    eps = 1.  # offset applied during log transform

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

    dims = (0, 1)
    selection_radius = 0.6

    # recreate Sven's original (unweighted) plot
    fig, ax, _ = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        annotate_features_radius=0.4,
        include_weighting=False,
        scale=10.
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated_unweighted.png" % dims), dpi=200)


    fig, ax, res = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        scale=0.05
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d.png" % dims), dpi=200)

    fig, ax, res = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        annotate_features_radius=selection_radius,
        scale=0.05
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated.png" % dims), dpi=200)

    if publish_plotly:
        selection_radius = 0.6
        size_scaling = [
            [0.1, 0.6],
            [2., 10.]
        ]
        feature_text = dat_with_gs['Gene Symbol']
        sample_text = dat.columns
        rad = (np.array(zip(*res['feature_data'])) ** 2).sum(axis=1) ** .5
        to_annotate = rad > selection_radius
        p1 = generate_plotly_plot(
            res,
            filename="pca_biplot_dims_%d-%d" % dims,
            feature_size_scaling=size_scaling,
            feature_text=feature_text,
            sample_text=sample_text,
            sample_colours=sample_colours,
            sample_markers=sample_markers,
            feature_text_mask=~to_annotate
        )

    dims = (1, 2)
    selection_radius = 0.3

    fig, ax, res = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        scale=0.05
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d.png" % dims), dpi=200)

    fig, ax, res = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        annotate_features_radius=selection_radius,
        scale=0.05
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated.png" % dims), dpi=200)

    if publish_plotly:
        selection_radius = 0.3
        size_scaling = [
            [0.1, 0.3],
            [2., 10.]
        ]
        rad = (np.array(zip(*res['feature_data'])) ** 2).sum(axis=1) ** .5
        to_annotate = rad > selection_radius
        p2 = generate_plotly_plot(
            res,
            filename="pca_biplot_dims_%d-%d" % dims,
            feature_size_scaling=size_scaling,
            feature_text=feature_text,
            sample_text=sample_text,
            sample_colours=sample_colours,
            sample_markers=sample_markers,
            feature_text_mask=~to_annotate
        )

    dims = (2, 3)
    selection_radius = 0.25

    fig, ax, res = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        scale=0.05
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d.png" % dims), dpi=200)

    fig, ax, res = plot_biplot(
        dat,
        obj.meta,
        dims,
        scatter_colours,
        scatter_markers,
        annotate_features_radius=selection_radius,
        scale=0.05
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated.png" % dims), dpi=200)

    if publish_plotly:
        selection_radius = 0.25
        size_scaling = [
            [0.1, 0.3],
            [2., 10.]
        ]
        rad = (np.array(zip(*res['feature_data'])) ** 2).sum(axis=1) ** .5
        to_annotate = rad > selection_radius
        p3 = generate_plotly_plot(
            res,
            filename="pca_biplot_dims_%d-%d" % dims,
            feature_size_scaling=size_scaling,
            feature_text=feature_text,
            sample_text=sample_text,
            sample_colours=sample_colours,
            sample_markers=sample_markers,
            feature_text_mask=~to_annotate
        )

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



    # TODO
#     tt = ax.text(
#         for_scatter[0].loc[e],
#         for_scatter[1].loc[e],
#         txt,
#         color=c,
#         weight=weight
#     )
#     texts.append(tt)
# adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'), ax=ax)