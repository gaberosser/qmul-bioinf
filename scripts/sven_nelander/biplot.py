import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import pandas as pd
import itertools
import collections
import os

from rnaseq import loader, filter, general
from scripts.hgic_final import consts
from plotting import common, pca, _plotly
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
        auto_open=False
):
    x, y = res['feature_data']

    msize_in = (x ** 2 + y ** 2) ** .5

    tmp = (msize_in - feature_size_scaling[0][0]) / (feature_size_scaling[0][1] - feature_size_scaling[0][0])
    tmp[tmp < 0] = 0
    tmp[tmp > 1] = 1
    msize = tmp * (feature_size_scaling[1][1] - feature_size_scaling[1][0]) + feature_size_scaling[1][0]

    feat_trace = go.Scatter(
        x=x,
        y=y,
        mode='markers',
        text=feature_text,
        marker={
            'size': msize,
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

    layout = go.Layout(showlegend=False)
    fig = go.Figure(data=[feat_trace, sample_trace], layout=layout)
    p = py.plot(fig, filename=filename, auto_open=auto_open)

    return p


def plot_biplot(
        dat,
        meta,
        dims,
        scatter_colours,
        scatter_markers,
        annotate_features_radius=None,
        **kwargs
):
    """

    :param dat:
    :param meta: pd.DataFrame, must have columns entitled `type` and `patient_id`
    :param dims:
    :param scatter_colours:
    :param scatter_markers:
    :param **kwargs: Passed to pca.biplot()
    :return:
    """
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

    common.add_custom_legend(ax, legend_dict, loc_outside=True)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)

    if annotate_features_radius is not None:
        # annotate most influential genes
        selected = pca.highlight_biplot_features(feat_x, feat_y, annotate_features_radius, ax)
        genes_selected = dat.index[selected]
        symbols_selected = references.ensembl_to_gene_symbol(genes_selected)

        # add gene symbol annotations
        for ix, gs in zip(np.where(selected)[0], symbols_selected):
            if not pd.isnull(gs):
                ax.text(feat_x[ix], feat_y[ix], gs)

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
        annotate_features_radius=0.6,
        scale=0.05
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated.png" % dims), dpi=200)

    size_scaling = [
        [0.1, 0.6],
        [2., 10.]
    ]
    feature_text = dat_with_gs['Gene Symbol']
    sample_text = dat.columns
    p1 = generate_plotly_plot(
        res,
        filename="pca_biplot_dims_%d-%d" % dims,
        feature_size_scaling=size_scaling,
        feature_text=feature_text,
        sample_text=sample_text,
        sample_colours=sample_colours,
        sample_markers=sample_markers,
    )

    dims = (1, 2)

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
        annotate_features_radius=0.3,
        scale=0.05
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated.png" % dims), dpi=200)

    size_scaling = [
        [0.1, 0.3],
        [2., 10.]
    ]
    p2 = generate_plotly_plot(
        res,
        filename="pca_biplot_dims_%d-%d" % dims,
        feature_size_scaling=size_scaling,
        feature_text=feature_text,
        sample_text=sample_text,
        sample_colours=sample_colours,
        sample_markers=sample_markers,
    )

    dims = (2, 3)

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
        annotate_features_radius=0.3,
        scale=0.05
    )
    fig.savefig(os.path.join(outdir, "pca_biplot_dims_%d-%d_annotated.png" % dims), dpi=200)

    size_scaling = [
        [0.1, 0.3],
        [2., 10.]
    ]
    p3 = generate_plotly_plot(
        res,
        filename="pca_biplot_dims_%d-%d" % dims,
        feature_size_scaling=size_scaling,
        feature_text=feature_text,
        sample_text=sample_text,
        sample_colours=sample_colours,
        sample_markers=sample_markers,
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

    ###### experimental: plotting with Plotly ######
