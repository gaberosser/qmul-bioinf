from methylation import loader, process
from plotting import clustering, pca, common
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils.output import unique_output_dir
from stats import transformations
import os
from settings import LOCAL_DATA_DIR


def plot_pca(
        dat,
        colour_subgroups,
        p=None,
        components=(0, 1),
        marker_subgroups=None,
        ax=None,
        colour_map=None,
        marker_map=None,
        marker_size=40,
        legend='outside',
        default_colour='gray',
        default_marker='o',
        edge_colour='k',
):
    if p is None:
        p = pca.PCA()
        pca_data = p.fit_transform(dat.transpose())
    else:
        pca_data = p.transform(dat.transpose())
    variance_explained = p.explained_variance_ratio_ * 100.

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')

    cidx, clabels = colour_subgroups.factorize()

    m_label_same_as_c = False
    if marker_subgroups is None:
        m_label_same_as_c = True
        marker_subgroups = colour_subgroups

    if colour_map is None:
        cmap = pca.get_cmap(len(clabels))
        colour_map = dict([
            (k, cmap(i)) for i, k in enumerate(clabels)
        ])

    midx, mlabels = marker_subgroups.factorize()

    if marker_map is None:
        marker_map = dict([(k, 'o') for k in clabels])

    for i, l in enumerate(clabels):
        for i2, l2 in enumerate(mlabels):
            c = colour_map.get(l, default_colour)
            m = marker_map.get(l2, default_marker)
            j = (cidx == i) & (midx == i2)
            if j.sum() != 0:
                if m_label_same_as_c:
                    lbl = l
                else:
                    lbl = None
                ax.scatter(
                    pca_data[j, components[0]],
                    pca_data[j, components[1]],
                    c=c,
                    s=marker_size,
                    label=lbl,
                    marker=m,
                    edgecolor=edge_colour
                )

    # set legend manually if it requires two groups
    if not m_label_same_as_c:
        pass
        #TODO

    ax.set_xlabel("PCA component %s (%.1f%%)" % (components[0] + 1, variance_explained[components[0]]))
    ax.set_ylabel("PCA component %s (%.1f%%)" % (components[1] + 1, variance_explained[components[1]]))

    if (legend is not None) and (legend != False):
        if legend == 'outside':
            common.legend_outside_axes(ax)
        elif isinstance(legend, str):
            ax.legend(loc=legend)
        else:
            ax.legend()

    return p, ax


def hc_plot_dendrogram_vary_n_gene(
        data,
        row_colours,
        mad = None,
        n_ftr=(1000, 2000, 3000, 5000, 10000),
        metric='correlation'
):
    """
    For each value in n_gene_arr, plot a dendrogram showing the result of hierarchical clustering of the data using
    that many genes (selected in descending MAD order)
    :param data: Cols are samples, rows are genes (or similar)
    :param row_colours: As passed to dendrogram routine
    :param n_gene_arr: The values to test
    :return:
    """
    if mad is None:
        mad = transformations.median_absolute_deviation(data).sort_values(ascending=False)
    fig_dict = {}
    for ng in n_ftr:
        the_dat = data.loc[mad.index[:ng]]
        d = clustering.dendrogram_with_colours(
            the_dat,
            row_colours,
            fig_kws={'figsize': (5.5, 10)},
            vertical=False,
            metric=metric
        )
        fig_dict[ng] = d
    return fig_dict



if __name__ == "__main__":
    norm_method = 'pbc'
    outdir = unique_output_dir("methylation_insc_characterisation", reuse_empty=True)
    # load 12 patients iNSC, 4 iPSC
    pids = ['017', '018', '019', '030', '031', '026', '044', '049', '050', '052', '054', '061']

    patient_obj = loader.load_by_patient(pids, norm_method=norm_method)

    refs = [
        ('Kim et al.', loader.gse38216(norm_method=norm_method, samples=['H9 ESC 1', 'H9 ESC 2', 'H9 NPC 1', 'H9 NPC 2'])),
        ('Morey et al.', loader.gse67283(norm_method=norm_method, samples=['NPC_wt'])),
        ('Zimmerlin et al.', loader.gse65214(norm_method=norm_method)),
        ('Encode EPIC', loader.encode_epic(norm_method=norm_method)),
        ('Encode 450k', loader.encode_450k(norm_method=norm_method)),
    ]

    to_aggr = [
        (r'H9 ESC [12]', 'H9 ESC'),
        (r'H9 NPC [12]', 'H9 NPC'),
    ]

    # TODO make this nicer
    r = refs[0][1]
    for srch, repl in to_aggr:
        idx = r.data.columns.str.contains(srch)
        new_col = r.data.loc[:, idx].mean(axis=1)
        r.data = r.data.loc[:, ~idx]
        r.data.insert(r.data.shape[1], repl, new_col)
        new_meta_row = r.meta.loc[idx].iloc[0]
        new_meta_row.name = repl
        r.meta = r.meta.loc[~idx]
        r.meta = r.meta.append(new_meta_row)

    # rename based on study
    for nm, t in refs:
        new_idx = ["%s (%s)" % (i, nm) for i in t.meta.index]
        t.meta.index = new_idx
        t.data.columns = new_idx

    ref_obj = loader.loader.MultipleBatchLoader([t[1] for t in refs])

    obj = loader.loader.MultipleBatchLoader([patient_obj, ref_obj])

    bdat = obj.data
    mdat = process.m_from_beta(bdat)

    # PCA plots
    # 1. By batch and cell type
    subgroups = obj.batch_id
    cmap = pd.Series(
        common.COLOUR_BREWERS[len(subgroups.unique())],
        index=subgroups.unique()
    )
    m_subgroups = obj.meta.type
    mmap = pd.Series(
        common.FILLED_MARKERS[len(m_subgroups.unique())],
        index=m_subgroups.unique()
    )

    plot_pca(
        mdat,
        subgroups,
        marker_subgroups=m_subgroups,
        marker_map=mmap
    )


    bmad = transformations.median_absolute_deviation(bdat).sort_values(ascending=False)
    mmad = transformations.median_absolute_deviation(mdat).sort_values(ascending=False)

    row_colours_all = pd.DataFrame('gray', index=mdat.columns, columns=[''])

    row_colours_all.loc[row_colours_all.index.str.contains(r'GBM')] = '#fff89e'
    row_colours_all.loc[row_colours_all.index.str.contains(r'NS27Z')] = '#fff89e'

    row_colours_all.loc[row_colours_all.index.str.contains('H9')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('H7')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('H1')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('ESO3')] = '#ff7777'

    row_colours_all.loc[row_colours_all.index.str.contains(r'NSC')] = 'blue'
    row_colours_all.loc[row_colours_all.index.str.contains(r'NPC')] = 'blue'

    row_colours_all.loc[row_colours_all.index.str.contains('neuron')] = '#ccebc5'
    row_colours_all.loc[row_colours_all.index.str.contains(r'[Aa]strocyte')] = '#e78ac3'

    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_NSC')] = '#7fc97f'  # green
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_IPSC')] = '#fdc086'  # orange
    row_colours_all.loc[row_colours_all.index.str.contains(r'GIBCO')] = '#96daff'  # light blue

    n_ftr = [1000, 2000, 3000, 5000, 10000, 20000, 50000, 100000, 1000000]
    plt_dict = hc_plot_dendrogram_vary_n_gene(mdat, row_colours_all, mad=mmad, n_ftr=n_ftr)

    for ng, x in plt_dict.items():
        if ng > bdat.shape[0]:
            fname = "dendrogram_M_corr_all.{ext}"
        else:
            fname = "dendrogram_M_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)

    plt_dict = hc_plot_dendrogram_vary_n_gene(bdat, row_colours_all, mad=bmad, n_ftr=n_ftr)
    for ng, x in plt_dict.items():
        if ng > bdat.shape[0]:
            fname = "dendrogram_B_corr_all.{ext}"
        else:
            fname = "dendrogram_B_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)

    ## no GBM, no FB, no astrocytes, various numbers of probes, beta values

    idx = (
        (~mdat.columns.str.contains('GBM')) &
        (~mdat.columns.str.contains('Astrocyte')) &
        (~mdat.columns.str.contains('_FB')) &
        (~mdat.columns.str.contains('DURA030_NSC_N9_P2'))
    )

    bdat_nogbm = bdat.loc[:, idx]

    bmad_nogbm = transformations.median_absolute_deviation(bdat_nogbm)
    row_colours_nogbm = row_colours_all.loc[bdat_nogbm.columns]

    plt_dict = hc_plot_dendrogram_vary_n_gene(bdat_nogbm, row_colours_nogbm, mad=bmad_nogbm, n_ftr=n_ftr)
    for ng, x in plt_dict.items():
        if ng > bdat.shape[0]:
            fname = "nogbm_dendrogram_B_corr_all.{ext}"
        else:
            fname = "nogbm_dendrogram_B_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)

    plt.close('all')

    # heatmap for 3000 probes
    ng = 3000
    cm = clustering.plot_clustermap(
        bdat_nogbm.loc[bmad_nogbm.index[:ng]],
        # cmap='RdBu_r',
        cmap='RdYlBu_r',
        metric='correlation',
        col_colors=row_colours_nogbm
    )
    cm.savefig(os.path.join(outdir, "nogbm_clustermap_B_corr_top%d_by_mad.png" % ng), dpi=200)
    cm.savefig(os.path.join(outdir, "nogbm_clustermap_B_corr_top%d_by_mad.tiff" % ng), dpi=200)

    ## our samples only, all probes, beta values

    idx = (
        (bdat.columns.str.contains('GBM')) |
        (bdat.columns.str.contains('DURA'))
    )
    bdat_ours = bdat.loc[:, idx]
    bmad_ours = transformations.median_absolute_deviation(bdat_ours)

    row_colours_ours = row_colours_all.loc[bdat_ours.columns]
    x = clustering.dendrogram_with_colours(
                                    bdat_ours,
                                    row_colours_ours,
                                    fig_kws={'figsize': (5.5, 10)},
                                    vertical=False,
                                    metric='correlation'
    )
    fname = "ours_dendrogram_B_corr_all.{ext}"
    x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
    x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
