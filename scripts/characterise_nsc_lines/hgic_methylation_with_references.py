from methylation import loader, process
from plotting import clustering, pca, common, scatter
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils.output import unique_output_dir
from stats import transformations, basic
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
        **kwargs
):
    if p is None:
        p = pca.PCA()
        pca_data = p.fit_transform(dat.transpose())
    else:
        pca_data = p.transform(dat.transpose())
    variance_explained = p.explained_variance_ratio_ * 100.

    ax = scatter.scatter_with_colour_and_markers(
        pca_data[:, components],
        colour_subgroups=colour_subgroups,
        colour_map=colour_map,
        marker_subgroups=marker_subgroups,
        marker_map=marker_map,
        ax=ax,
        **kwargs
    )

    ax.set_xlabel("PCA component %s (%.1f%%)" % (components[0] + 1, variance_explained[components[0]]))
    ax.set_ylabel("PCA component %s (%.1f%%)" % (components[1] + 1, variance_explained[components[1]]))

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
    norm_method = 'bmiq'

    outdir = unique_output_dir("methylation_insc_characterisation", reuse_empty=True)
    # load 12 patients iNSC, 4 iPSC
    pids = ['017', '018', '019', '030', '031', '026', '044', '049', '050', '052', '054', '061']

    patient_obj = loader.load_by_patient(pids, norm_method=norm_method)
    ix = patient_obj.meta.type != 'iOPC'
    patient_obj.meta = patient_obj.meta.loc[ix]
    patient_obj.batch_id = patient_obj.batch_id[ix]
    patient_obj.data = patient_obj.data.loc[:, patient_obj.meta.index]

    nazor_ldr = loader.gse31848(norm_method=norm_method)
    ix = nazor_ldr.meta.index.str.contains(r'(ES__WA)|(iPS__HDF)')
    ix = ix & (~nazor_ldr.meta.index.str.contains(r'HDF51IPS7'))  # this sample is an outlier, so remove it now
    nazor_ldr.filter_samples(ix)

    # Zhou et al.: lots of samples here, but we'll only keep 2 x ESC lines
    zhou_ldr = loader.gse92462_450k(norm_method=norm_method)
    ix = zhou_ldr.meta.index.str.contains(r'^H[19]ES')
    zhou_ldr.filter_samples(ix)


    hip_epic_ldr = loader.hipsci(norm_method=norm_method, n_sample=12, array_type='epic')
    ## FIXME: this is required to avoid a BUG where the meta column gets renamed to batch_1 in all other loaders
    hip_epic_ldr.meta.drop('batch', axis=1, inplace=True)
    # hip_450k_meta, hip_450k_data = loader.hipsci(norm_method=norm_method, n_sample=30, array_type='450k')

    # Kurscheid et al.: only interested in the GSC spheroids here
    # kurs_ldr = loader.gse60274(norm_method=norm_method, samples=['LN-2207GS', 'LN-2540GS', 'LN-2669GS', 'LN-2683GS'])

    # Weltner et al. (E-MTAB-6194)
    e6194_ldr = loader.e_mtab_6194(norm_method=norm_method)
    ix = ~e6194_ldr.meta.cell_line.isin([
        'NA07057',
        'HCT116',
        'HEL46.11',
        'CCD-1112Sk (CRL-2429)'
    ])
    e6194_ldr.filter_samples(ix)

    refs = [
        ('Kim et al.', loader.gse38216(norm_method=norm_method, samples=['H9 ESC 1', 'H9 ESC 2', 'H9 NPC 1', 'H9 NPC 2'])),
        ('Morey et al.', loader.gse67283(norm_method=norm_method, samples=['NPC_wt'])),
        ('Zimmerlin et al.', loader.gse65214(norm_method=norm_method)),
        ('Nazor et al.', nazor_ldr),
        ('Encode EPIC', loader.encode_epic(norm_method=norm_method, samples=['GM23338', 'H9 NPC', 'H7 hESC', 'Astrocyte'])),
        ('Encode 450k', loader.encode_450k(norm_method=norm_method)),
        ('Zhou et al.', zhou_ldr),
        ('Hipsci EPIC', hip_epic_ldr),
        ('Weltner et al.', e6194_ldr)
        # ('Kurscheid et al.', kurs_ldr)
    ]

    ref_name_map = {
        'GSE38216': 'Kim et al.',
        'GSE67283': 'Morey et al.',
        'GSE65214': 'Zimmerlin et al.',
        'GSE31848': 'Nazor et al.',
        'GSE92462_450K': 'Zhou et al.',
        'GSE60274': 'Kurscheid et al.',
        'E-MTAB-6194': 'Weltner et al.'
    }

    to_aggr = [
        (r'H9 ESC [12]', 'H9 ESC'),
        (r'H9 NPC [12]', 'H9 NPC'),
        ('ES__WA09_', 'H9 ESC'),
        ('ES__WA07_', 'H7 ESC'),
        ('22_H9_', 'H9 ESC'),
        ('HEL139.2_', 'HEL139.2'),
        ('HEL139.5_', 'HEL139.5'),
        ('HEL139.8_', 'HEL139.8'),
        ('HEL140.1_', 'HEL140.1'),
        ('HEL141.1_', 'HEL141.1'),

    ]
    for i in range(1, 15):
        to_aggr.append(
            (r'iPS__HDF51IPS%d_' % i, 'iPS_HDF51_%d' % i)
        )

    ref_obj = loader.loader.MultipleBatchLoader([t[1] for t in refs])

    for srch, repl in to_aggr:
        idx = ref_obj.meta.index.str.contains(srch)
        if idx.sum() == 0:
            print "Warning: search string %s matches no samples." % srch
            continue

        # define new data and meta entries
        new_col = ref_obj.data.loc[:, idx].mean(axis=1)
        new_meta_row = ref_obj.meta.loc[idx].iloc[0]
        new_meta_row.name = repl

        # filter out old entries
        ref_obj.filter_samples(~idx)

        # add new entries: data, meta and batch_id
        ref_obj.data.insert(ref_obj.data.shape[1], repl, new_col, allow_duplicates=True)
        ref_obj.meta = ref_obj.meta.append(new_meta_row)
        ref_obj.batch_id = ref_obj.batch_id.append(pd.Series({repl: new_meta_row.batch}))

    # rename to include publication / reference
    new_index = np.array(ref_obj.meta.index.tolist()).astype(object)  # need to cast this or numpy truncates it later
    for b in ref_obj.meta.batch.unique():
        if b in ref_name_map:
            suff = ref_name_map[b]
        else:
            suff = b
        ix = ref_obj.meta.batch == b
        old_names = ref_obj.meta.index[ix]
        new_names = ["%s (%s)" % (t, suff) for t in old_names]
        new_index[ix] = new_names

    ref_obj.meta.index = new_index
    ref_obj.data.columns = new_index
    ref_obj.batch_id.index = new_index

    obj = loader.loader.MultipleBatchLoader([patient_obj, ref_obj])

    bdat = obj.data
    mdat = process.m_from_beta(bdat)

    # PCA plot (by batch and cell type)

    colour_subgroups = obj.batch_id

    m_subgroups = obj.meta.type
    mmap = pd.Series(
        common.FILLED_MARKERS[len(m_subgroups.unique())],
        index=m_subgroups.unique()
    )

    # shorten batch name
    colour_subgroups = colour_subgroups.replace('2016-12-19_ucl_genomics', '2016-12-19')

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    p, ax = plot_pca(
        mdat,
        colour_subgroups,
        marker_subgroups=m_subgroups,
        marker_map=mmap,
        ax=ax
    )
    ax.figure.subplots_adjust(left=0.1, right=0.8)
    ax.figure.savefig(os.path.join(outdir, "pca_plot_batch_cell_type_all.png"), dpi=200)

    # ECDF plot (separate for cell types) to show batch effects

    cell_types = ['GBM', 'iNSC']
    cmap = common.get_best_cmap(len(colour_subgroups.unique()))
    colour_map = dict(zip(colour_subgroups.unique(), cmap))

    xi = np.linspace(-8, 8, 200)
    for ct in cell_types:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        this_idx = obj.meta.type == ct
        this_dat = mdat.loc[:, this_idx]
        this_batch = colour_subgroups[this_idx]
        labels_included = dict([(k, False) for k in colour_subgroups.unique()])

        for i in range(this_dat.shape[1]):
            func = basic.ecdf_func(this_dat.iloc[:, i])
            yi = func(xi)
            csg = this_batch.iloc[i]
            if labels_included[csg]:
                lbl = None
            else:
                lbl = csg
                labels_included[csg] = True
            ax.plot(xi, yi, c=colour_map[csg], label=lbl)

        ax.set_xlabel('M value')
        ax.set_ylabel('ECDF')
        ax.set_title(ct)
        common.legend_outside_axes(ax)
        fig.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.99)
        fig.savefig(os.path.join(outdir, "ecdf_batches_%s.png" % ct), dpi=200)

    # normalise and repeat PCA
    mdat_qn = transformations.quantile_normalisation(mdat)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    p_qn, ax = plot_pca(
        mdat_qn,
        colour_subgroups,
        marker_subgroups=m_subgroups,
        marker_map=mmap,
        ax=ax
    )
    ax.figure.subplots_adjust(left=0.1, right=0.8)
    ax.figure.savefig(os.path.join(outdir, "pca_plot_batch_cell_type_all_qn.png"), dpi=200)

    # now try without GBM

    # (a) no QN
    ix = ~m_subgroups.str.contains('GBM')
    mdat_nogbm = mdat.loc[:, ix]
    m_subgroups_nogbm = m_subgroups.loc[ix]
    colour_subgroups_nogbm = colour_subgroups.loc[ix]
    p_nogbm, ax = plot_pca(
        mdat_nogbm,
        colour_subgroups_nogbm,
        marker_subgroups=m_subgroups_nogbm,
        marker_map=mmap
    )
    ax.figure.subplots_adjust(left=0.1, right=0.8)
    ax.figure.savefig(os.path.join(outdir, "pca_plot_batch_cell_type_nogbm.png"), dpi=200)

    mdat_qn_nogbm = mdat_qn.loc[:, ix]
    p_qn_nogbm, ax = plot_pca(
        mdat_qn_nogbm,
        colour_subgroups_nogbm,
        marker_subgroups=m_subgroups_nogbm,
        marker_map=mmap
    )
    ax.figure.subplots_adjust(left=0.1, right=0.8)
    ax.figure.savefig(os.path.join(outdir, "pca_plot_batch_cell_type_nogbm_qn.png"), dpi=200)


    # plot dendrograms

    # bmad = transformations.median_absolute_deviation(bdat).sort_values(ascending=False)
    # mmad = transformations.median_absolute_deviation(mdat).sort_values(ascending=False)

    row_colours_all = pd.DataFrame('white', index=mdat.columns, columns=[''])

    row_colours_all.loc[row_colours_all.index.str.contains(r'_FB_')] = '#fff89e'  # yellow

    row_colours_all.loc[row_colours_all.index.str.contains(r'GBM')] = '#e6e6e6'  # light gray
    row_colours_all.loc[row_colours_all.index.str.contains(r'NS27Z')] = '#4d4d4d'  # dark grey

    row_colours_all.loc[row_colours_all.index.str.contains('H9')] = '#ff7777'  # light red
    row_colours_all.loc[row_colours_all.index.str.contains('H7')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('H1')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('ESO3')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('GM23338')] = '#ff7777'

    row_colours_all.loc[row_colours_all.index.str.contains(r'NSC')] = '#ccffcc'  # pale green
    row_colours_all.loc[row_colours_all.index.str.contains(r'NPC')] = '#ccffcc'
    row_colours_all.loc[row_colours_all.index.str.contains(r'GIBCO')] = '#006600'  # dark green

    row_colours_all.loc[row_colours_all.index.str.contains(r'iPS_')] = '#990000'  # dark red
    row_colours_all.loc[row_colours_all.index.str.contains('HPSI')] = '#990000'
    row_colours_all.loc[row_colours_all.index.str.contains('HEL1')] = '#990000'

    row_colours_all.loc[row_colours_all.index.str.contains('neuron')] = '#ccebc5'
    row_colours_all.loc[row_colours_all.index.str.contains(r'[Aa]strocyte')] = '#e78ac3'  # pink

    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_NSC')] = '#7fc97f'  # green
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_IPSC')] = '#fdc086'  # orange


    # by M value, keeping only variable probes

    # min_delta_m = 14
    # m_range = mdat.max(axis=1) - mdat.min(axis=1)
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.hist(m_range.loc[m_range < 20], 300)
    # ax.set_xlabel('M range')
    # ax.axvline(min_delta_m, c='k', ls='--', lw=2.)
    # fig.tight_layout()
    # fig.savefig(os.path.join(outdir, "m_range_hist.png"), dpi=200)
    #
    # ix = m_range > min_delta_m
    # the_dat = mdat.loc[ix]

    n_ftr = [1000, 2000, 3000, 5000, 10000, 20000, int(1e7)]

    the_dat = mdat

    plt_dict = hc_plot_dendrogram_vary_n_gene(the_dat, row_colours_all, n_ftr=n_ftr)
    for ng, x in plt_dict.items():
        if ng > the_dat.shape[0]:
            fname = "dendrogram_M_corr_all.{ext}"
        else:
            fname = "dendrogram_M_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)
    plt.close('all')

    # heatmap: use clustering from n=20000 probes (M vals), but only show top 500 most variable between clusters
    clust_n_ftr = 20000
    n_probe_to_show = 500
    lkg = plt_dict[clust_n_ftr]['linkage']
    this_mad = transformations.median_absolute_deviation(mdat).sort_values(ascending=False)
    this_dat = mdat.loc[this_mad.index[:n_probe_to_show]]

    # heatmap for 3000 probes
    cm = clustering.plot_clustermap(
        this_dat,
        cmap='RdYlBu_r',
        metric='correlation',
        col_colors=row_colours_all,
        col_linkage=lkg,
        vmin=-10,
        vmax=10,
        figsize=(11.8, 10.)
    )
    cm.gs.update(bottom=0.25, right=0.99)

    cm.savefig(os.path.join(outdir, "clustermap_M_corr_linkage%d_heatmap%d.png" % (clust_n_ftr, n_probe_to_show)), dpi=200)
    cm.savefig(os.path.join(outdir, "clustermap_M_corr_linkage%d_heatmap%d.tiff" % (clust_n_ftr, n_probe_to_show)), dpi=200)

    # by beta value, keeping only variable probes
    # min_delta_beta = 0.6
    # b_range = bdat.max(axis=1) - bdat.min(axis=1)
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.hist(b_range, 300)
    # ax.set_xlabel('Beta range')
    # ax.axvline(min_delta_beta, c='k', ls='--', lw=2.)
    # fig.tight_layout()
    # fig.savefig(os.path.join(outdir, "beta_range_hist.png"), dpi=200)
    #
    # ix = b_range > min_delta_beta
    # the_dat = bdat.loc[ix]

    the_dat = bdat

    plt_dict = hc_plot_dendrogram_vary_n_gene(the_dat, row_colours_all, n_ftr=n_ftr)

    for ng, x in plt_dict.items():
        if ng > bdat.shape[0]:
            fname = "dendrogram_B_corr_all.{ext}"
        else:
            fname = "dendrogram_B_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)
    plt.close('all')

    ## no GBM, no FB, no astrocytes, various numbers of probes, beta values

    idx = (
        (~mdat.columns.str.contains('GBM')) &
        (~mdat.columns.str.contains('GSC')) &
        (~mdat.columns.str.contains('LN-')) &
        (~mdat.columns.str.contains('Astrocyte')) &
        (~mdat.columns.str.contains('_FB'))
        # I had been removing one of the iNSC030 lines, having branded it an outlier
        # can't see the evidence for this any more?
        # (~mdat.columns.str.contains('DURA030_NSC_N9_P2'))
    )

    bdat_nogbm = bdat.loc[:, idx]

    # bmad_nogbm = transformations.median_absolute_deviation(bdat_nogbm)
    row_colours_nogbm = row_colours_all.loc[bdat_nogbm.columns]

    plt_dict = hc_plot_dendrogram_vary_n_gene(bdat_nogbm, row_colours_nogbm, n_ftr=n_ftr)
    for ng, x in plt_dict.items():
        if ng > bdat.shape[0]:
            fname = "nogbm_dendrogram_B_corr_all.{ext}"
        else:
            fname = "nogbm_dendrogram_B_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)
    plt.close('all')

    mdat_nogbm = mdat.loc[:, idx]
    plt_dict = hc_plot_dendrogram_vary_n_gene(mdat_nogbm, row_colours_nogbm, n_ftr=n_ftr)
    for ng, x in plt_dict.items():
        if ng > bdat.shape[0]:
            fname = "nogbm_dendrogram_M_corr_all.{ext}"
        else:
            fname = "nogbm_dendrogram_M_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)
    plt.close('all')

    # heatmap: use clustering from n=20000 probes (M vals), but only show top 500 most variable between clusters

    clust_n_ftr = 20000
    n_probe_to_show = 500
    lkg = plt_dict[clust_n_ftr]['linkage']
    this_mad = transformations.median_absolute_deviation(mdat_nogbm).sort_values(ascending=False)
    this_dat = mdat_nogbm.loc[this_mad.index[:n_probe_to_show]]

    # heatmap for 3000 probes
    cm = clustering.plot_clustermap(
        this_dat,
        cmap='RdYlBu_r',
        metric='correlation',
        col_colors=row_colours_nogbm,
        col_linkage=lkg,
        vmin=-10,
        vmax=10
    )
    cm.gs.update(bottom=0.25, right=0.99)

    cm.savefig(os.path.join(outdir, "nogbm_clustermap_M_corr_linkage%d_heatmap%d.png" % (clust_n_ftr, n_probe_to_show)), dpi=200)
    cm.savefig(os.path.join(outdir, "nogbm_clustermap_M_corr_linkage%d_heatmap%d.tiff" % (clust_n_ftr, n_probe_to_show)), dpi=200)

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
