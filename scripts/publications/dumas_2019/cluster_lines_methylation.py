"""
Added on 4th Sept 2018

Modified version of scripts.assess_reprogramming.dna_methylation

Load DNA methylation array data and produce PCA and hierarchical clustering representations.
"""

from methylation import loader, process
from plotting import clustering, pca, common, scatter
import pandas as pd
import collections
from matplotlib import pyplot as plt
import seaborn as sns
from utils import output
from stats import transformations, basic
import os, sys
import collections
from scripts.hgic_final import consts
from scripts.publications.dumas_2019.cluster_lines_rnaseq import format_clustermap


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
    # norm_method = 'swan'
    n_hipsci = 12
    # qn_method = 'median'
    qn_method = None

    outdir = output.unique_output_dir()
    # load 12 patients iNSC, 4 iPSC
    pids = consts.PIDS

    # we'll list our samples explicitly to avoid results changing in future
    our_samples = [
        'DURA019_FB_P7',
        'DURA019_IPSC_N8C_P13',
        'DURA030_FB_P8',
        'DURA030_IPSC_N16B6_P13',
        'DURA031_FB_P7',
        'DURA031_IPSC_N44B_P10',
        'DURA017_FB_P7',
        'DURA050_IPSC_N12_P5',
        'DURA050_FB_P7',
        'DURA054_IPSC_N3C_P11',
        'DURA054_FB_P5',
    ]

    patient_obj = loader.load_by_patient(pids, norm_method=norm_method, samples=our_samples)

    nazor_ldr = loader.load_reference('GSE31848', norm_method=norm_method)
    ix = nazor_ldr.meta.index.str.contains(r'(ES__WA)|(iPS__HDF)')
    ix = ix & (~nazor_ldr.meta.index.str.contains(r'HDF51IPS7'))  # this sample is an outlier, so remove it now
    nazor_ldr.filter_samples(ix)

    # Zhou et al.: lots of samples here, but we'll only keep 2 x ESC lines
    zhou_ldr = loader.load_reference('GSE92462_450K', norm_method=norm_method)
    ix = zhou_ldr.meta.index.str.contains(r'^H[19]ES')
    zhou_ldr.filter_samples(ix)

    hip_epic_ldr = loader.hipsci(norm_method=norm_method, n_sample=n_hipsci, array_type='epic')
    ## FIXME: this is required to avoid a BUG where the meta column gets renamed to batch_1 in all other loaders
    hip_epic_ldr.meta.drop('batch', axis=1, inplace=True)
    hip_epic_ldr.batch_id = 'HipSci'

    # Weltner et al. (E-MTAB-6194)
    e6194_ldr = loader.load_reference('E-MTAB-6194', norm_method=norm_method)
    ix = ~e6194_ldr.meta.cell_line.isin([
        'NA07057',
        'HCT116',
        'HEL46.11',
        'CCD-1112Sk (CRL-2429)'
    ])
    e6194_ldr.filter_samples(ix)

    refs = collections.OrderedDict([
        ('Kim et al.', loader.load_reference('GSE38216', norm_method=norm_method, samples=['H9 ESC 1', 'H9 ESC 2'])),
        ('Zimmerlin et al.', loader.load_reference('GSE65214', norm_method=norm_method)),
        ('Nazor et al.', nazor_ldr),
        ('Encode EPIC', loader.load_reference('ENCODE_EPIC', norm_method=norm_method, samples=['GM23338', 'H7 hESC', 'Astrocyte'])),
        ('Encode 450k', loader.load_reference('ENCODE_450k', norm_method=norm_method)),
        ('Zhou et al.', zhou_ldr),
        ('Hipsci EPIC', hip_epic_ldr),
        ('Weltner et al.', e6194_ldr)
    ])
    refs['Encode EPIC'].batch_id = 'ENCODE'
    refs['Encode 450k'].batch_id = 'ENCODE'

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

    ref_obj = loader.loader.MultipleBatchLoader(refs.values())

    for srch, repl in to_aggr:
        ref_obj.aggregate_by_pattern(srch, repl)

    # rename to include publication / reference
    for k, v in ref_name_map.items():
        ix = ref_obj.meta.batch == k
        if ix.sum() > 0:
            ref_obj.meta.loc[ix, 'batch'] = v
    ref_obj.rename_with_attributes(existing_attr='batch')

    # restrict references to relevant samples
    ix = ref_obj.meta.type.isin(['iPSC', 'ESC', 'PSC', 'FB'])
    ref_obj.filter_samples(ix)

    obj = loader.loader.MultipleBatchLoader([patient_obj, ref_obj])

    # tidy up batch IDs
    obj.meta.loc[obj.meta.batch.isnull(), 'batch'] = obj.meta.loc[obj.meta.batch.isnull(), 'batch_1']
    obj.meta.batch = obj.meta.batch.str.replace('2016-12-19_ucl_genomics', '2016-12-19')

    # the only batch names without letters are ours
    obj.meta.loc[~obj.meta.batch.str.contains(r'[A-Z]'), 'batch'] = 'This study'

    # manually rename our samples for nicer plotting
    obj.meta.index = obj.meta.index.str.replace(r'DURA(?P<n>[0-9]{3})_(?P<typ>[^_]*).*', r'\g<typ>\g<n>')
    obj.data.columns = obj.meta.index

    bdat = obj.data
    mdat = process.m_from_beta(bdat)

    if qn_method is not None:
        mdat = transformations.quantile_normalisation(mdat, method=qn_method)

    # PCA plot (by batch and cell type)
    colour_subgroups = obj.meta.batch
    c_sub_sorted = sorted(colour_subgroups.unique(), key=lambda x: 'A' if x == 'This study' else x)

    cmap = collections.OrderedDict(zip(
        c_sub_sorted,
        common.get_best_cmap(len(c_sub_sorted)),
    ))

    m_subgroups = obj.meta.type
    subgroups_sorted = sorted(m_subgroups.unique(), key=lambda x: x[0] if x[0] != 'i' else x[1])
    mmap = pd.Series(
        common.FILLED_MARKERS[len(m_subgroups.unique())],
        index=subgroups_sorted
    )

    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot(111)
    p, ax = plot_pca(
        mdat,
        colour_subgroups,
        colour_map=cmap,
        marker_subgroups=m_subgroups,
        marker_map=mmap,
        ax=ax
    )
    # increase legend font size
    ax.figure.subplots_adjust(left=0.12, right=0.75, top=0.98)
    ax.figure.savefig(os.path.join(outdir, "pca_ipsc_esc_fb.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "pca_ipsc_esc_fb.tiff"), dpi=200)

    # plot dendrograms

    cell_line_colours = {
        'FB': '#fff89e',  # yellow
        'GBM (this study)': '#e6e6e6',  # light gray
        'GBM': '#4d4d4d',  # dark grey
        'ESC': '#ff7777',  # light red
        'iPSC': '#990000',  # dark red
        'iPSC (this study)': '#fdc086',  # orange
        'NSC': '#006600',  # dark green
        'iNSC (this study)': '#7fc97f',  # green
    }

    row_colours_all = pd.DataFrame('white', index=mdat.columns, columns=['Cell type', 'Study'])
    row_colours_all.loc[obj.meta.type == 'FB', 'Cell type'] = cell_line_colours['FB']
    row_colours_all.loc[obj.meta.type == 'iNSC', 'Cell type'] = cell_line_colours['iNSC (this study)']
    row_colours_all.loc[obj.meta.type == 'NSC', 'Cell type'] = cell_line_colours['NSC']
    row_colours_all.loc[obj.meta.type == 'iPSC', 'Cell type'] = cell_line_colours['iPSC (this study)']
    row_colours_all.loc[obj.meta.type == 'ESC', 'Cell type'] = cell_line_colours['ESC']
    row_colours_all.loc[obj.meta.type == 'PSC', 'Cell type'] = cell_line_colours['ESC']

    # specific cases that are not iPSC
    row_colours_all.loc[row_colours_all.index.str.contains('GM23338'), 'Cell type'] = cell_line_colours['ESC']
    row_colours_all.loc[row_colours_all.index.str.contains(r'iPS_', 'Cell type')] = cell_line_colours['iPSC']
    row_colours_all.loc[row_colours_all.index.str.contains('HPSI'), 'Cell type'] = cell_line_colours['iPSC']
    row_colours_all.loc[row_colours_all.index.str.contains('HEL1'), 'Cell type'] = cell_line_colours['iPSC']

    study_cmap = common.get_best_cmap(len(obj.meta.batch.unique()))
    study_colours = dict(zip(
        obj.meta.batch.unique(),
        study_cmap
    ))
    row_colours_all.loc[obj.meta.index, 'Study'] = obj.meta.batch.apply(study_colours.get)

    # hierarchical clustering by M value, keeping only variable probes
    clust_n_ftr = 20000
    n_probe_to_show = 2000

    the_dat = mdat

    plt_dict = hc_plot_dendrogram_vary_n_gene(the_dat, row_colours_all, n_ftr=[clust_n_ftr])
    for ng, x in plt_dict.items():
        fname = "dendrogram_M_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)

    # heatmap: use clustering from n=20000 probes (M vals), but show fewer probes values

    # pick these using MAD
    this_mad = transformations.median_absolute_deviation(mdat).sort_values(ascending=False)
    this_dat = mdat.loc[this_mad.index[:n_probe_to_show]]

    leg_entry = {
        'class': 'patch',
        'edgecolor': 'k',
        'linewidth': 1.,
    }

    lkg = plt_dict[clust_n_ftr]['linkage']

    leg_dict = collections.OrderedDict()
    for k, cmap in zip(['Cell type', 'Study'], [cell_line_colours, study_colours]):
        this_leg_dict = collections.OrderedDict()
        for j in sorted(cmap):
            if cmap[j] in row_colours_all[k].values:
                this_leg_dict[j] = dict(leg_entry)
                this_leg_dict[j].update({'facecolor': cmap[j]})
        leg_dict[k] = this_leg_dict

    gc = clustering.plot_clustermap(
        this_dat,
        cmap='RdYlBu_r',
        metric='correlation',
        col_colors=row_colours_all,
        col_linkage=lkg,
        vmin=-10,
        vmax=10
    )
    common.add_custom_legend(gc.ax_heatmap, leg_dict, loc_outside=True, fontsize=14)
    gc.fig.set_size_inches((9.4, 8.))
    # reduce cbar ticks
    cbar = gc.ax_heatmap.collections[0].colorbar
    cbar.set_ticks([-8, 0, 8])

    plt.setp(gc.cax.yaxis.get_ticklabels(), fontsize=12)
    plt.setp(gc.ax_col_colors.yaxis.get_ticklabels(), fontsize=12)
    gc.cax.yaxis.set_ticks_position('left')
    gc.cax.set_title('M value')

    gc.gs.set_height_ratios([0.16, 0.04, 0.08, 0.8])
    gc.gs.update(right=0.76, left=0.04, bottom=0.3, wspace=0., top=0.96)

    gc.savefig(os.path.join(outdir, "clustermap_ipsc_esc_fb.png"), dpi=200)
    gc.savefig(os.path.join(outdir, "clustermap_ipsc_esc_fb.tiff"), dpi=200)
