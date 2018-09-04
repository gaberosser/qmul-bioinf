"""
Added on 4th Sept 2018

Modified version of scripts.assess_reprogramming.dna_methylation

Load DNA methylation array data and produce PCA and hierarchical clustering representations.
"""

from methylation import loader, process
from plotting import clustering, pca, common, scatter
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils import output
from stats import transformations, basic
import os, sys
import collections
import consts


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
    n_hipsci = 12
    # qn_method = 'median'
    qn_method = None

    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    outdir = output.unique_output_dir(script_name)
    # load 12 patients iNSC, 4 iPSC
    pids = consts.PIDS

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

    hip_epic_ldr = loader.hipsci(norm_method=norm_method, n_sample=n_hipsci, array_type='epic')
    ## FIXME: this is required to avoid a BUG where the meta column gets renamed to batch_1 in all other loaders
    hip_epic_ldr.meta.drop('batch', axis=1, inplace=True)

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
        ref_obj.aggregate_by_pattern(srch, repl)

    # rename to include publication / reference
    for k, v in ref_name_map.items():
        ix = ref_obj.meta.batch == k
        if ix.sum() > 0:
            ref_obj.meta.loc[ix, 'batch'] = v
    ref_obj.rename_with_attributes(existing_attr='batch')

    obj = loader.loader.MultipleBatchLoader([patient_obj, ref_obj])

    # remove a few
    ix = obj.meta.type != 'astrocyte'
    obj.filter_samples(ix)

    ix = ~obj.meta.index.str.contains('GBM')
    obj.filter_samples(ix)

    ix = obj.meta.index != 'H9 NPC (Encode EPIC)'
    obj.filter_samples(ix)

    bdat = obj.data
    mdat = process.m_from_beta(bdat)

    if qn_method is not None:
        mdat = transformations.quantile_normalisation(mdat, method=qn_method)

    # tidy up batch IDs
    obj.meta.batch = obj.meta.batch.str.replace('2016-12-19_ucl_genomics', '2016-12-19')
    # the only batch names without letters are ours
    obj.meta.loc[~obj.meta.batch.str.contains(r'[A-Z]'), 'batch'] = 'This study'

    # PCA plot (by batch and cell type)

    colour_subgroups = obj.meta.batch

    m_subgroups = obj.meta.type
    mmap = pd.Series(
        common.FILLED_MARKERS[len(m_subgroups.unique())],
        index=m_subgroups.unique()
    )

    fig = plt.figure()
    ax = fig.add_subplot(111)
    p, ax = plot_pca(
        mdat,
        colour_subgroups,
        marker_subgroups=m_subgroups,
        marker_map=mmap,
        ax=ax
    )
    ax.figure.subplots_adjust(left=0.1, right=0.8, top=0.98)
    ax.figure.savefig(os.path.join(outdir, "pca_ipsc_esc_nsc_fb.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir, "pca_ipsc_esc_nsc_fb.tiff"), dpi=200)

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

    row_colours_all = pd.DataFrame('white', index=mdat.columns, columns=[''])

    row_colours_all.loc[row_colours_all.index.str.contains(r'_FB_')] = cell_line_colours['FB']

    row_colours_all.loc[row_colours_all.index.str.contains(r'GBM')] = cell_line_colours['GBM (this study)']
    row_colours_all.loc[row_colours_all.index.str.contains(r'NS27Z')] = cell_line_colours['GBM']

    row_colours_all.loc[row_colours_all.index.str.contains('H9')] = cell_line_colours['ESC']
    row_colours_all.loc[row_colours_all.index.str.contains('H7')] = cell_line_colours['ESC']
    row_colours_all.loc[row_colours_all.index.str.contains('H1')] = cell_line_colours['ESC']
    row_colours_all.loc[row_colours_all.index.str.contains('ESO3')] = cell_line_colours['ESC']
    row_colours_all.loc[row_colours_all.index.str.contains('GM23338')] = cell_line_colours['ESC']

    row_colours_all.loc[row_colours_all.index.str.contains(r'NSC')] = cell_line_colours['NSC']
    row_colours_all.loc[row_colours_all.index.str.contains(r'NPC')] = cell_line_colours['NSC']
    row_colours_all.loc[row_colours_all.index.str.contains(r'GIBCO')] = cell_line_colours['NSC']

    row_colours_all.loc[row_colours_all.index.str.contains(r'iPS_')] = cell_line_colours['iPSC']
    row_colours_all.loc[row_colours_all.index.str.contains('HPSI')] = cell_line_colours['iPSC']
    row_colours_all.loc[row_colours_all.index.str.contains('HEL1')] = cell_line_colours['iPSC']

    row_colours_all.loc[row_colours_all.index.str.contains('neuron')] = '#ccebc5'
    row_colours_all.loc[row_colours_all.index.str.contains(r'[Aa]strocyte')] = '#e78ac3'  # pink

    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_NSC')] = cell_line_colours['iNSC (this study)']  # green
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_IPSC')] = cell_line_colours['iPSC (this study)']  # orange

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

    lkg = plt_dict[clust_n_ftr]['linkage']
    leg_dict = {
        'Cell type': collections.OrderedDict([
            (k, cell_line_colours[k]) for k in sorted(cell_line_colours)
        ])
    }
    leg_dict['Cell type'].pop('GBM')
    leg_dict['Cell type'].pop('GBM (this study)')

    cm = clustering.plot_clustermap(
        this_dat,
        cmap='RdYlBu_r',
        metric='correlation',
        col_colors=row_colours_all,
        col_linkage=lkg,
        vmin=-10,
        vmax=10
    )
    clustering.add_legend(leg_dict, cm.ax_heatmap, loc='right')
    cm.gs.update(bottom=0.25, right=0.85, left=0.03)

    cm.savefig(os.path.join(outdir, "clustermap_ipsc_esc_nsc_fb.png"), dpi=200)
    cm.savefig(os.path.join(outdir, "clustermap_ipsc_esc_nsc_fb.tiff"), dpi=200)
