"""
Added on 4th Sept 2018

Modified version of scripts.assess_reprogramming.gene_expression

Load RNA-Seq gene expression data and produce PCA and hierarchical clustering representations.
"""

import collections
import os
import sys
from copy import copy, deepcopy

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy as hc
from sklearn.decomposition import pca

from plotting import clustering, common, scatter
from rnaseq import loader
from scripts.hgic_final import consts
from stats import transformations
from utils import output, reference_genomes


class SetMe(object):
    """
    This is useful to determine which fields need to be set
    """


def load_refs(ref_dict, **load_kwds):
    ref_objs_arr = []

    for k, v in ref_dict.items():
        the_kwds = copy(load_kwds)
        for k1, v1 in the_kwds.items():
            if v1 is SetMe:
                the_kwds[k1] = v.get(k1)
        the_obj = loader.load_references(k, **the_kwds)
        the_obj.batch_id = v['batch']
        ref_objs_arr.append(the_obj)

    if len(ref_objs_arr) == 1:
        return ref_objs_arr[0]
    else:
        return loader.loader.MultipleBatchLoader(ref_objs_arr)


def filter_loader(obj, min_val=1, n_above_min=3):
    dat = obj.data
    idx = (dat > min_val).sum(axis=1) > n_above_min
    dat = dat.loc[idx]
    obj.data = dat
    return obj


def construct_colour_array_legend_studies(meta):
    ct_leg = {
        'FB': '#fff89e',
        'iPSC (this study)': '#fdc086',  # orange
        'iPSC': '#990000',  # dark red
        'ESC': '#ff7777',  # light red
    }

    cc = pd.DataFrame('gray', index=meta.index, columns=['Cell type', 'Study'])

    cols_incl = ['FB', 'ESC']
    for t in cols_incl:
        cc.loc[meta['type'] == t, 'Cell type'] = ct_leg[t]

    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'].str.contains('WTCHG')), 'Cell type'] = ct_leg['iPSC (this study)']
    cc.loc[(meta['type'] == 'iPSC') & (~meta['batch'].str.contains('WTCHG')), 'Cell type'] = ct_leg['iPSC']

    all_batches = meta.batch.copy()
    # to keep individual batches, comment out this next line
    all_batches[all_batches.str.contains('wtchg')] = 'This study'
    batches = all_batches.unique()
    n_study = len(batches)
    study_colours = common.get_best_cmap(n_study)
    studies = {}
    for i, c in enumerate(study_colours):
        cc.loc[all_batches == batches[i], 'Study'] = c
        studies[batches[i]] = c

    all_colours = cc.loc[:, 'Cell type'].unique()
    for k in ct_leg.keys():
        if ct_leg[k] not in all_colours:
            ct_leg.pop(k)

    # legend dictionary
    leg_dict = {
        'Cell type': ct_leg,
        'Study': studies
    }

    return cc, studies, leg_dict


def plot_dendrogram(
        obj_arr,
        n_by_mad=None,
        qn_method=None,
        eps=0.01,
        min_val=1,
        n_above_min=3,
        vertical=False,
        figsize=(7, 8),
        **kwargs
):
    if len(obj_arr) > 1:
        the_obj = loader.MultipleBatchLoader(obj_arr)
    else:
        the_obj = obj_arr[0]

    the_obj = filter_loader(the_obj, min_val=min_val, n_above_min=n_above_min)
    dat = np.log2(the_obj.data + eps)
    if qn_method is not None:
        dat = transformations.quantile_normalisation(dat, method=qn_method)

    if n_by_mad is not None:
        mad = transformations.median_absolute_deviation(dat).sort_values(ascending=False)
        dat = dat.loc[mad.index[:n_by_mad]]

    cc, st, leg_dict = construct_colour_array_legend_studies(the_obj.meta)

    dend = clustering.dendrogram_with_colours(
        dat,
        cc,
        vertical=vertical,
        legend_labels=leg_dict,
        fig_kws={'figsize': figsize},
        **kwargs
    )

    return dend


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


def format_clustermap(gc, title=r'$\log_{10}(\mathrm{TPM})$'):
    """
    Separate formatting process so that we can call it from the methylation side too
    :param gc:
    :return:
    """
    gc.fig.set_size_inches((9.4, 8.))
    # reduce cbar ticks
    cbar = gc.ax_heatmap.collections[0].colorbar
    cbar.set_ticks([-5, 0, 5])

    plt.setp(gc.cax.yaxis.get_ticklabels(), fontsize=12)
    plt.setp(gc.ax_col_colors.yaxis.get_ticklabels(), fontsize=12)
    gc.cax.yaxis.set_ticks_position('left')
    gc.cax.set_title(title)

    gc.gs.set_height_ratios([0.16, 0.04, 0.08, 0.8])
    gc.gs.update(bottom=0.3, right=0.73, left=0.05, wspace=0.0, top=0.95)


def plot_clustermap(
    obj,
    quantile_norm,
    method='average',
    metric='correlation',
    n_gene_by_mad=5000,
    n_gene_for_heatmap=500,
    fmin=0.05,
    fmax=0.95,
    eps=0.01,
    cell_line_colours=None
):
    if cell_line_colours is None:
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

    the_dat = np.log2(obj.data + eps)

    if quantile_norm is not None:
        the_dat = transformations.quantile_normalisation(the_dat, method=quantile_norm)
    the_mad = transformations.median_absolute_deviation(the_dat).sort_values(ascending=False)
    cc, st, leg_dict = construct_colour_array_legend_studies(obj.meta)

    # linkage
    lkg = hc.linkage(
        the_dat.loc[the_mad.index[:n_gene_by_mad]].transpose(),
        method=method,
        metric=metric,
    )

    # ref line colours
    for k, v in cell_line_colours.items():
        cc.loc[obj.meta.type == k, 'Cell type'] = v
    # our line colours
    cc.loc[obj.meta.batch.str.contains('wtchg') & (obj.meta.type == 'iPSC'), 'Cell type'] = \
    cell_line_colours['iPSC (this study)']

    # get appropriate clims
    the_dat = the_dat.loc[the_mad.index[:n_for_heatmap]]
    the_dat_flat = np.sort(the_dat.values.flatten())
    vmin = the_dat_flat[int(len(the_dat_flat) * fmin)] - 0.5
    vmax = the_dat_flat[int(len(the_dat_flat) * fmax)] + 0.5

    gc = clustering.plot_clustermap(
        the_dat.loc[the_mad.index[:n_gene_for_heatmap]],
        cmap='RdBu_r',
        col_linkage=lkg,
        col_colors=cc,
        vmin=vmin,
        vmax=vmax,
    )

    leg_entry = {
        'class': 'patch',
        'edgecolor': 'k',
        'linewidth': 1.,
    }
    leg_dict2 = collections.OrderedDict()
    leg_dict2['Cell type'] = collections.OrderedDict()

    for k in sorted(cell_line_colours):
        if k.replace(' (this study)', '') in obj.meta.type.unique():
            leg_dict2['Cell type'][k] = dict(leg_entry)
            leg_dict2['Cell type'][k].update({'facecolor': cell_line_colours[k]})

    leg_dict2['Study'] = {}
    for k, v in leg_dict['Study'].items():
        leg_dict2['Study'][k] = dict(leg_entry)
        leg_dict2['Study'][k].update({'facecolor': v})

    common.add_custom_legend(gc.ax_heatmap, leg_dict2, loc_outside=True, fontsize=14)
    format_clustermap(gc)

    return gc


if __name__ == '__main__':
    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    outdir = output.unique_output_dir(script_name)

    pids = consts.PIDS
    min_val = 1
    n_above_min = 3
    n_gene_by_mad = 5000
    n_for_heatmap = 500
    source = 'salmon'
    n_hipsci = 10  # only plot a limited number to avoid flooding the plots
    eps = 0.01  # offset to use when applying log transform

    # Switch QN on here ('mean' or 'median')
    quantile_norm = 'median'

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

    ref_dict = {
        'GSE80732': {'batch': 'Yang et al.', 'strandedness': 'u', 'alignment_subdir': 'human/trimgalore'},
        # 'GSE63577': {'batch': 'Marthandan et al.', 'strandedness': 'TODO: run STAR'},
        'GSE107965': {'batch': 'Yilmaz et al.', 'strandedness': 'TODO: run STAR'},
        'encode_roadmap/ENCSR000EYP': {'batch': 'ENCODE Wold', 'strandedness': 'u'},
        'encode_roadmap/ENCSR000COU': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        'encode_roadmap/ENCSR490SQH': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR000CRJ': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},  # low qual; outlier
        'encode_roadmap/ENCSR670WQY': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        'encode_roadmap/ENCSR043RSE': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
    }
    to_aggr = [
        (r'Fibroblasts_control_rep[123]', 'FB control'),
        (r'H1-hESC rep (1_1|2_1|3|4)', 'H1 hESC'),
        (r'H1-hESC rep [12]', 'H1 hESC'),
        (r'H1-[12]', 'H1 hESC'),
        (r'H1-hESC(|_2)$', 'H1 hESC'),
        (r'H7-hESC rep [12]', 'H7 hESC'),
        (r'hESCs_control_rep[123]', 'CSES9 hESC'),
    ]

    # Ruiz 9 gene signature - should distinguish ESC and iPSC
    gene_sign = ['PTPRT', 'TMEM132C', 'TMEM132D', 'TCERG1L', 'DPP6', 'FAM19A5', 'RBFOX1', 'CSMD1', 'C22orf34']
    gene_sign_ens = reference_genomes.gene_symbol_to_ensembl(gene_sign)

    load_kwds = {'source': source, 'alignment_subdir': SetMe}
    if source == 'salmon':
        units = 'tpm'
        load_kwds['units'] = 'tpm'
    if source == 'star':
        # set strandedness as a cue to import for each
        load_kwds['strandedness'] = SetMe

    # restrict samples manually to avoid changes going forwards
    our_samples = consts.S1_RNASEQ_SAMPLES_IPSC + consts.S1_RNASEQ_SAMPLES_FB

    # our data (everything)

    obj = loader.load_by_patient(pids, source=source)
    # obj.filter_by_sample_name(our_samples)

    # HipSci data
    hip_obj = loader.hipsci_ipsc(aggregate_to_gene=True)
    hip_obj.meta.insert(3, 'batch', hip_obj.batch_id)
    # hip_obj.meta.insert(3, 'batch', 'HipSci')

    # reduce the number in a (repeatably) random fashion
    rs = np.random.RandomState(42)  # set the seed so we always get the same samples
    keep = np.zeros(hip_obj.meta.shape[0]).astype(bool)
    idx = hip_obj.meta.index.tolist()
    rs.shuffle(idx)
    idx = idx[:n_hipsci]
    hip_obj.meta = hip_obj.meta.loc[idx]
    hip_obj.data = hip_obj.data.loc[:, idx]

    # References without NSC
    ref_obj = load_refs(ref_dict, **load_kwds)

    # discard irrelevant samples
    to_discard = [
        '_rapamycin_',
        'HFF_PD16', 'HFF_PD46', 'HFF_PD64', 'HFF_PD74',
        'BJ_OLD',
        'IMR90_O',
        'MRC_5_PD52', 'MRC_5_PD62', 'MRC_5_PD72',
        'WI_38_O',
        '-EPS-',
    ]
    for td in to_discard:
        the_idx = ~ref_obj.data.columns.str.contains(td)
        ref_obj.filter_samples(the_idx)

    # fill in missing cell types
    ref_obj.meta.loc[ref_obj.meta.type.isnull(), 'type'] = ref_obj.meta.loc[ref_obj.meta.type.isnull(), 'cell type']

    for srch, repl in to_aggr:
        ref_obj.aggregate_by_pattern(srch, repl)

    ref_obj.rename_with_attributes(existing_attr='batch')

    # 1a. iPSC, ESC, FB
    obj1 = copy(obj)
    ix = obj1.meta.type.isin(['iPSC', 'FB'])
    obj1.filter_samples(ix)

    the_obj = loader.loader.MultipleBatchLoader([obj1, ref_obj])
    # manually rename our samples for nicer plotting
    the_obj.meta.index = the_obj.meta.index.str.replace(r'DURA(?P<n>[0-9]{3})_(?P<typ>[^_]*).*', r'\g<typ>\g<n>')
    the_obj.data.columns = the_obj.meta.index

    gc = plot_clustermap(
        the_obj,
        quantile_norm,
        n_gene_by_mad=n_gene_by_mad,
        n_gene_for_heatmap=n_for_heatmap,
        cell_line_colours=cell_line_colours,
        fmin=0.02,
        fmax=0.98
    )

    plt.setp(gc.ax_heatmap.get_xticklabels(), fontsize=13)
    gc.gs.update(bottom=0.33)

    gc.savefig(os.path.join(outdir, "clustermap_ipsc_esc_fb.png"), dpi=200)
    gc.savefig(os.path.join(outdir, "clustermap_ipsc_esc_fb.tiff"), dpi=200)

    # 2. HipSci, iPSC, ESC, FB, iNSC
    obj2 = deepcopy(hip_obj)
    obj2.meta.index = ["%s (HipSci)" % t for t in obj2.meta.index]
    obj2.data.columns = obj2.meta.index
    the_obj = loader.MultipleBatchLoader([obj1, ref_obj, obj2])
    # manually rename samples for nicer plotting
    the_obj.meta.index = the_obj.meta.index.str.replace(r'DURA(?P<n>[0-9]{3})_(?P<typ>[^_]*).*', r'\g<typ>\g<n>')
    the_obj.data.columns = the_obj.meta.index

    gc = plot_clustermap(
        the_obj,
        quantile_norm,
        n_gene_by_mad=n_gene_by_mad,
        n_gene_for_heatmap=n_for_heatmap,
        cell_line_colours=cell_line_colours,
        fmin=0.02,
        fmax=0.98
    )

    gc.savefig(os.path.join(outdir, "clustermap_ipsc_hipsci_esc_fb.png"), dpi=200)
    gc.savefig(os.path.join(outdir, "clustermap_ipsc_hipsci_esc_fb.tiff"), dpi=200)

    # 3. PCA with all samples
    the_obj = loader.MultipleBatchLoader([obj1, ref_obj, obj2])
    the_dat = np.log2(the_obj.data + eps)
    if quantile_norm is not None:
        the_dat = transformations.quantile_normalisation(the_dat, method=quantile_norm)
    studies = the_obj.meta.batch.copy()
    # to keep individual batches, comment out this next line
    studies[studies.str.contains('wtchg')] = 'This study'

    p, ax = plot_pca(the_dat, the_obj.meta.type, marker_subgroups=studies)
    ax.figure.subplots_adjust(right=0.78)
    ax.figure.savefig(os.path.join(outdir, "pca_ipsc_esc_fb_hipsci%d.png" % n_hipsci), dpi=200)

