"""
Added on 4th Sept 2018

Modified version of scripts.assess_reprogramming.gene_expression

Load RNA-Seq gene expression data and produce PCA and hierarchical clustering representations.
"""


from rnaseq import loader
from plotting import clustering, common, scatter
from stats import transformations
import pandas as pd
import collections
import numpy as np
import os, sys
from utils import output
import references
from copy import copy
from sklearn.decomposition import pca
import consts


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
        'iPSC (this study)': 'blue',
        'iPSC': '#96daff',
        'ESC': 'green',
        'iNSC (this study)': '#9e3900',
        'iNSC': '#db7b00',
        'NSC': '#f4b342',
        'Fetal NSC': '#ffaf47',
    }

    cc = pd.DataFrame('gray', index=meta.index, columns=['Cell type', 'Study'])

    cols_incl = ['FB', 'ESC', 'NSC']
    for t in cols_incl:
        cc.loc[meta['type'] == t, 'Cell type'] = ct_leg[t]

    cc.loc[meta['type'] == 'NPC', 'Cell type'] = ct_leg['NSC']
    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'].str.contains('WTCHG')), 'Cell type'] = ct_leg['iPSC (this study)']
    cc.loc[(meta['type'] == 'iPSC') & (~meta['batch'].str.contains('WTCHG')), 'Cell type'] = ct_leg['iPSC']
    cc.loc[(meta['type'] == 'iNSC') & (meta['batch'].str.contains('WTCHG')), 'Cell type'] = ct_leg['iNSC (this study)']  # chestnut
    cc.loc[(meta['type'] == 'iNSC') & (~meta['batch'].str.contains('WTCHG')), 'Cell type'] = ct_leg['iNSC']  # chestnut
    cc.loc[(meta['type'] == 'NSC') & (meta.index.str.contains('fetal')), 'Cell type'] = ct_leg['Fetal NSC']  # light orange

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


if __name__ == '__main__':
    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    outdir = output.unique_output_dir(script_name)

    pids = consts.PIDS
    min_val = 1
    n_above_min = 3
    n_gene_by_mad = 5000
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

    nsc_ref_dict = {
        # 'E-MTAB-3867': {'batch': 'Caren et al.', 'strandedness': '?'},
        'GSE61794': {'batch': 'Duan et al.', 'strandedness': '?'},
        # 'GSE64882': {'batch': 'Shahbazi et al.', 'strandedness': '?'},
        # 'encode_roadmap/ENCSR244ISQ': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR291IZK': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR572EET': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR977XUX': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
    }

    to_aggr_nsc = [
        (r'H9_NSC_[12]', 'H9 NSC'),
        # (r'Pollard NSC [12]', 'Fetal NSC'),
    ]

    # Ruiz 9 gene signature - should distinguish ESC and iPSC
    gene_sign = ['PTPRT', 'TMEM132C', 'TMEM132D', 'TCERG1L', 'DPP6', 'FAM19A5', 'RBFOX1', 'CSMD1', 'C22orf34']
    gene_sign_ens = references.gene_symbol_to_ensembl(gene_sign)

    load_kwds = {'source': source, 'alignment_subdir': SetMe}
    if source == 'salmon':
        units = 'tpm'
        load_kwds['units'] = 'tpm'
    if source == 'star':
        # set strandedness as a cue to import for each
        load_kwds['strandedness'] = SetMe

    # restrict samples manually to avoid changes going forwards
    our_samples = consts.S1_RNASEQ_SAMPLES_INSC + consts.S1_RNASEQ_SAMPLES_IPSC + consts.S1_RNASEQ_SAMPLES_FB + ['GIBCO_NSC_P4']

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

    # References with NSC
    nsc_ref_obj = load_refs(nsc_ref_dict, **load_kwds)
    # discard
    to_discard = [
        'INSC fibroblast',
    ]
    for td in to_discard:
        the_idx = ~nsc_ref_obj.data.columns.str.contains(td)
        nsc_ref_obj.filter_samples(the_idx)

    for srch, repl in to_aggr_nsc:
        nsc_ref_obj.aggregate_by_pattern(srch, repl)

    if 'batch' not in nsc_ref_obj.meta.columns:
        batch = pd.Series(nsc_ref_obj.batch_id, index=nsc_ref_obj.meta.index)
        nsc_ref_obj.meta.insert(nsc_ref_obj.meta.shape[1], 'batch', batch)
    nsc_ref_obj.rename_with_attributes(existing_attr='batch')

    # fill in missing cell types
    # only necessary if some of the references have a differently named column
    if 'cell type' in nsc_ref_obj.meta.columns:
        nsc_ref_obj.meta.loc[nsc_ref_obj.meta.type.isnull(), 'type'] = nsc_ref_obj.meta.loc[nsc_ref_obj.meta.type.isnull(), 'cell type']

    # 1a. iPSC, ESC, FB, iNSC
    obj1 = copy(obj)
    ix = obj1.meta.type.isin(['iPSC', 'FB', 'iNSC', 'NSC'])
    obj1.filter_samples(ix)

    dend = plot_dendrogram(
        [obj1, ref_obj, nsc_ref_obj],
        vertical=False,
        figsize=(7, 14),
        qn_method=quantile_norm,
        n_by_mad=n_gene_by_mad
    )
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_nsc.png"), dpi=200)

    # 1b. Heatmap from clustering result of (1a)
    n_for_heatmap = 500
    the_obj = loader.MultipleBatchLoader([obj1, ref_obj, nsc_ref_obj])
    the_dat = np.log2(the_obj.data + eps)
    if quantile_norm is not None:
        the_dat = transformations.quantile_normalisation(the_dat, method=quantile_norm)
    the_mad = transformations.median_absolute_deviation(the_dat).sort_values(ascending=False)
    cc, st, leg_dict = construct_colour_array_legend_studies(the_obj.meta)

    # ref line colours
    for k, v in cell_line_colours.items():
        cc.loc[the_obj.meta.type == k, 'Cell type'] = v
    # our line colours
    cc.loc[the_obj.meta.batch.str.contains('wtchg') & (the_obj.meta.type == 'iNSC'), 'Cell type'] = cell_line_colours['iNSC (this study)']
    cc.loc[the_obj.meta.batch.str.contains('wtchg') & (the_obj.meta.type == 'iPSC'), 'Cell type'] = cell_line_colours['iPSC (this study)']

    # get appropriate clims
    the_dat = the_dat.loc[the_mad.index[:n_for_heatmap]]
    the_dat_flat = np.sort(the_dat.values.flatten())
    fmin = 0.05
    fmax = 0.95
    vmin = the_dat_flat[int(len(the_dat_flat) * fmin)] - 0.5
    vmax = the_dat_flat[int(len(the_dat_flat) * fmax)] + 0.5

    gc = clustering.plot_clustermap(
        the_dat.loc[the_mad.index[:n_for_heatmap]],
        cmap='RdBu_r',
        col_linkage=dend['linkage'],
        col_colors=cc,
        vmin=vmin,
        vmax=vmax
    )

    leg_entry = {
        'class': 'patch',
        'edgecolor': 'k',
        'linewidth': 1.,
    }
    leg_dict2 = collections.OrderedDict()
    leg_dict2['Cell type'] = collections.OrderedDict()

    for k in sorted(cell_line_colours):
        if k.replace(' (this study)', '') in the_obj.meta.type.unique():
            leg_dict2['Cell type'][k] = dict(leg_entry)
            leg_dict2['Cell type'][k].update({'facecolor': cell_line_colours[k]})

    leg_dict2['Study'] = {}
    for k, v in leg_dict['Study'].items():
        leg_dict2['Study'][k] = dict(leg_entry)
        leg_dict2['Study'][k].update({'facecolor': v})

    common.add_custom_legend(gc.ax_heatmap, leg_dict2, loc_outside=True, fontsize=14)
    gc.fig.set_size_inches((10.9, 8.))
    gc.gs.update(bottom=0.3, right=0.77, left=0.01, wspace=0.03)

    # clustering.add_legend(leg_dict, gc.ax_heatmap, loc='right')
    gc.savefig(os.path.join(outdir, "clustermap_ipsc_esc_fb_nsc.png"), dpi=200)
    gc.savefig(os.path.join(outdir, "clustermap_ipsc_esc_fb_nsc.tiff"), dpi=200)

    # 2. HipSci, iPSC, ESC, FB, iNSC
    dend = plot_dendrogram(
        [obj1, ref_obj, nsc_ref_obj, hip_obj],
        vertical=False,
        figsize=(7, 12),
        qn_method=quantile_norm,
        n_by_mad=n_gene_by_mad
    )
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_nsc_hipsci%d.png" % n_hipsci), dpi=200)

    # 3. PCA with all samples
    the_obj = loader.MultipleBatchLoader([obj1, ref_obj, nsc_ref_obj, hip_obj])
    the_dat = np.log2(the_obj.data + eps)
    if quantile_norm is not None:
        the_dat = transformations.quantile_normalisation(the_dat, method=quantile_norm)
    studies = the_obj.meta.batch.copy()
    # to keep individual batches, comment out this next line
    studies[studies.str.contains('wtchg')] = 'This study'

    p, ax = plot_pca(the_dat, the_obj.meta.type, marker_subgroups=studies)
    ax.figure.subplots_adjust(right=0.78)
    ax.figure.savefig(os.path.join(outdir, "pca_ipsc_esc_fb_nsc_hipsci%d.png" % n_hipsci), dpi=200)

