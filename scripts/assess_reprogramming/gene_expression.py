from rnaseq import loader
from plotting import clustering, common, scatter
from stats import transformations
import pandas as pd
import numpy as np
from copy import copy
import os
from utils import output
import references
from scipy.stats import zscore
from copy import copy
from sklearn.decomposition import pca


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


def combine_filter(dat_arr, meta_arr, min_val=1, n_above_min=3):
    dat = pd.concat(dat_arr, axis=1).dropna(axis=0)
    meta = pd.concat(meta_arr, axis=0).loc[dat.columns]

    # combine `cell type` and `type` fields
    meta['type'].loc[meta['type'].isnull()] = meta['cell type'].loc[meta['type'].isnull()]

    # merge the batch for our data
    meta.loc[meta['batch'].str.contains('wtchg'), 'batch'] = 'WTCHG'

    # filter
    dat = dat.loc[(dat > min_val).sum(axis=1) > n_above_min]

    return dat, meta


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
    pids = ['017', '018', '019', '030', '031', '049', '050', '054', '061', '026', '052']
    min_val = 1
    n_above_min = 3
    n_gene_by_mad = 5000
    source = 'salmon'
    n_hipsci = 30  # only plot a limited number to avoid flooding the plots
    eps = 0.01  # offset to use when applying log transform

    # Switch QN on here ('mean' or 'median')
    quantile_norm = 'median'

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

    outdir = output.unique_output_dir("assess_reprogramming")

    load_kwds = {'source': source, 'alignment_subdir': SetMe}
    if source == 'salmon':
        units = 'tpm'
        load_kwds['units'] = 'tpm'
    if source == 'star':
        # set strandedness as a cue to import for each
        load_kwds['strandedness'] = SetMe

    # our data (everything)
    obj = loader.load_by_patient(pids, source=source)

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

    nsc_ref_obj.rename_with_attributes(existing_attr='batch')

    # fill in missing cell types
    # only necessary if some of the references have a differently named column
    if 'cell type' in nsc_ref_obj.meta.columns:
        nsc_ref_obj.meta.loc[nsc_ref_obj.meta.type.isnull(), 'type'] = nsc_ref_obj.meta.loc[nsc_ref_obj.meta.type.isnull(), 'cell type']

    # 1. Our data: iPSC and FB only
    obj1 = copy(obj)
    ix = obj1.meta.type.isin(['iPSC', 'FB'])
    obj1.filter_samples(ix)

    obj2 = ref_obj

    dend = plot_dendrogram([obj1, obj2], qn_method=quantile_norm)
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb.png"), dpi=200)

    # 2. Only iPSC and ESC
    obj1 = copy(obj)
    ix = obj1.meta.type.isin(['iPSC'])
    obj1.filter_samples(ix)

    obj2 = copy(ref_obj)
    ix = obj2.meta.type == 'ESC'
    obj2.filter_samples(ix)

    dend = plot_dendrogram([obj1, obj2], qn_method=quantile_norm)
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc.png"), dpi=200)

    # 3. iPSC, ESC, Ruiz signature (only)
    the_obj = loader.MultipleBatchLoader([obj1, obj2])
    dat_r_z = pd.DataFrame(np.log2(the_obj.data + eps))
    dat_r_z = dat_r_z.reindex(gene_sign_ens.values).dropna()
    for r in dat_r_z.index:
        dat_r_z.loc[r] = zscore(dat_r_z.loc[r])

    dat_r_z.index = gene_sign_ens.index[gene_sign_ens.isin(dat_r_z.index)]

    cg = clustering.plot_clustermap(dat_r_z, show_gene_labels=True, cmap='RdBu_r')
    cg.gs.update(bottom=0.2)
    cg.savefig(os.path.join(outdir, "clustermap_ruiz_ipsc_esc_ztrans.png"), dpi=200)

    # 4. HipSci, iPSC, ESC, FB
    obj1 = copy(obj)
    ix = obj1.meta.type.isin(['iPSC', 'FB'])
    obj1.filter_samples(ix)

    dend = plot_dendrogram([obj1, ref_obj, hip_obj], qn_method=quantile_norm)
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_with_hipsci%d.png" % n_hipsci), dpi=200)

    # 5. HipSci, iPSC, ESC
    obj1 = copy(obj)
    ix = obj1.meta.type.isin(['iPSC'])
    obj1.filter_samples(ix)

    obj2 = copy(ref_obj)
    ix = obj2.meta.type == 'ESC'
    obj2.filter_samples(ix)

    dend = plot_dendrogram([obj1, obj2, hip_obj], qn_method=quantile_norm)
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_with_hipsci%d.png" % n_hipsci), dpi=200)

    # 6. HipSci, iPSC, ESC Ruiz signature (only)
    the_obj = loader.MultipleBatchLoader([obj1, obj2, hip_obj])

    dat_r_z = pd.DataFrame(np.log2(the_obj.data + 1))
    dat_r_z = dat_r_z.reindex(gene_sign_ens.values).dropna()
    for r in dat_r_z.index:
        dat_r_z.loc[r] = zscore(dat_r_z.loc[r])

    dat_r_z.index = gene_sign_ens.index[gene_sign_ens.isin(dat_r_z.index)]

    cg = clustering.plot_clustermap(dat_r_z, show_gene_labels=True, cmap='RdBu_r', )
    cg.gs.update(bottom=0.2)
    cg.savefig(os.path.join(outdir, "clustermap_ruiz_ipsc_esc_with_hipsci%d.png" % n_hipsci), dpi=200)

    # 7. iPSC, ESC, FB, iNSC
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

    # 8. HipSci, iPSC, ESC, FB, iNSC
    dend = plot_dendrogram(
        [obj1, ref_obj, nsc_ref_obj, hip_obj],
        vertical=False,
        figsize=(7, 12),
        qn_method=quantile_norm,
        n_by_mad=n_gene_by_mad
    )
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_nsc_hipsci%d.png" % n_hipsci), dpi=200)

    # 9. PCA with all samples
    the_obj = loader.MultipleBatchLoader([obj1, ref_obj, nsc_ref_obj, hip_obj])
    the_dat = np.log2(the_obj.data + eps)
    if quantile_norm is not None:
        the_dat = transformations.quantile_normalisation(the_dat, method=quantile_norm)
    studies = the_obj.meta.batch.copy()
    # to keep individual batches, comment out this next line
    studies[studies.str.contains('wtchg')] = 'This study'

    p, ax = plot_pca(the_dat, the_obj.meta.type, marker_subgroups=studies)
    ax.figure.subplots_adjust(right=0.78)
    ax.figure.savefig(os.path.join(outdir, "pca_all_samples.png"), dpi=200)

