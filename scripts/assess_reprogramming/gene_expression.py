from rnaseq import loader
from plotting import clustering, common
from stats import transformations
import pandas as pd
import numpy as np
import copy
import os
from utils import output
import references
from scipy.stats import zscore


class SetMe(object):
    """
    This is useful to determine which fields need to be set
    """


def load_refs(ref_dict, **load_kwds):
    ref_objs_arr = []

    for k, v in ref_dict.items():
        the_kwds = copy.copy(load_kwds)
        for k1, v1 in the_kwds.items():
            if v1 is SetMe:
                the_kwds[k1] = v.get(k1)
        the_obj = loader.load_references(k, **the_kwds)
        the_obj.batch_id = v['batch']
        ref_objs_arr.append(the_obj)

    return loader.loader.MultipleBatchLoader(ref_objs_arr)


def combine_filter(dat_arr, meta_arr, min_val=1, n_above_min=3):
    dat = pd.concat(dat_arr, axis=1).dropna(axis=0)
    meta = pd.concat(meta_arr, axis=0).loc[dat.columns]

    # combine `cell type` and `type` fields
    meta['type'].loc[meta['type'].isnull()] = meta['cell type'].loc[meta['type'].isnull()]

    # merge the batch for our data
    meta.loc[meta['batch'].str.contains('wtchg'), 'batch'] = 'WTCHG'

    # filter and QN
    dat = dat.loc[(dat > min_val).sum(axis=1) > n_above_min]

    return dat, meta


def construct_colour_array_legend_studies(meta):
    studies = {}
    cc = pd.DataFrame('gray', index=meta.index, columns=['Cell type', 'Study'])

    cc.loc[meta['type'] == 'FB', 'Cell type'] = '#fff89e'
    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'].str.contains('WTCHG')), 'Cell type'] = 'blue'
    cc.loc[(meta['type'] == 'iPSC') & (~meta['batch'].str.contains('WTCHG')), 'Cell type'] = '#96daff'
    cc.loc[meta['type'] == 'ESC', 'Cell type'] = 'green'
    cc.loc[meta['type'] == 'EPS', 'Cell type'] = '#7fc97f'
    cc.loc[(meta['type'] == 'iNSC') & (meta['batch'].str.contains('WTCHG')), 'Cell type'] = '#9e3900'  # chestnut
    cc.loc[meta['type'] == 'iNSC', 'Cell type'] = '#db7b00'  # orange
    cc.loc[meta['type'] == 'NPC', 'Cell type'] = '#db7b00'  # orange
    cc.loc[meta['type'] == 'NSC', 'Cell type'] = '#db7b00'  # orange
    cc.loc[(meta['type'] == 'NSC') & (meta.index.str.contains('fetal')), 'Cell type'] = '#ffaf47'  # light orange

    batches = meta.batch.unique()
    n_study = len(batches)
    study_colours = common.COLOUR_BREWERS[n_study]
    for i, c in enumerate(study_colours):
        cc.loc[meta['batch'] == batches[i], 'Study'] = c
        studies[batches[i]] = c

    return cc, studies


if __name__ == '__main__':
    pids = ['017', '018', '019', '030', '031', '049', '050', '054', '061', '026', '052']
    min_val = 1
    n_above_min = 3
    source = 'salmon'
    n_hipsci = 30  # only plot a limited number to avoid flooding the plots

    ref_dict = {
        'GSE80732': {'batch': 'Yang et al.', 'strandedness': 'u', 'alignment_subdir': 'human/trimgalore'},
        'GSE63577': {'batch': 'Marthandan et al.', 'strandedness': 'TODO: run STAR'},
        'GSE107965': {'batch': 'Yilmaz et al.', 'strandedness': 'TODO: run STAR'},
        'encode_roadmap/ENCSR000EYP': {'batch': 'ENCODE Wold', 'strandedness': 'u'},
        'encode_roadmap/ENCSR000COU': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        'encode_roadmap/ENCSR490SQH': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR000CRJ': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},  # low qual; outlier
        'encode_roadmap/ENCSR670WQY': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        'encode_roadmap/ENCSR043RSE': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
    }

    nsc_ref_dict = {
        'E-MTAB-3867': {'batch': 'Caren et al.', 'strandedness': '?'},
        'GSE61794': {'batch': 'Duan et al.', 'strandedness': '?'},
        'GSE64882': {'batch': 'Shahbazi et al.', 'strandedness': '?'},
        'encode_roadmap/ENCSR244ISQ': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        'encode_roadmap/ENCSR291IZK': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        'encode_roadmap/ENCSR572EET': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        'encode_roadmap/ENCSR977XUX': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
    }

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

    dat_hip, meta_hip = loader.hipsci_ipsc(aggregate_to_gene=True)
    meta_hip.insert(3, 'batch', 'HipSci')

    # reduce the number in a (repeatably) random fashion
    rs = np.random.RandomState(42)  # set the seed so we always get the same samples
    idx = np.arange(dat_hip.shape[1])
    rs.shuffle(idx)
    dat_hip_reduced = dat_hip.iloc[:, idx[:n_hipsci]]
    meta_hip_reduced = meta_hip.loc[dat_hip_reduced.columns]

    ref_obj = load_refs(ref_dict, **load_kwds)
    dat_ref = ref_obj.data.copy()
    meta_ref = ref_obj.meta.copy()

    # our data - more selective
    dat_home = obj.data.loc[:, obj.meta.type.isin(['iPSC', 'FB'])]

    # discard irrelevant samples

    to_discard = [
        '_rapamycin_',
        'HFF_PD16', 'HFF_PD46', 'HFF_PD64', 'HFF_PD74',
        'BJ_OLD',
        'IMR90_O',
        'MRC_5_PD52', 'MRC_5_PD62', 'MRC_5_PD72',
        'WI_38_O',
    ]
    for td in to_discard:
        the_idx = ~dat_ref.columns.str.contains(td)
        dat_ref = dat_ref.loc[:, the_idx]
        meta_ref = meta_ref.loc[the_idx]

    dat, meta = combine_filter(
        (dat_home, dat_ref),
        (obj.meta, ref_obj.meta),
        min_val=min_val,
        n_above_min=n_above_min,
    )
    dat_qn = transformations.quantile_normalisation(np.log2(dat + 1))

    cc, st = construct_colour_array_legend_studies(meta)

    leg_dict = {
        'Cell type': {
            'FB': '#fff89e',
            'iPSC (this study)': 'blue',
            'iPSC': '#96daff',
            'ESC': 'green',
            'Enhanced PSC': '#7fc97f',
        },
        'Study': st,
    }

    dend = clustering.dendrogram_with_colours(dat_qn, cc, vertical=True, legend_labels=leg_dict, fig_kws={'figsize': [14, 6]})
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb.png"), dpi=200)

    # Only iPSC and ESC

    dat = dat.loc[:, meta.type != 'FB']
    meta = meta.loc[dat.columns]
    dat_qn = transformations.quantile_normalisation(np.log2(dat + 1))

    cc, st = construct_colour_array_legend_studies(meta)

    leg_dict = {
        'Cell type': {
            'iPSC (this study)': 'blue',
            'iPSC': '#96daff',
            'ESC': 'green',
            'Enhanced PSC': '#7fc97f',
        },
        'Study': st,
    }

    dend = clustering.dendrogram_with_colours(dat_qn, cc, vertical=True, legend_labels=leg_dict, fig_kws={'figsize': [14, 6]})
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc.png"), dpi=200)

    # Ruiz signature (only)

    dat_r_z = np.log2(dat + 1).loc[gene_sign_ens.values].dropna().apply(zscore, axis=1)
    dat_r_z.index = gene_sign_ens.index[gene_sign_ens.isin(dat_r_z.index)]

    cg = clustering.plot_clustermap(dat_r_z, show_gene_labels=True, cmap='RdBu_r')
    cg.gs.update(bottom=0.2)
    cg.savefig(os.path.join(outdir, "clustermap_ruiz_ipsc_esc_ztrans.png"), dpi=200)

    # now add some HiPSCi and repeat

    dat, meta = combine_filter(
        (dat_home, dat_ref, dat_hip_reduced),
        (obj.meta, ref_obj.meta, meta_hip_reduced),
        min_val=min_val,
        n_above_min=n_above_min,
    )
    dat_qn = transformations.quantile_normalisation(np.log2(dat + 1))

    cc, st = construct_colour_array_legend_studies(meta)
    leg_dict = {
        'Cell type': {
            'FB': '#fff89e',
            'iPSC (this study)': 'blue',
            'iPSC': '#96daff',
            'ESC': 'green',
            'Enhanced PSC': '#7fc97f',
        },
        'Study': st,
    }

    dend = clustering.dendrogram_with_colours(
        dat_qn,
        cc,
        vertical=True,
        legend_labels=leg_dict,
        fig_kws={'figsize': [14, 6]},
        # show_labels=False
    )
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_with_hipsci%d.png" % n_hipsci), dpi=200)

    # remove FB and repeat
    dat = dat.loc[:, meta.type != 'FB']
    meta = meta.loc[dat.columns]
    dat_qn = transformations.quantile_normalisation(np.log2(dat + 1))

    cc, st = construct_colour_array_legend_studies(meta)

    leg_dict = {
        'Cell type': {
            'iPSC (this study)': 'blue',
            'iPSC': '#96daff',
            'ESC': 'green',
            'Enhanced PSC': '#7fc97f',
        },
        'Study': st,
    }

    dend = clustering.dendrogram_with_colours(dat_qn, cc, vertical=True, legend_labels=leg_dict, fig_kws={'figsize': [14, 6]})
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_with_hipsci%d.png" % n_hipsci), dpi=200)

    # Ruiz signature (only)

    dat_r_z = np.log2(dat + 1).loc[gene_sign_ens.values].dropna().apply(zscore, axis=1)
    dat_r_z.index = gene_sign_ens.index[gene_sign_ens.isin(dat_r_z.index)]

    cg = clustering.plot_clustermap(dat_r_z, show_gene_labels=True, cmap='RdBu_r', )
    cg.gs.update(bottom=0.2)
    cg.savefig(os.path.join(outdir, "clustermap_ruiz_ipsc_esc_with_hipsci%d.png" % n_hipsci), dpi=200)

    # Let's bring in our NSC and ref NSC now

    dat_home = obj.data.loc[:, obj.meta.type.isin(['iPSC', 'FB', 'NSC', 'iNSC'])]

    # load additional NSC ref data, previously withheld
    nsc_ref_obj = load_refs(nsc_ref_dict, **load_kwds)
    dat_nsc_ref = nsc_ref_obj.data

    dat, meta = combine_filter(
        (dat_home, dat_ref, dat_nsc_ref),
        (obj.meta, ref_obj.meta, nsc_ref_obj.meta),
        min_val=min_val,
        n_above_min=n_above_min
    )
    dat_qn = transformations.quantile_normalisation(np.log2(dat + 1))

    cc, st = construct_colour_array_legend_studies(meta)

    leg_dict = {
        'Cell type': {
            'FB': '#fff89e',
            'iPSC (this study)': 'blue',
            'iPSC': '#96daff',
            'ESC': 'green',
            'Enhanced PSC': '#7fc97f',
            'iNSC (this study)': '#9e3900',
            'iNSC': '#db7b00',
            'Fetal NSC': '#ffaf47'
        },
        'Study': st,
    }

    dend = clustering.dendrogram_with_colours(
        dat_qn,
        cc,
        vertical=True,
        legend_labels=leg_dict,
        fig_kws={'figsize': [16, 6]},
        # show_labels=False
    )

    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_nsc.png"), dpi=200)

    # all of that AND HipSci samples

    dat, meta = combine_filter(
        (dat_home, dat_ref, dat_nsc_ref, dat_hip_reduced),
        (obj.meta, ref_obj.meta, nsc_ref_obj.meta, meta_hip_reduced),
        min_val=min_val,
        n_above_min=n_above_min
    )
    dat_qn = transformations.quantile_normalisation(np.log2(dat + 1))

    cc, st = construct_colour_array_legend_studies(meta)

    leg_dict = {
        'Cell type': {
            'FB': '#fff89e',
            'iPSC (this study)': 'blue',
            'iPSC': '#96daff',
            'ESC': 'green',
            'Enhanced PSC': '#7fc97f',
            'iNSC (this study)': '#9e3900',
            'iNSC': '#db7b00',
            'Fetal NSC': '#ffaf47'
        },
        'Study': st,
    }

    dend = clustering.dendrogram_with_colours(
        dat_qn,
        cc,
        vertical=True,
        legend_labels=leg_dict,
        fig_kws={'figsize': [16, 6]},
    )

    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_nsc_hipsci%d.png" % n_hipsci), dpi=200)
