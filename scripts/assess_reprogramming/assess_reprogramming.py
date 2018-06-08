from rnaseq import loader, general, filter
from plotting import clustering, common
from stats import transformations
import pandas as pd
import numpy as np
import copy
import os
from utils import output


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


if __name__ == '__main__':
    pids = ['017', '018', '019', '030', '031', '049', '050', '054', '061', '026', '052']
    min_val = 1
    source = 'salmon'

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

    outdir = output.unique_output_dir("assess_reprogramming")

    load_kwds = {'source': source, 'alignment_subdir': SetMe}
    if source == 'salmon':
        units = 'tpm'
        load_kwds['units'] = 'tpm'
    if source == 'star':
        # set strandedness as a cue to import for each
        load_kwds['strandedness'] = SetMe

    obj = loader.load_by_patient(pids, source=source)
    dat_home = obj.data.loc[:, obj.meta.type.isin(['iPSC', 'FB'])]

    dat_hip, meta_hip = loader.hipsci_ipsc(aggregate_to_gene=True)
    meta_hip.insert(3, 'batch', 'HipSci')

    ref_obj = load_refs(ref_dict, **load_kwds)
    dat_ref = ref_obj.data

    # discard irrelevant samples

    to_discard = [
        '_rapamycin_',
    ]
    for td in to_discard:
        the_idx = ~dat_ref.columns.str.contains(td)
        dat_ref = dat_ref.loc[:, the_idx]

    dat = pd.concat((dat_home, dat_ref), axis=1).dropna(axis=0)
    meta = pd.concat((obj.meta, ref_obj.meta), axis=0).loc[dat.columns]
    meta['type'].loc[meta['type'].isnull()] = meta['cell type'].loc[meta['type'].isnull()]
    # merge the batch for our data
    meta.loc[meta['batch'].str.contains('wtchg'), 'batch'] = 'WTCHG'

    dat = dat.loc[(dat > min_val).sum(axis=1) > 6]
    # dat_qn = transformations.quantile_normalisation(dat)
    dat_qn_log = transformations.quantile_normalisation(np.log2(dat + 1))

    cc = pd.DataFrame('gray', index=dat.columns, columns=['Cell type', 'Study'])

    cc.loc[meta['type'] == 'FB', 'Cell type'] = '#fff89e'
    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'].str.contains('WTCHG')), 'Cell type'] = 'blue'
    cc.loc[(meta['type'] == 'iPSC') & (~meta['batch'].str.contains('WTCHG')), 'Cell type'] = '#96daff'
    cc.loc[meta['type'] == 'ESC', 'Cell type'] = 'green'
    cc.loc[meta['type'] == 'EPS', 'Cell type'] = '#7fc97f'

    leg_dict = {
        'Cell type': {
            'FB': '#fff89e',
            'iPSC (this study)': 'blue',
            'iPSC': '#96daff',
            'ESC': 'green',
            'Enhanced PSC': '#7fc97f',
        },
        'Study': {},
    }

    batches = meta.batch.unique()
    n_study = len(batches)
    study_colours = common.COLOUR_BREWERS[n_study]
    for i, c in enumerate(study_colours):
        cc.loc[meta['batch'] == batches[i], 'Study'] = c
        leg_dict['Study'][batches[i]] = c

    dend = clustering.dendrogram_with_colours(dat_qn_log, cc, vertical=True, legend_labels=leg_dict, fig_kws={'figsize': [14, 6]})
    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb.png"), dpi=200)

    # now add some HiPSCi and repeat
    n_hipsci = 50
    rs = np.random.RandomState(42)
    idx = np.arange(dat_hip.shape[1])
    rs.shuffle(idx)
    this_dat_hip = dat_hip.iloc[:, idx[:n_hipsci]]
    this_meta_hip = meta_hip.loc[this_dat_hip.columns]

    dat = pd.concat((dat_home, dat_ref, this_dat_hip), axis=1).dropna(axis=0)
    meta = pd.concat((obj.meta, ref_obj.meta, this_meta_hip), axis=0).loc[dat.columns]
    meta['type'].loc[meta['type'].isnull()] = meta['cell type'].loc[meta['type'].isnull()]
    # merge the batch for our data
    meta.loc[meta['batch'].str.contains('wtchg'), 'batch'] = 'WTCHG'

    dat = dat.loc[(dat > min_val).sum(axis=1) > 6]
    dat_qn_log = transformations.quantile_normalisation(np.log2(dat + 1))

    cc = pd.DataFrame('gray', index=dat.columns, columns=['Cell type', 'Study'])

    cc.loc[meta['type'] == 'FB', 'Cell type'] = '#fff89e'
    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'].str.contains('WTCHG')), 'Cell type'] = 'blue'
    cc.loc[(meta['type'] == 'iPSC') & (~meta['batch'].str.contains('WTCHG')), 'Cell type'] = '#96daff'
    cc.loc[meta['type'] == 'ESC', 'Cell type'] = 'green'
    cc.loc[meta['type'] == 'EPS', 'Cell type'] = '#7fc97f'

    batches = meta.batch.unique()
    n_study = len(batches)
    study_colours = common.COLOUR_BREWERS[n_study]
    for i, c in enumerate(study_colours):
        cc.loc[meta['batch'] == batches[i], 'Study'] = c
        leg_dict['Study'][batches[i]] = c

    dend = clustering.dendrogram_with_colours(
        dat_qn_log,
        cc,
        vertical=True,
        legend_labels=leg_dict,
        fig_kws={'figsize': [14, 6]},
        show_labels=False
    )

    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_with_hipsci%d.png" % n_hipsci), dpi=200)

    # Let's bring in our NSC and ref NSC now?

    # reload our data, because we want the iNSC and Gibco included
    obj = loader.load_by_patient(pids, source=source)
    dat_home = obj.data.loc[:, obj.meta.type.isin(['iPSC', 'FB', 'NSC', 'iNSC'])]

    # load additional NSC ref data, previously withheld
    nsc_ref_obj = load_refs(nsc_ref_dict, **load_kwds)
    dat_nsc_ref = nsc_ref_obj.data

    dat = pd.concat((dat_home, dat_ref, dat_nsc_ref), axis=1).dropna(axis=0)
    meta = pd.concat((obj.meta, ref_obj.meta, nsc_ref_obj.meta), axis=0).loc[dat.columns]
    meta['type'].loc[meta['type'].isnull()] = meta['cell type'].loc[meta['type'].isnull()]
    # merge the batch for our data
    meta.loc[meta['batch'].str.contains('wtchg'), 'batch'] = 'WTCHG'

    dat = dat.loc[(dat > min_val).sum(axis=1) > 6]
    dat_qn_log = transformations.quantile_normalisation(np.log2(dat + 1))

    leg_dict['Cell type'].update({
        'iNSC (this study)': '#9e3900',
        'iNSC': '#db7b00',
        'Fetal NSC': '#ffaf47'
    })

    cc = pd.DataFrame('gray', index=dat.columns, columns=['Cell type', 'Study'])

    cc.loc[meta['type'] == 'FB', 'Cell type'] = '#fff89e'
    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'].str.contains('WTCHG')), 'Cell type'] = 'blue'
    cc.loc[(meta['type'] == 'iPSC') & (~meta['batch'].str.contains('WTCHG')), 'Cell type'] = '#96daff'
    cc.loc[meta['type'] == 'ESC', 'Cell type'] = 'green'
    cc.loc[meta['type'] == 'EPS', 'Cell type'] = '#7fc97f'
    cc.loc[meta['type'] == 'fetal NSC', 'Cell type'] = '#ffaf47'  # light orange
    cc.loc[(meta['type'] == 'iNSC') & (meta['batch'].str.contains('WTCHG')), 'Cell type'] = '#9e3900'  # chestnut
    cc.loc[meta['type'] == 'NPC', 'Cell type'] = '#db7b00'  # orange
    cc.loc[meta['type'] == 'NSC', 'Cell type'] = '#db7b00'  # orange

    batches = meta.batch.unique()
    n_study = len(batches)
    study_colours = common.COLOUR_BREWERS[n_study]
    for i, c in enumerate(study_colours):
        cc.loc[meta['batch'] == batches[i], 'Study'] = c
        leg_dict['Study'][batches[i]] = c

    leg_dict = {
        'Cell type': {
            'FB': '#fff89e',
            'iPSC (this study)': 'blue',
            'iPSC': '#96daff',
            'ESC': 'green',
            'Enhanced PSC': '#7fc97f',
        },
        'Study': {},
    }

    batches = meta.batch.unique()
    n_study = len(batches)
    study_colours = common.COLOUR_BREWERS[n_study]
    for i, c in enumerate(study_colours):
        cc.loc[meta['batch'] == batches[i], 'Study'] = c
        leg_dict['Study'][batches[i]] = c

    dend = clustering.dendrogram_with_colours(
        dat_qn_log,
        cc,
        vertical=True,
        legend_labels=leg_dict,
        fig_kws={'figsize': [16, 6]},
        # show_labels=False
    )

    dend['fig'].savefig(os.path.join(outdir, "cluster_ipsc_esc_fb_nsc.png"), dpi=200)