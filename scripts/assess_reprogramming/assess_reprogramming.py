from rnaseq import loader, general, filter
from plotting import clustering, common
from stats import transformations
import pandas as pd
import numpy as np
import copy


class SetMe(object):
    """
    This is useful to determine which fields need to be set
    """


if __name__ == '__main__':
    pids = ['017', '018', '019', '030', '031', '049', '050', '054', '061', '026', '052']
    min_val = 1
    source = 'salmon'
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

    ref_dict = {
        'GSE80732': {'batch': 'Yang et al.', 'strandedness': 'u', 'alignment_subdir': 'human/trimgalore'},
        # 'GSE38993': {'batch': 'Kelley and Rinn', 'strandedness': 'u'},
        'encode_roadmap/ENCSR000EYP': {'batch': 'ENCODE Wold', 'strandedness': 'u'},
        'encode_roadmap/ENCSR000COU': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        'encode_roadmap/ENCSR490SQH': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        'encode_roadmap/ENCSR670WQY': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        'encode_roadmap/ENCSR043RSE': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
    }

    ref_objs_arr = []

    for k, v in ref_dict.items():
        the_kwds = copy.copy(load_kwds)
        for k1, v1 in the_kwds.items():
            if v1 is SetMe:
                the_kwds[k1] = v.get(k1)
        the_obj = loader.load_references(k, **the_kwds)
        the_obj.batch_id = v['batch']
        ref_objs_arr.append(the_obj)

    ref_obj = loader.loader.MultipleBatchLoader(ref_objs_arr)
    dat_ref = ref_obj.data

    # discard irrelevant samples

    to_discard = [
        # 'INSC fibroblast',
        # 'fetal NSC',
        'H1 NSC',
    ]
    for td in to_discard:
        the_idx = ~dat_ref.columns.str.contains(td)
        dat_ref = dat_ref.loc[:, the_idx]

    dat = pd.concat((dat_home, dat_ref), axis=1).dropna(axis=0)
    meta = pd.concat((obj.meta, ref_obj.meta), axis=0).loc[dat.columns]
    meta['type'].loc[meta['type'].isnull()] = meta['cell type'].loc[meta['type'].isnull()]

    dat = dat.loc[(dat > min_val).sum(axis=1) > 6]
    dat_qn = transformations.quantile_normalisation(dat)
    dat_qn_log = transformations.quantile_normalisation(np.log2(dat + 1))

    cc = pd.DataFrame('gray', index=dat.columns, columns=['Cell type'])
    cc.loc[meta['type'] == 'FB'] = '#fff89e'
    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'].str.contains('wtchg'))] = 'blue'
    cc.loc[(meta['type'] == 'iPSC') & (~meta['batch'].str.contains('wtchg'))] = '#96daff'
    cc.loc[meta['type'] == 'ESC'] = 'green'
    cc.loc[meta['type'] == 'EPS'] = '#7fc97f'

    dend = clustering.dendrogram_with_colours(dat_qn_log, cc, vertical=False)

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

    dat = dat.loc[(dat > min_val).sum(axis=1) > 6]
    dat_qn_log = transformations.quantile_normalisation(np.log2(dat + 1))

    cc = pd.DataFrame('gray', index=dat.columns, columns=['Cell type', 'Study'])

    cc.loc[meta['type'] == 'FB', 'Cell type'] = '#fff89e'
    cc.loc[(meta['type'] == 'iPSC') & (meta['batch'].str.contains('wtchg')), 'Cell type'] = 'blue'
    cc.loc[(meta['type'] == 'iPSC') & (~meta['batch'].str.contains('wtchg')), 'Cell type'] = '#96daff'
    cc.loc[meta['type'] == 'ESC', 'Cell type'] = 'green'
    cc.loc[meta['type'] == 'EPS', 'Cell type'] = '#7fc97f'

    batches = meta.batch.unique()
    n_study = len(batches)
    study_colours = common.COLOUR_BREWERS[n_study]
    for i, c in enumerate(study_colours):
        cc.loc[meta['batch'] == batches[i], 'Study'] = c

    dend = clustering.dendrogram_with_colours(dat_qn_log, cc, vertical=False)