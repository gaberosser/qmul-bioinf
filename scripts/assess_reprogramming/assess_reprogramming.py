from rnaseq import loader, general, filter
from plotting import clustering
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

    ref_dict = {
        'GSE80732': {'batch': 'Yang et al.', 'strandedness': 'u', 'alignment_subdir': 'human/trimgalore'},
        'GSE38993': {'batch': 'Kelley and Rinn', 'strandedness': 'u'},
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
        # 'H1 NSC',
    ]
    for td in to_discard:
        the_idx = ~ref.columns.str.contains(td)
        ref = ref.loc[:, the_idx]
        batches = batches.loc[the_idx]
        labels = labels.loc[the_idx]

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

    # now add all HiPSCi and repeat - don't plot sample names for clarity
    dat = pd.concat((dat_home, dat_ref, dat_hip), axis=1).dropna(axis=0)
    dat = dat.loc[(dat > min_val).sum(axis=1) > 6]
    dat_qn = transformations.quantile_normalisation(dat)