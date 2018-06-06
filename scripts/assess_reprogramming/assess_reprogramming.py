from rnaseq import loader, general, filter
from plotting import clustering
from stats import transformations
import pandas as pd
import numpy as np


if __name__ == '__main__':
    pids = ['017', '018', '019', '030', '031', '049', '050', '054', '061', '026', '052']
    min_val = 1
    source = 'star'
    if source == 'salmon':
        units = 'tpm'
    else:
        units = None

    obj = loader.load_by_patient(pids, source=source)
    dat_home = obj.data.loc[:, obj.meta.type.isin(['iPSC', 'FB'])]

    ref_dats = [
        (None, 'Yang et al.', loader.load_references('GSE80732', source=source, units=units),),
        (None, 'Shahbazi et al.', loader.load_references('GSE64882', source=source, units=units),),
        (None, 'Kelley and Rinn', loader.load_references('GSE38993', source=source, units=units),),
        (None, 'ENCODE Wold', loader.load_references('encode_roadmap/ENCSR000EYP', source=source, units=units),),
        (['H1 PSC Costello_PSB'], 'ENCODE Costello', loader.load_references('encode_roadmap/ENCSR950PSB', source=source, units=units),),
        (['ENCSR000COU rep 1', 'ENCSR000COU rep 2'], 'ENCODE Gingeras', loader.load_references('encode_roadmap/ENCSR000COU', source=source, units=units),),
        (None, 'ENCODE Gingeras', loader.load_references('encode_roadmap/ENCSR490SQH', source=source, units=units),),
        (['H1 PSC Ecker_WQY'], 'ENCODE Ecker', loader.load_references('encode_roadmap/ENCSR670WQY', source=source, units=units),),
        (['H1 PSC Ecker_RSE'], 'ENCODE Ecker', loader.load_references('encode_roadmap/ENCSR043RSE', source=source, units=units),),
    ]

    ref_arr = []
    ref_labels = []
    ref_cols = []
    batches = []

    for i, r in enumerate(ref_dats):
        the_dat = r[2].data
        if r[0] is not None:
            the_cols = r[0]
        else:
            the_cols = the_dat.columns
        ref_cols.extend(the_cols)
        ref_labels.extend(the_dat.columns)
        ref_arr.append(the_dat)
        batches.extend([r[1]] * the_dat.shape[1])

    ref = pd.concat(ref_arr, axis=1)
    ref.columns = ref_cols
    batches = pd.Series(batches, index=ref_cols)
    labels = pd.Series(ref_labels, ref_cols)
    # ref.index = ref.index.str.replace(r'.[0-9]+$', '')

    # discard Barres irrelevant samples
    # discard immortalised cell line
    # discard fibroblasts (careful, need to be selective here)

    to_discard = [
        'INSC fibroblast',
        'fetal NSC',
        'H1 NSC',
    ]
    for td in to_discard:
        the_idx = ~ref.columns.str.contains(td)
        ref = ref.loc[:, the_idx]
        batches = batches.loc[the_idx]
        labels = labels.loc[the_idx]

    dat = pd.concat((dat_home, ref), axis=1).dropna(axis=0)
    dat = dat.loc[(dat > min_val).sum(axis=1) > 6] + 1
    dat_qn = transformations.quantile_normalisation(dat)
    dat = np.log2(dat + 1)
    cc = pd.DataFrame('k', index=dat.columns, columns=['foo'])
