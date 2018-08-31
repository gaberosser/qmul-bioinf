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


if __name__ == '__main__':
    pids = ['017', '018', '019', '030', '031', '049', '050', '054', '061', '026', '052']
    min_val = 1
    n_above_min = 3
    eps = 0.01  # offset to use when applying log transform

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
    # to_aggr = [
    #     (r'Fibroblasts_control_rep[123]', 'FB control'),
    #     (r'H1-hESC rep (1_1|2_1|3|4)', 'H1 hESC'),
    #     (r'H1-hESC rep [12]', 'H1 hESC'),
    #     (r'H1-[12]', 'H1 hESC'),
    #     (r'H1-hESC(|_2)$', 'H1 hESC'),
    #     (r'H7-hESC rep [12]', 'H7 hESC'),
    #     (r'hESCs_control_rep[123]', 'CSES9 hESC'),
    # ]

    nsc_ref_dict = {
        # 'E-MTAB-3867': {'batch': 'Caren et al.', 'strandedness': '?'},
        'GSE61794': {'batch': 'Duan et al.', 'strandedness': '?'},
        # 'GSE64882': {'batch': 'Shahbazi et al.', 'strandedness': '?'},
        # 'encode_roadmap/ENCSR244ISQ': {'batch': 'ENCODE Gingeras', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR291IZK': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR572EET': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        # 'encode_roadmap/ENCSR977XUX': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
    }

    # to_aggr_nsc = [
    #     (r'H9_NSC_[12]', 'H9 NSC'),
    #     # (r'Pollard NSC [12]', 'Fetal NSC'),
    # ]

    outdir = output.unique_output_dir("assess_reprogramming_de")

    load_kwds = {
        'source': 'star',
        'alignment_subdir': SetMe,
        'strandedness': SetMe,
    }

    # our data (everything)
    obj = loader.load_by_patient(pids, source='star')

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