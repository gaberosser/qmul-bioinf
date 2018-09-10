from rnaseq import loader, differential_expression, filter
from plotting import common, venn
import pandas as pd
import numpy as np
import os
from utils import output
from matplotlib import pyplot as plt
import seaborn as sns
from copy import copy


class SetMe(object):
    """
    This is useful to determine which fields need to be set
    """


def de_grouped_dispersion(dat, groups, comparisons, min_cpm=1., **de_params):
    cpm = dat.divide(dat.sum(), axis=1) * 1e6
    keep = (cpm > min_cpm).sum(axis=1) > 0
    dat = dat.loc[keep]

    res = differential_expression.run_multiple_de(
        the_data=dat,
        the_groups=groups,
        comparisons=comparisons,
        **de_params
    )

    for the_comparison, this_res in res.items():
        try:
            the_genes = this_res.index
            the_groups = groups[groups.isin(the_comparison)]
            the_cpm = cpm.loc[the_genes, the_groups.index]
            keep = filter.filter_cpm_by_group(
                the_cpm,
                groups,
                min_cpm=min_cpm
            )
            res[the_comparison] = this_res.loc[keep]
        except Exception as exc:
            print repr(exc)

    return res


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
    min_cpm = 1.
    n_above_min = 3
    eps = 0.01  # offset to use when applying log transform
    n_hipsci = 12

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'QLGLM'
    }

    ref_dict = {
        'encode_roadmap/ENCSR000EYP': {'batch': 'ENCODE Wold', 'strandedness': 'u'},
        # 'encode_roadmap/ENCSR670WQY': {'batch': 'ENCODE Ecker', 'strandedness': 'r'},
        'GSE62772': {'batch': 'Cacchiarelli et al.', 'strandedness': 'u'},
        'GSE97265': {'batch': 'Kogut et al.', 'strandedness': 'r'},
    }

    # discard irrelevant samples
    to_discard = [
        '_rapamycin_',
        'HFF_PD16', 'HFF_PD46', 'HFF_PD64', 'HFF_PD74',
        'BJ_OLD',
        'IMR90_O',
        'MRC_5_PD52', 'MRC_5_PD62', 'MRC_5_PD72',
        'WI_38_O',
        '-EPS-',
        'F50S', 'I50S', 'FN1',
        'DOX', 'hiF_', 'hiF-T_', 'BJ_P12', '_MTG_', '_VTN_', 'hIPSC_P15_Rep1'
    ]

    # TODO: update this appropriately
    to_aggr = [
        (r'H1-hESC rep (1_1|2_1|3|4)', 'H1 hESC'),
        (r'H1-hESC rep [12]', 'H1 hESC'),
        (r'H1-[12]', 'H1 hESC'),
        (r'H1-hESC(|_2)$', 'H1 hESC'),
        (r'H7-hESC rep [12]', 'H7 hESC'),
        (r'hESCs_control_rep[123]', 'CSES9 hESC'),
        ('IN2-', 'IN2'),
        ('I50-', 'I50'),
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
    # ix = obj.meta.type.isin(['iPSC', 'FB'])
    # obj.filter_samples(ix)

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
    for td in to_discard:
        the_idx = ~ref_obj.data.columns.str.contains(td)
        ref_obj.filter_samples(the_idx)

    # fill in missing cell types
    ref_obj.meta.loc[ref_obj.meta.type.isnull(), 'type'] = ref_obj.meta.loc[ref_obj.meta.type.isnull(), 'cell type']
    ref_obj.rename_with_attributes(existing_attr='batch')

    for srch, repl in to_aggr:
        ref_obj.aggregate_by_pattern(srch, repl)

    # assemble the data we need for comparison:
    # - ESC (2 studies)
    # - iPSC (ours, GSE97265)
    # - FB (ours, GSE97265)
    ref_labels = ['cacchiarelli', 'encode']

    ix = obj.meta.type.isin(['iPSC', 'FB'])
    dat1 = obj.data.loc[:, ix]

    ix = ref_obj.meta.batch.str.contains('Kogut') & ref_obj.meta.type.isin(['iPSC', 'FB'])
    dat2 = ref_obj.data.loc[:, ix]

    ix = ref_obj.meta.index.str.contains('Cacchiarelli') & (ref_obj.meta.type == 'ESC')
    dat_esc_1 = ref_obj.data.loc[:, ix]

    ix = ref_obj.meta.index.str.contains('ENCODE') & (ref_obj.meta.type == 'ESC')
    dat_esc_2 = ref_obj.data.loc[:, ix]

    the_dat = pd.concat((dat1, dat2, dat_esc_1, dat_esc_2), axis=1)
    the_groups = pd.Series(index=the_dat.columns)
    the_comparisons = {}

    for pid in pids:
        k_fb = 'FB_%s_ours' % pid
        k_ipsc = 'iPSC_%s_ours' % pid
        ix_fb = the_groups.index.str.contains('DURA%s_FB' % pid)
        ix_ipsc = the_groups.index.str.contains('DURA%s_IPSC' % pid)
        the_groups.loc[ix_fb] = k_fb
        the_groups.loc[ix_ipsc] = k_ipsc
        if (ix_fb.sum() > 0) & (ix_ipsc.sum() > 0):
            the_comparisons[(k_ipsc, k_fb)] = "%s - %s" % (k_ipsc, k_fb)
            for r in ref_labels:
                the_comparisons[(k_ipsc, 'ESC_%s' % r)] = "%s - ESC_%s" % (k_ipsc, r)

    for s in ['N2', '50']:
        k_fb = 'FB_%s_kogut' % s
        k_ipsc = 'iPSC_%s_kogut' % s
        ix_fb = the_groups.index.str.contains('F%s' % s) & the_groups.index.str.contains('Kogut')
        ix_ipsc = the_groups.index.str.contains('I%s' % s) & the_groups.index.str.contains('Kogut')
        the_groups.loc[ix_fb] = k_fb
        the_groups.loc[ix_ipsc] = k_ipsc
        if (ix_fb.sum() > 0) & (ix_ipsc.sum() > 0):
            the_comparisons[(k_ipsc, k_fb)] = "%s - %s" % (k_ipsc, k_fb)
            for r in ref_labels:
                the_comparisons[(k_ipsc, 'ESC_%s' % r)] = "%s - ESC_%s" % (k_ipsc, r)

    the_groups.loc[the_groups.index.str.contains('Cacchiarelli')] = 'ESC_cacchiarelli'
    the_groups.loc[the_groups.index.str.contains('ENCODE')] = 'ESC_encode'

    de_res_full = de_grouped_dispersion(
        the_dat,
        the_groups,
        the_comparisons,
        min_cpm=min_cpm,
        return_full=True,
        **de_params
    )
    # rename the keys to simplify
    # de_res_full_s1 = dict([(pid, de_res_full[("GBM%s" % pid, "iNSC%s" % pid)]) for pid in pids])

    # extract only significant DE genes
    de_res_sign = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full.items()])

    # Venn diagrams (two ESC studies)
    s1 = ['_'.join(t) for t in zip(pids, ['ours'] * len(pids))]
    s2 = ['_'.join(t) for t in zip(['N2', '50'], ['kogut'] * 2)]
    fig, axs = plt.subplots(ncols=3, nrows=3)
    i = 0
    for s in s1 + s2:
        this_arr = []
        for r in ref_labels:
            k = ('iPSC_%s' % s, 'ESC_%s' % r)
            if k in de_res_sign:
                this_arr.append(de_res_sign[k])
        if len(this_arr):
            print "Found comparison %s" % s
            venn.venn_diagram(*[t.index for t in this_arr], ax=axs.flat[i])
            axs.flat[i].set_title(s)
            i += 1
        else:
            print "No comparison %s" % s
