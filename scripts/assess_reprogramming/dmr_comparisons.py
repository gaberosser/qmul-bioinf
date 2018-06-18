from methylation import loader, dmr, process
from plotting import venn
import pandas as pd
from utils import output, setops
import numpy as np
import multiprocessing as mp
from matplotlib import pyplot as plt
import seaborn as sns
import os


def load_methylation(pids, ref_names=None, norm_method='swan', ref_name_filter=None, units='beta'):
    """
    Load and prepare the Illumina methylation data
    """
    # patient data
    obj = loader.load_by_patient(pids, norm_method=norm_method)
    anno = loader.load_illumina_methylationepic_annotation()

    # reference data
    if ref_names is not None:
        ref_obj = loader.load_reference(ref_names, norm_method=norm_method)
        if ref_name_filter is not None:
            ref_obj.filter_by_sample_name(ref_name_filter, exact=True)
        obj = loader.loader.MultipleBatchLoader([obj, ref_obj])

    me_data = obj.data.dropna()
    if units == 'm':
        me_data = process.m_from_beta(me_data)

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(me_data.index)

    anno = anno.loc[common_probes]
    # dmr.add_merged_probe_classes(anno)
    me_data = me_data.loc[common_probes]
    obj.data = me_data

    return obj, anno


def compute_dmr_clusters(anno, dmr_params):
    clusters = []
    cid = 0

    for cc in anno.CHR.unique():
        coords = anno.loc[anno.CHR == cc, 'MAPINFO'].sort_values()
        this_clust = dmr.identify_cluster(coords, dmr_params['n_min'], dmr_params['d_max'])

        for cl in this_clust.values():
            clusters.append(
                dmr.ProbeCluster(cl, anno, cluster_id=cid, chr=cc)
            )
            cid += 1
    return dmr.DmrResults(clusters=clusters, anno=anno)


def run_dmr(me_meta, me_data, dmr_clusters, str_contains_dict, type1='iPSC', type2='FB', **dmr_params):
    """
    Run DMR as type1 - type2.
    :param me_meta:
    :param me_data:
    :param dmr_clusters:
    :param str_contains_dict: Dictionary. Keys will be used to identify results. Values are 2-tuples with the strings
    used in the meta.index.str.contains() call that separates samples.
    :param type1: Refers to the `type` column in meta.
    :param type2:
    :param dmr_params:
    :return:
    """
    res = {}

    for k, (c1, c2) in str_contains_dict.items():
        this = dmr_clusters.copy()
        the_idx1 = me_meta.index.str.contains(c1) & (me_meta.loc[:, 'type'] == type1)
        the_idx2 = me_meta.index.str.contains(c2) & (me_meta.loc[:, 'type'] == type2)
        the_idx = the_idx1 | the_idx2
        the_groups = me_meta.loc[the_idx, 'type'].values
        the_samples = me_meta.index[the_idx].groupby(the_groups)
        the_samples = [the_samples[type1], the_samples[type2]]

        this.test_clusters(
            me_data,
            samples=the_samples,
            n_jobs=dmr_params['n_jobs'],
            min_median_change=dmr_params['delta_m_min'],
            method=dmr_params['dmr_test_method'],
            alpha=dmr_params['alpha'],
            **dmr_params['test_kwargs']
        )
        res[k] = this
    return dmr.DmrResultCollection(**res)


def combine_data_meta(data_arr, meta_arr, units='beta'):
    if len(data_arr) != len(meta_arr):
        raise ValueError("data_arr and meta_arr must have the same size")

    # include all probes again
    dat = pd.concat(
        data_arr,
        axis=1,
        join='inner'
    )
    meta = pd.concat(
        meta_arr,
        axis=0,
        join='outer',
        sort=True
    )
    if units.lower() == 'm':
        # convert to M values
        dat = process.m_from_beta(dat)

    # drop any infinite valued probes (should be very few)
    inft = (~np.isfinite(dat)).sum(axis=1) > 0

    if inft.any():
        dat = dat.loc[~inft]
        print "Dropped %d probes with infinite M values" % inft.sum()

    return meta, dat


if __name__ == "__main__":
    outdir = output.unique_output_dir("assess_reprogramming_dmr")
    pids = ['019', '030', '031', '050', '054']

    # these are the only two norming methods available in all data sets
    # norm_method = 'raw'
    norm_method = 'bmiq'
    # norm_method = 'pbc'

    dmr_params = {
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'alpha': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': mp.cpu_count(),
    }

    # our data
    me_obj, anno = load_methylation(pids, norm_method=norm_method)
    our_data = me_obj.data
    our_meta = me_obj.meta

    # discard unneeded samples
    our_meta = our_meta.loc[our_meta.type.isin(['iPSC', 'FB'])]
    our_data = our_data.loc[:, our_meta.index]
    our_meta.loc[:, 'batch'] = 'Our data'
    our_meta.insert(1, 'array_type', 'EPIC')

    # ref data

    # Encode EPIC data
    encode_epic_obj = loader.load_reference(
        'encode_epic',
        norm_method=norm_method,
        samples=['H7 hESC', 'GM23248', 'GM23338', 'IMR-90']
    )
    encode_epic_obj.meta.insert(1, 'array_type', 'EPIC')
    encode_epic_obj.batch_id = 'Encode EPIC'

    # Encode 450K data
    encode_450k_obj = loader.encode_450k(norm_method=norm_method, samples=['H1 hESC'])
    encode_450k_obj.meta.insert(1, 'array_type', '450K')

    # E-MTAB-6194
    e6194_obj = loader.load_reference(
        'E-MTAB-6194',
        norm_method=norm_method,
    )
    discard_ix = e6194_obj.meta.cell_line.isin([
        'NA07057',
        'HCT116',
        'HEL46.11',
    ])
    e6194_obj.meta = e6194_obj.meta.loc[~discard_ix]
    e6194_obj.data = e6194_obj.data[e6194_obj.meta.index]
    e6194_obj.meta.insert(1, 'array_type', 'EPIC')

    refs = [
        encode_epic_obj,
        encode_450k_obj,
        e6194_obj
    ]

    # HipSci data
    hip_epic_meta, hip_epic_data = loader.hipsci(norm_method=norm_method, n_sample=12, array_type='epic')

    meta, dat_m = combine_data_meta(
        (our_data, hip_epic_data, e6194_obj.data, encode_epic_obj.data),
        (our_meta, hip_epic_meta, e6194_obj.meta, encode_epic_obj.meta),
        units='m'
    )

    this_anno = anno.loc[dat_m.index]
    dmr_clusters = compute_dmr_clusters(this_anno, dmr_params)

    ipsc_ref_names_6194 = ['HEL139', 'HEL140', 'HEL141']
    esc_ref_names = ['H7', 'H9']

    # 1. iPSC vs FB

    # 1a) Our iPSC vs matched parental FB

    comparisons = {}
    for pid in pids:
        comparisons[pid] = ("DURA%s" % pid, "DURA%s" % pid)

    dmr_res_our_ipsc_vs_our_fb = run_dmr(
        meta,
        dat_m,
        dmr_clusters,
        comparisons,
        type1='iPSC',
        type2='FB',
        **dmr_params
    )
    dmr_res_our_ipsc_vs_our_fb.to_pickle(
        os.path.join(outdir, "our_ipsc_vs_our_fb.pkl"),
        include_annotation=False
    )

    # 1b) HipSci iPSC vs our FB

    comparisons = {}
    for pid in pids:
        for hid in hip_epic_meta.index:
            comparisons["%s-%s" % (pid, hid)] = (hid, "DURA%s" % pid)

    dmr_res_hipsci_vs_our_fb = run_dmr(
        meta,
        dat_m,
        dmr_clusters,
        comparisons,
        type1='iPSC',
        type2='FB',
        **dmr_params
    )
    dmr_res_hipsci_vs_our_fb.to_pickle(
        os.path.join(outdir, "hipsci_vs_our_fb.pkl"),
        include_annotation=False
    )

    # 1c) E-MTAB-6194 iPSC vs our FB

    comparisons = {}
    for pid in pids:
        for r in ipsc_ref_names_6194:
            comparisons["%s-%s" % (r, pid)] = (r, "DURA%s" % pid)

    dmr_res_e6194_ipsc_vs_our_fb = run_dmr(
        meta,
        dat_m,
        dmr_clusters,
        comparisons,
        type1='iPSC',
        type2='FB',
        **dmr_params
    )
    dmr_res_e6194_ipsc_vs_our_fb.to_pickle(
        os.path.join(outdir, "e6194_ipsc_vs_our_fb.pkl"),
        include_annotation=False
    )

    # 2. iPSC vs ESC

    # 2a) Our iPSC vs 2 x EPIC reference ESC (Encode, E-MTAB-6194)

    comparisons = {}
    for pid in pids:
        for r in esc_ref_names:
            comparisons["%s-%s" % (pid, r)] = ("DURA%s" % pid, r)

    dmr_res_our_ipsc_vs_esc = run_dmr(
        meta,
        dat_m,
        dmr_clusters,
        comparisons,
        type1='iPSC',
        type2='ESC',
        **dmr_params
    )
    dmr_res_our_ipsc_vs_esc.to_pickle(
        os.path.join(outdir, "our_ipsc_vs_esc.pkl"),
        include_annotation=False
    )

    # 2b) HipSci iPSC vs 2 x EPIC reference ESC (Encode, E-MTAB-6194)

    comparisons = {}
    for hid in hip_epic_meta.index:
        for r in esc_ref_names:
            comparisons["%s-%s" % (hid, r)] = (hid, r)

    dmr_res_hipsci_vs_esc = run_dmr(
        meta,
        dat_m,
        dmr_clusters,
        comparisons,
        type1='iPSC',
        type2='ESC',
        **dmr_params
    )
    dmr_res_hipsci_vs_esc.to_pickle(
        os.path.join(outdir, "hipsci_vs_esc.pkl"),
        include_annotation=False
    )

    # 2c) E-MTAB-6194 iPSC vs 2 x EPIC reference ESC (Encode, E-MTAB-6194)

    comparisons = {}
    for hid in ipsc_ref_names_6194:
        for r in esc_ref_names:
            comparisons["%s-%s" % (hid, r)] = (hid, r)

    dmr_res_e6914_ipsc_vs_esc = run_dmr(
        meta,
        dat_m,
        dmr_clusters,
        comparisons,
        type1='iPSC',
        type2='ESC',
        **dmr_params
    )
    dmr_res_e6914_ipsc_vs_esc.to_pickle(
        os.path.join(outdir, "e6194_ipsc_vs_esc.pkl"),
        include_annotation=False
    )

    # 3. iPSC vs iPSC
    # for these comparisons, we need to modify the types to distinguish the batches
    meta.loc[meta.batch.str.contains('HipSci').fillna(False), 'type'] = 'iPSC_HipSci'
    meta.loc[meta.index.str.contains('HEL1'), 'type'] = 'iPSC_E6194'

    # 3a) Our iPSC vs HipSci iPSC

    comparisons = {}
    for hid in hip_epic_meta.index:
        for pid in pids:
            comparisons["%s-%s" % (pid, hid)] = ("DURA%s" % pid, hid)
    dmr_res_hipsci_vs_our_ipsc = run_dmr(
        meta,
        dat_m,
        dmr_clusters,
        comparisons,
        type1='iPSC',
        type2='iPSC_HipSci',
        **dmr_params
    )
    dmr_res_hipsci_vs_our_ipsc.to_pickle(
        os.path.join(outdir, "hipsci_vs_our_ipsc.pkl"),
        include_annotation=False
    )

    # 3b) Our iPSC vs E-MTAB-6194 iPSC

    comparisons = {}
    for hid in ipsc_ref_names_6194:
        for pid in pids:
            comparisons["%s-%s" % (pid, hid)] = ("DURA%s" % pid, hid)
    dmr_res_e6194_vs_our_ipsc = run_dmr(
        meta,
        dat_m,
        dmr_clusters,
        comparisons,
        type1='iPSC',
        type2='iPSC_E6194',
        **dmr_params
    )
    dmr_res_e6194_vs_our_ipsc.to_pickle(
        os.path.join(outdir, "e6194_ipsc_vs_our_ipsc.pkl"),
        include_annotation=False
    )

    # 3c) HipSci vs E-MTAB-6194 iPSC

    comparisons = {}
    for hid in hip_epic_meta.index:
        for r in ipsc_ref_names_6194:
            comparisons["%s-%s" % (hid, r)] = (hid, r)
    dmr_res_hipsci_vs_e6194 = run_dmr(
        meta,
        dat_m,
        dmr_clusters,
        comparisons,
        type1='iPSC_HipSci',
        type2='iPSC_E6194',
        **dmr_params
    )
    dmr_res_hipsci_vs_e6194.to_pickle(
        os.path.join(outdir, "hipsci_vs_e6194_ipsc.pkl"),
        include_annotation=False
    )

    # Analyse these results
    # TODO: refactor from here

    # # for each PID, define the core DMRs (shared by both ref comparisons)
    # core_dmr_our_ipsc_ref_esc = {}
    # for pid in pids:
    #     this_sets = []
    #     for r in ['H9', 'H7']:
    #         this_res = dmr_res_s2['%s-%s' % (pid, r)]
    #         this_sets.append(this_res.results_significant.keys())
    #     core_cids = setops.reduce_intersection(*this_sets)
    #     tbls = []
    #     for r in ['H9', 'H7']:
    #         this_res = dmr_res_s2['%s-%s' % (pid, r)]
    #         ## FIXME: this ugly hack is necessary if classes are not defined (make it not so) to run to_table
    #         this_res._classes = []
    #         this_tbl = this_res.to_table(include='significant', skip_geneless=False).loc[core_cids]
    #         this_tbl.columns = ["%s_%s" % (t, r) for t in this_tbl.columns]
    #         tbls.append(this_tbl)
    #     this_comb = pd.concat(tbls, axis=1)
    #     for col in ['chr', 'genes', 'median_1']:
    #         this_comb.insert(0, col, this_comb['%s_H9' % col])
    #         this_comb.drop('%s_H9' % col, axis=1, inplace=True)
    #         this_comb.drop('%s_H7' % col, axis=1, inplace=True)
    #     core_dmr_our_ipsc_ref_esc[pid] = this_comb
    #
    # # for each PID, plot the Venn diag
    # fig, axs = plt.subplots(nrows=len(pids), figsize=(3, 11))
    # for i, pid in enumerate(pids):
    #     ax = axs[i]
    #     this_sets = []
    #     for r in ['H9', 'H7']:
    #         this_res = dmr_res_s2['%s-%s' % (pid, r)]
    #         this_sets.append(this_res.results_significant.keys())
    #     venn.venn_diagram(*this_sets, set_labels=['H9', 'H7'], ax=ax)
    #     ax.set_title(pid)
    # fig.tight_layout()
    # fig.savefig(os.path.join(outdir, "venn_overlap_dmrs_our_ipsc_vs_ref.png"), dpi=200)
    #
    # # starting with core DMRs, split into hyper and hypo and check agreement between refs
    # our_ipsc_ref_esc_direction = {}
    # for pid, v in core_dmr_our_ipsc_ref_esc.items():
    #     disagree_ix = np.sign(v.median_delta_H7) != np.sign(v.median_delta_H9)
    #     n_disagree = disagree_ix.sum()
    #     print "Patient %s. Of the %d DMRs (iPSC - ref. ESC), %d do not agree in direction." % (
    #         pid, v.shape[0], n_disagree
    #     )
    #     this_med_delta = v.loc[~disagree_ix].median_delta_H9
    #     our_ipsc_ref_esc_direction[pid] = {
    #         'hypo': (this_med_delta < 0).sum(),
    #         'hyper': (this_med_delta > 0).sum(),
    #     }
    # our_ipsc_ref_esc_direction = pd.DataFrame.from_dict(our_ipsc_ref_esc_direction)
    # ax = our_ipsc_ref_esc_direction.transpose().plot.bar()
    # ax.figure.savefig(os.path.join(outdir, "our_ipsc_esc_ref_dmr_direction.png"), dpi=200)
    #
    #
    # # DMRs: iPSC (HipSci) vs ESC
    # # Let's use H7 and H9 by Zimmerlin for this purpose
    # res_hipsci_esc = {}
    # suff = ' hESC_Zimmerlin et al.'
    # hip_ids = hip_meta.index[:12]
    # for pid in hip_ids:
    #     for r in ['H9', 'H7']:
    #         this = dmr_clusters.copy()
    #         the_idx1 = meta.index.str.contains(pid) & (meta.loc[:, 'type'] == 'iPSC')
    #         the_idx2 = meta.index == (r + suff)
    #         the_idx = the_idx1 | the_idx2
    #         the_groups = meta.loc[the_idx, 'type'].values
    #         the_samples = meta.index[the_idx].groupby(the_groups)
    #         the_samples = [the_samples['iPSC'], the_samples['ESC']]
    #
    #         this.test_clusters(
    #             dat_m,
    #             samples=the_samples,
    #             n_jobs=dmr_params['n_jobs'],
    #             min_median_change=dmr_params['delta_m_min'],
    #             method=dmr_params['dmr_test_method'],
    #             alpha=dmr_params['alpha'],
    #             **dmr_params['test_kwargs']
    #         )
    #         res_hipsci_esc["%s-%s" % (pid, r)] = this
    # dmr_res_hipsci_esc = dmr.DmrResultCollection(**res_hipsci_esc)
    #
    # # for each PID, define the core DMRs (shared by both ref comparisons)
    # core_dmr_hipsci_ref_esc = {}
    # for pid in hip_ids:
    #     this_sets = []
    #     for r in ['H9', 'H7']:
    #         this_res = dmr_res_hipsci_esc['%s-%s' % (pid, r)]
    #         this_sets.append(this_res.results_significant.keys())
    #     core_cids = setops.reduce_intersection(*this_sets)
    #     tbls = []
    #     for r in ['H9', 'H7']:
    #         this_res = dmr_res_hipsci_esc['%s-%s' % (pid, r)]
    #         ## FIXME: this ugly hack is necessary if classes are not defined (make it not so) to run to_table
    #         this_res._classes = []
    #         this_tbl = this_res.to_table(include='significant', skip_geneless=False).loc[core_cids]
    #         this_tbl.columns = ["%s_%s" % (t, r) for t in this_tbl.columns]
    #         tbls.append(this_tbl)
    #     this_comb = pd.concat(tbls, axis=1)
    #     for col in ['chr', 'genes', 'median_1']:
    #         this_comb.insert(0, col, this_comb['%s_H9' % col])
    #         this_comb.drop('%s_H9' % col, axis=1, inplace=True)
    #         this_comb.drop('%s_H7' % col, axis=1, inplace=True)
    #     core_dmr_hipsci_ref_esc[pid] = this_comb
    #
    # # for each PID, plot the Venn diag
    # fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(9, 8))
    # for i, pid in enumerate(hip_ids):
    #     ax = axs.flat[i]
    #     this_sets = []
    #     for r in ['H9', 'H7']:
    #         this_res = dmr_res_hipsci_esc['%s-%s' % (pid, r)]
    #         this_sets.append(this_res.results_significant.keys())
    #     venn.venn_diagram(*this_sets, set_labels=['H9', 'H7'], ax=ax)
    #     ax.set_title(pid)
    # fig.tight_layout()
    # fig.savefig(os.path.join(outdir, "venn_overlap_dmrs_hipsci_ipsc_vs_ref.png"), dpi=200)
    #
    # # starting with core DMRs, split into hyper and hypo and check agreement between refs
    # hipsci_ref_esc_direction = {}
    # for pid, v in core_dmr_hipsci_ref_esc.items():
    #     disagree_ix = np.sign(v.median_delta_H7) != np.sign(v.median_delta_H9)
    #     n_disagree = disagree_ix.sum()
    #     print "HipSci %s. Of the %d DMRs (iPSC - ref. ESC), %d do not agree in direction." % (
    #         pid, v.shape[0], n_disagree
    #     )
    #     this_med_delta = v.loc[~disagree_ix].median_delta_H9
    #     hipsci_ref_esc_direction[pid] = {
    #         'hypo': (this_med_delta < 0).sum(),
    #         'hyper': (this_med_delta > 0).sum(),
    #     }
    # hipsci_ref_esc_direction = pd.DataFrame.from_dict(hipsci_ref_esc_direction)
    # ax = hipsci_ref_esc_direction.transpose().plot.bar()
    # ax.figure.tight_layout()
    # ax.figure.savefig(os.path.join(outdir, "hipsci_esc_ref_dmr_direction.png"), dpi=200)
