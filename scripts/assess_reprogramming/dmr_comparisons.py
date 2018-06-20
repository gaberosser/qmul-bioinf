from methylation import loader, dmr, process
from plotting import venn
import pandas as pd
from utils import output, setops
import numpy as np
import multiprocessing as mp
from matplotlib import pyplot as plt
import seaborn as sns
import os
from scipy import stats
from settings import OUTPUT_DIR


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


def run_dmr_set(fn, meta, dat_m, dmr_clusters, anno, comparisons, type1, type2, dmr_params):
    if not os.path.exists(fn):
        res = run_dmr(
            meta,
            dat_m,
            dmr_clusters,
            comparisons,
            type1='iPSC',
            type2='FB',
            **dmr_params
        )
        res.to_pickle(
            fn,
            include_annotation=False
        )
    else:
        res = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    return res


def get_dmr_number_direction(res_obj):
    out = {}
    for k in res_obj.keys():
        this_res = pd.DataFrame(res_obj[k].results_significant)
        ix = this_res.loc['median_change'] < 0
        out[k] = {
            'Hypermethylated': (~ix).sum(),
            'Hypomethylated': ix.sum(),
        }
    return out


def get_dmr_cid_direction(res_obj):
    out = {}
    for k in res_obj.keys():
        this_res = pd.DataFrame(res_obj[k].results_significant)
        ix = this_res.loc['median_change'] < 0
        out[k] = {
            'Hypermethylated': this_res.loc[:, ~ix].columns,
            'Hypomethylated': this_res.loc[:, ix].columns,
        }
    return out


if __name__ == "__main__":

    outdir = output.unique_output_dir("assess_reprogramming_dmr")
    indir = os.path.join(OUTPUT_DIR, "assess_reprogramming_dmr")
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
    ipsc_ref_names_6194_n1 = ['HEL139.2_p17', 'HEL139.5_p14', 'HEL139.8_p13', 'HEL140.1_p12', 'HEL141.1_p11']
    fb_ref_name_6194_n1 = '27_HFF_p7'
    esc_ref_names = ['H7', 'H9']

    # 1. iPSC vs FB

    # 1a-i) Our iPSC vs our (matched) FB
    fn = os.path.join(indir, "our_ipsc_vs_our_fb.pkl")
    comparisons = {}
    for pid in pids:
        comparisons[pid] = ("DURA%s" % pid, "DURA%s" % pid)
    dmr_res_our_ipsc_vs_our_fb = run_dmr_set(
        fn, meta, dat_m, dmr_clusters, anno, comparisons, 'iPSC', 'FB', dmr_params
    )

    # 1a-ii) Our iPSC vs our FB (all possible comparisons)
    fn = os.path.join(indir, "our_ipsc_vs_our_fb_all.pkl")
    comparisons = {}
    for pid in pids:
        for pid2 in pids:
            comparisons["%s-%s" % (pid, pid2)] = ("DURA%s" % pid, "DURA%s" % pid2)
    dmr_res_our_ipsc_vs_our_fb_all = run_dmr_set(
        fn, meta, dat_m, dmr_clusters, anno, comparisons, 'iPSC', 'FB', dmr_params
    )

    # 1b) HipSci iPSC vs our FB
    fn = os.path.join(indir, "hipsci_vs_our_fb.pkl")
    comparisons = {}
    for pid in pids:
        for hid in hip_epic_meta.index:
            comparisons["%s-%s" % (pid, hid)] = (hid, "DURA%s" % pid)
    dmr_res_hipsci_vs_our_fb = run_dmr_set(
        fn, meta, dat_m, dmr_clusters, anno, comparisons, 'iPSC', 'FB', dmr_params
    )

    # 1c-i) E-MTAB-6194 iPSC (all replicates) vs our FB
    fn = os.path.join(indir, "e6194_ipsc_vs_our_fb.pkl")
    comparisons = {}
    for pid in pids:
        for r in ipsc_ref_names_6194:
            comparisons["%s-%s" % (r, pid)] = (r, "DURA%s" % pid)

    dmr_res_e6194_ipsc_vs_our_fb = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'FB',
        dmr_params
    )

    # 1c-ii) E-MTAB-6194 iPSC (n=1) vs our FB
    fn = os.path.join(indir, "e6194_ipsc_n1_vs_our_fb.pkl")
    comparisons = {}
    for pid in pids:
        for r in ipsc_ref_names_6194_n1:
            comparisons["%s-%s" % (r, pid)] = (r, "DURA%s" % pid)

    dmr_res_e6194_ipsc_n1_vs_our_fb = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'FB',
        dmr_params
    )

    # 1d-i) E-MTAB-6194 (all replicates) vs matched FB (only those iPSC derived from CRL2429)
    fn = os.path.join(indir, "e6194_ipsc_vs_e6194_fb.pkl")
    comparisons = {}
    for r in ['HEL140', 'HEL141']:
        comparisons[r] = (r, "HFF_p7")

    dmr_res_e6194_ipsc_vs_e6194_fb = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'FB',
        dmr_params
    )

    # 1d-ii) E-MTAB-6194 (n=1) vs matched FB (n=1)
    fn = os.path.join(indir, "e6194_ipsc_n1_vs_e6194_fb_n1.pkl")
    comparisons = {}
    for r in ipsc_ref_names_6194_n1:
        comparisons[r] = (r, fb_ref_name_6194_n1)

    dmr_res_e6194_ipsc_n1_vs_e6194_fb_n1 = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'FB',
        dmr_params
    )

    # 1e-i) Our iPSC vs E-MTAB-6194 FB (CRL2429, all replicates)
    fn = os.path.join(indir, "our_ipsc_vs_e6194_fb.pkl")
    comparisons = {}
    for pid in pids:
        comparisons["%s-HFF" % pid] = ("DURA%s" % pid, "HFF_p7")

    dmr_res_our_ipsc_vs_e6194_fb = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'FB',
        dmr_params
    )

    # 1e-ii) Our iPSC vs E-MTAB-6194 FB (CRL2429, n=1)
    fn = os.path.join(indir, "our_ipsc_vs_e6194_fb_n1.pkl")
    comparisons = {}
    for pid in pids:
        comparisons["%s-HFF" % pid] = ("DURA%s" % pid, fb_ref_name_6194_n1)

    dmr_res_our_ipsc_vs_e6194_fb_n1 = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'FB',
        dmr_params
    )

    # 1f-i) HipSci iPSC vs E-MTAB-6194 FB (CRL2429, all replicates)
    fn = os.path.join(indir, "hipsci_vs_e6194_fb.pkl")
    comparisons = {}
    for hid in hip_epic_meta.index:
        comparisons["%s-HFF" % hid] = ("%s" % hid, "HFF_p7")

    dmr_res_hipsci_vs_e6194_fb= run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'FB',
        dmr_params
    )

    # 1f-ii) HipSci iPSC vs E-MTAB-6194 FB (CRL2429, n=1)
    fn = os.path.join(indir, "hipsci_vs_e6194_fb_n1.pkl")
    comparisons = {}
    for hid in hip_epic_meta.index:
        comparisons["%s-HFF" % hid] = ("%s" % hid, fb_ref_name_6194_n1)

    dmr_res_hipsci_vs_e6194_fb_n1 = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'FB',
        dmr_params
    )

    # 2. iPSC vs ESC

    # 2a) Our iPSC vs 2 x EPIC reference ESC (Encode, E-MTAB-6194)
    fn = os.path.join(indir, "our_ipsc_vs_esc.pkl")
    comparisons = {}
    for pid in pids:
        for r in esc_ref_names:
            comparisons["%s-%s" % (pid, r)] = ("DURA%s" % pid, r)

    dmr_res_our_ipsc_vs_esc = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'ESC',
        dmr_params
    )

    # 2b) HipSci iPSC vs 2 x EPIC reference ESC (Encode, E-MTAB-6194)
    fn = os.path.join(indir, "hipsci_vs_esc.pkl")
    comparisons = {}
    for hid in hip_epic_meta.index:
        for r in esc_ref_names:
            comparisons["%s-%s" % (hid, r)] = (hid, r)

    dmr_res_hipsci_vs_esc = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'ESC',
        dmr_params
    )

    # 2c) E-MTAB-6194 iPSC vs 2 x EPIC reference ESC (Encode, E-MTAB-6194)
    fn = os.path.join(indir, "e6194_ipsc_vs_esc.pkl")
    comparisons = {}
    for hid in ipsc_ref_names_6194:
        for r in esc_ref_names:
            comparisons["%s-%s" % (hid, r)] = (hid, r)

    dmr_res_e6914_ipsc_vs_esc = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'ESC',
        dmr_params
    )

    # 3. iPSC vs iPSC
    # for these comparisons, we need to modify the types to distinguish the batches
    meta.loc[meta.batch.str.contains('HipSci').fillna(False), 'type'] = 'iPSC_HipSci'
    meta.loc[meta.index.str.contains('HEL1'), 'type'] = 'iPSC_E6194'

    # 3a) Our iPSC vs HipSci iPSC
    fn = os.path.join(indir, "hipsci_vs_our_ipsc.pkl")
    comparisons = {}
    for hid in hip_epic_meta.index:
        for pid in pids:
            comparisons["%s-%s" % (pid, hid)] = ("DURA%s" % pid, hid)
    dmr_res_hipsci_vs_our_ipsc = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'iPSC_HipSci',
        dmr_params
    )

    # 3b) Our iPSC vs E-MTAB-6194 iPSC
    fn = os.path.join(indir, "e6194_ipsc_vs_our_ipsc.pkl")
    comparisons = {}
    for hid in ipsc_ref_names_6194:
        for pid in pids:
            comparisons["%s-%s" % (pid, hid)] = ("DURA%s" % pid, hid)
    dmr_res_e6194_vs_our_ipsc = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC',
        'iPSC_E6194',
        dmr_params
    )

    # 3c) HipSci vs E-MTAB-6194 iPSC
    fn = os.path.join(indir, "hipsci_vs_e6194_ipsc.pkl")
    comparisons = {}
    for hid in hip_epic_meta.index:
        for r in ipsc_ref_names_6194:
            comparisons["%s-%s" % (hid, r)] = (hid, r)
    dmr_res_hipsci_vs_e6194_ipsc = run_dmr_set(
        fn,
        meta,
        dat_m,
        dmr_clusters,
        anno,
        comparisons,
        'iPSC_HipSci',
        'iPSC_E6194',
        dmr_params
    )

    # Analyse these results

    colours = {
        'Hypermethylated': '#e09191',  # red
        'Hypomethylated': '#91e097',  # green
    }

    colour_by_pid = {
        '019': '#7fc97f',
        '030': '#beaed4',
        '031': '#fdc086',
        '050': '#ffff99',
        '054': '#386cb0',
    }

    def jitter_points(x, spacing_buffer=1, jitter_step=0.02):
        curr = 0
        prev = -1e9
        out = np.zeros_like(x, dtype=float)
        for i in np.argsort(x):
            if np.abs(x[i] - prev) > spacing_buffer:
                curr = max(0., curr - jitter_step)
            else:
                curr += jitter_step
            out[i] = curr
            prev = x[i]
        return out



    # Our iPSC vs our FB: is there any obvious increased similarity between matched and unmatached comparisons?

    t = pd.DataFrame(get_dmr_number_direction(dmr_res_our_ipsc_vs_our_fb_all)).transpose()

    matched = {}
    unmatched = {}

    for pid in pids:
        for pid2 in pids:
            if pid2 == pid:
                the_dict = matched
                loc = 0.
                scale = 0.
            else:
                the_dict = unmatched
                loc = 1.
                scale = 0.05

            y = t.loc["%s-%s" % (pid, pid2), 'Hypomethylated']
            the_dict.setdefault('hypo', {}).setdefault(pid, []).append(y)
            y = t.loc["%s-%s" % (pid, pid2), 'Hypermethylated']
            the_dict.setdefault('hyper', {}).setdefault(pid, []).append(y)

    spacings = {
        'hypo': 2,
        'hyper': 7,
    }
    fig, axs = plt.subplots(1, 2, sharex=True)
    y = {}
    c = {}
    for pid in pids:
        for i, k in enumerate(['hypo', 'hyper']):
            ax = axs[i]
            this_y = unmatched[k][pid]
            y.setdefault(k, []).extend(this_y)
            this_c = [colour_by_pid[pid]] * len(this_y)
            c.setdefault(k, []).extend(this_c)
            ax.scatter(
                [0.],
                matched[k][pid],
                marker='o',
                facecolor=colour_by_pid[pid],
                s=40,
                edgecolor='k',
                linewidths=1.0,
                label=pid
            )

    for i, k in enumerate(['hypo', 'hyper']):
        ax = axs[i]
        this_y = y[k]
        this_x = jitter_points(this_y, spacing_buffer=spacings[k], jitter_step=0.05) + 1.
        this_c = c[k]
        ax.scatter(
            this_x,
            this_y,
            marker='o',
            facecolors=this_c,
            s=40,
            edgecolor='k',
            linewidths=1.0,
        )

    for ax in axs:
        ax.set_xlim([-.5, 1.5])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Matched', 'Cross-compared'])
        ax.set_ylim([0, ax.get_ylim()[1]])

    axs[0].set_title("Hypomethylated DMRs")
    axs[0].set_ylabel("Number of DMRs")
    axs[1].set_title("Hypermethylated DMRs")
    axs[1].legend(loc='lower right', facecolor='w', framealpha=0.6, frameon=True)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "our_ipsc_vs_our_fb_number_dmrs.png"), dpi=200)

    # confidence intervals
    alpha = 0.95
    for pid in pids:
        for k in ['hypo', 'hyper']:
            arr = unmatched[k][pid]
            ci = stats.norm.interval(alpha, loc=np.mean(arr), scale=np.std(arr, ddof=1))
            if ci[0] <= matched[k][pid][0] <= ci[1]:
                pass
            else:
                print "%s %s %s -> %d" % (pid, k, str(ci), matched[k][pid][0])
                print "Patient %s, %s DMRs, reject H0 at %.1f%% level." % (pid, k, alpha * 100)


    # iPSC vs FB: numbers and direction

    # Nesting format: iNSC then FB

    # to_plot = {
    #     'Our data': {'Our data': dmr_res_our_ipsc_vs_our_fb, 'E-MTAB-6194': dmr_res_our_ipsc_vs_e6194_fb},
    #     'E-MTAB-6194': {'Our data': dmr_res_e6194_ipsc_vs_our_fb, 'E-MTAB-6194': dmr_res_e6194_ipsc_vs_e6194_fb},
    #     'HipSci': {'Our data': dmr_res_hipsci_vs_our_fb, 'E-MTAB-6194': dmr_res_hipsci_vs_e6194_fb}
    # }

    k_our_fb = 'Our data (n=%d)' % (our_meta.type == 'FB').sum()
    k_our_ipsc = 'Our data (n=%d)' % (our_meta.type == 'iPSC').sum()
    k_e6194_fb = 'E-MTAB-6194 (n=1)'
    k_e6194_ipsc = 'E-MTAB-6194 (n=%d)' % len(ipsc_ref_names_6194_n1)
    k_hipsci_ipsc = 'HipSci (n=%d)' % hip_epic_meta.shape[0]

    to_plot = {
        k_our_ipsc: {
            k_our_fb: dmr_res_our_ipsc_vs_our_fb,
            k_e6194_fb: dmr_res_our_ipsc_vs_e6194_fb_n1
        },
        k_e6194_ipsc: {
            k_our_fb: dmr_res_e6194_ipsc_n1_vs_our_fb,
            k_e6194_fb: dmr_res_e6194_ipsc_n1_vs_e6194_fb_n1
        },
        k_hipsci_ipsc: {
            k_our_fb: dmr_res_hipsci_vs_our_fb,
            k_e6194_fb: dmr_res_hipsci_vs_e6194_fb_n1}
    }


    fig_hypo, axs_hypo = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(3.5, 6))
    fig_hyper, axs_hyper = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(3.5, 6))
    # colours = ['#91e097', '#e09191']

    marker_styles = ['o', 's', '>', '<', '^', 'v', 'X', 'P', 'd', 'H', '*', 'p']

    medianprops = dict(linestyle='-', linewidth=2., color='k')

    ymax = {}

    for i, (k1, d) in enumerate(to_plot.items()):
        axs_hyper[i, 0].set_ylabel(k1)
        axs_hypo[i, 0].set_ylabel(k1)
        for j, (k2, obj) in enumerate(d.items()):
            t = pd.DataFrame(get_dmr_number_direction(obj)).transpose()
            for axs, k3 in zip([axs_hyper, axs_hypo], ['Hypermethylated', 'Hypomethylated']):
                this_y = t.loc[:, k3].values
                this_x = np.random.normal(scale=0.05, size=this_y.size) + 1
                ymax[k3] = max(ymax.get(k3, 0), this_y.max())
                axs[-1, j].set_xlabel(k2)
                ax = axs[i, j]
                bplot = ax.boxplot(
                    this_y,
                    vert=True,
                    patch_artist=True,
                    medianprops=medianprops,
                    widths=0.7
                )
                # for xx, yy, mm in zip(this_x, this_y, marker_styles):
                ax.scatter(this_x, this_y, facecolor='none', edgecolor='k', linewidths=1., s=40, marker='o', zorder=5)
                ax.set_xticks([])
                plt.setp(bplot['boxes'], facecolor=colours[k3])
                # for p, c in zip(bplot['boxes'], colours):
                #     p.set_facecolor(c)
    axs_hyper[0, 0].set_ylim([0, ymax['Hypermethylated'] * 1.1])
    axs_hypo[0, 0].set_ylim([0, ymax['Hypomethylated'] * 1.1])
    fig_hyper.tight_layout()
    fig_hypo.tight_layout()
    fig_hyper.savefig(os.path.join(outdir, "iPSC_vs_FB_number_dmr_hyper.png"), dpi=200)
    fig_hypo.savefig(os.path.join(outdir, "iPSC_vs_FB_number_dmr_hypo.png"), dpi=200)

    # iPSC vs FB: core DMRs in each comparison type

    members = {'hypo': {}, 'hyper': {}}
    for i, (k1, d) in enumerate(to_plot.items()):
        for j, (k2, obj) in enumerate(d.items()):
            t = get_dmr_cid_direction(obj)
            hypo_cid = setops.reduce_intersection(*[x['Hypomethylated'] for x in t.values()])
            hyper_cid = setops.reduce_intersection(*[x['Hypermethylated'] for x in t.values()])
            members['hypo'][(k1, k2)] = hypo_cid
            members['hyper'][(k1, k2)] = hyper_cid

    ks = sorted(members['hypo'].keys())
    upset_dat = [members['hypo'][k] for k in ks]
    upset_lbl = ["%s - %s" % t for t in ks]
    ups = venn.upset_set_size_plot(
        upset_dat, upset_lbl,
    )
    ups['gs'].update(wspace=0.7)
    ups['axes']['set_size'].set_xlabel('Number hypomethylated DMRs in iPSC')
    ups['axes']['main'].set_ylabel('Number DMRs')
    ups['figure'].savefig(os.path.join(outdir, "ipsc_vs_fb_core_hypo_dmr_upset.png"), dpi=200)

    ks = sorted(members['hyper'].keys())
    upset_dat = [members['hyper'][k] for k in ks]
    upset_lbl = ["%s - %s" % t for t in ks]
    ups = venn.upset_set_size_plot(
        upset_dat, upset_lbl,
    )
    ups['gs'].update(wspace=0.7)
    ups['axes']['set_size'].set_xlabel('Number hypermethylated DMRs in iPSC')
    ups['axes']['main'].set_ylabel('Number DMRs')
    ups['figure'].savefig(os.path.join(outdir, "ipsc_vs_fb_core_hyper_dmr_upset.png"), dpi=200)

    # iPSC vs ESC: numbers and direction

    to_plot = {
        'Our iPSC': dmr_res_our_ipsc_vs_esc,
        'E-MTAB-6194 iPSC': dmr_res_e6914_ipsc_vs_esc,
        'HipSci iPSC': dmr_res_hipsci_vs_esc
    }
    ax1 = None
    gs = plt.GridSpec(6, 3)
    fig = plt.figure(figsize=(5, 6))
    axs = np.empty((3, 3), dtype=object)
    colours = ['#91e097', '#e09191']
    medianprops = dict(linestyle='-', linewidth=2., color='k')



    for i, (k1, obj) in enumerate(to_plot.items()):
        ax_i = (2 * i)
        t = pd.DataFrame(get_dmr_number_direction(obj)).transpose()
        u = get_dmr_cid_direction(obj)
        u_hypo = {}
        u_hyper = {}
        for j, r in enumerate(esc_ref_names):
            # ax_ix = 3 * i + j + 1
            if ax1 is None:
                # ax = fig.add_subplot(3, 3, ax_ix)
                ax = fig.add_subplot(gs[ax_i:(ax_i + 2), j])
                ax1 = ax
            else:
                # ax = fig.add_subplot(3, 3, ax_ix, sharex=ax1, sharey=ax1)
                ax = fig.add_subplot(gs[ax_i:(ax_i + 2), j], sharex=ax1, sharey=ax1)

            if j != 0:
                plt.setp(
                    ax.yaxis.get_ticklabels(), visible=False
                )

            axs[i, j] = ax
            u_hypo[r] = setops.reduce_intersection(*[x['Hypomethylated'] for k, x in u.items() if r in k])
            u_hyper[r] = setops.reduce_intersection(*[x['Hypermethylated'] for k, x in u.items() if r in k])

            this_t = t.loc[t.index.str.contains(r)]
            # ax = axs[i, j]
            bplot = ax.boxplot(
                this_t.values,
                vert=True,
                patch_artist=True,
                medianprops=medianprops,
                widths=0.7
            )
            ax.set_xticks([])
            for p, c in zip(bplot['boxes'], colours):
                p.set_facecolor(c)

        # add third column with Venn diagrams
        ax = fig.add_subplot(gs[ax_i + 1, 2], facecolor='none', frame_on=False, xticks=[], yticks=[])
        v = venn.venn_diagram(
            *u_hypo.values(),
            ax=ax,
            set_labels=u_hypo.keys(),
            set_colors=[colours[0]] * 2
        )
        ax = fig.add_subplot(gs[ax_i, 2], facecolor='none', frame_on=False, xticks=[], yticks=[])
        v = venn.venn_diagram(
            *u_hyper.values(),
            ax=ax,
            set_labels=u_hyper.keys(),
            set_colors=[colours[1]] * 2,
            alpha=0.7
        )

    for i, (k1, obj) in enumerate(to_plot.items()):
        axs[i, 0].set_ylabel(k1)

    for j in range(len(esc_ref_names)):
        ax = axs[-1, j]
        ax.set_xlabel(esc_ref_names[j])

    gs.update(bottom=0.05, top=0.96, right=0.97, wspace=0.05)
    fig.savefig(os.path.join(outdir, "iPSC_vs_ESC_number_dmr.png"), dpi=200)

    # final piece of information: is there any overlap between the resultant core sets?
    core_dmrs_hypo = {}
    core_dmrs_hyper = {}
    for k1, obj in to_plot.items():
        u = get_dmr_cid_direction(obj)
        u_hypo = {}
        u_hyper = {}

        core_dmrs_hypo[k1] = setops.reduce_intersection(*[
            setops.reduce_intersection(*[x['Hypomethylated'] for k, x in u.items() if r in k])
            for r in esc_ref_names
            ])
        core_dmrs_hyper[k1] = setops.reduce_intersection(*[
            setops.reduce_intersection(*[x['Hypermethylated'] for k, x in u.items() if r in k])
            for r in esc_ref_names
            ])

    # outcome
    vs, vc = setops.venn_from_arrays(*core_dmrs_hypo.values())
    print "Hypomethylated core DMRs (hypo in both ESC comparisons). "
    print "Of the %d DMRs in our data, %d are shared with both HipSci and E-MTAB-6194" % (
        len(core_dmrs_hypo['Our iPSC']),
        vc['111']
    )


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
