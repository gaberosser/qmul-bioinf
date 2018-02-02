import multiprocessing as mp
import os
import pandas as pd
import numpy as np
from utils import output, setops
from methylation import dmr, process, loader as methylation_loader
from rnaseq import loader as rnaseq_loader, differential_expression
from integrator import rnaseq_methylationarray
from analysis import cross_comparison
from load_data import loader
from plotting import venn


def load_methylation(pids, norm_method='swan'):
    """
    Load and prepare the Illumina methylation data
    """
    obj = methylation_loader.load_by_patient(pids, norm_method=norm_method)
    me_data = obj.data
    me_data.dropna(inplace=True)
    me_data = process.m_from_beta(me_data)
    anno = methylation_loader.load_illumina_methylationepic_annotation()

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(me_data.index)
    anno = anno.loc[common_probes]
    dmr.add_merged_probe_classes(anno)
    me_data = me_data.loc[common_probes]
    obj.data = me_data

    return obj, anno


def load_rnaseq(pids, ref_names, ref_name_filter='NSC', discard_filter='IPSC'):
    # Load RNA-Seq from STAR
    obj = rnaseq_loader.load_by_patient(pids)

    # load additional references
    ref_objs = []
    for rn in ref_names:
        ref_obj = rnaseq_loader.load_references(rn)
        if ref_name_filter is not None:
            # only keep relevant references
            ref_obj.meta = ref_obj.meta.loc[ref_obj.meta.index.str.contains(ref_name_filter)]
            ref_obj.data = ref_obj.data.loc[:, ref_obj.meta.index]
        ref_objs.append(ref_obj)
    obj = loader.MultipleBatchLoader([obj] + ref_objs)

    if discard_filter is not None:
        if not hasattr(discard_filter, '__iter__'):
            discard_filter = [discard_filter]
        for d in discard_filter:
            obj.meta = obj.meta.loc[~obj.meta.index.str.contains(d)]
            obj.data = obj.data.loc[:, obj.meta.index]

    return obj


def dmr_results_hash(pids, dmr_params):
    hash_elements = tuple(sorted(pids)) + (
        dmr_params['d_max'],
        dmr_params['n_min'],
        dmr_params['delta_m_min'],
        dmr_params['dmr_test_method'],
    ) + tuple(dmr_params['test_kwargs'].items())
    return hash(hash_elements)


def expanded_core_sets(venn_set, subgroup_ind):
    """
    Compute the sets that belong to the 'expanded core', which comprises any combination of patients that bridges
    multiple subgroups.
    :param venn_set: As returned by setops.venn_from_arrays
    :param subgroup_ind: Dict, keys are the names of the subgroups, values are boolean arrays with the membership for
    each patient. It's important that the ordering here is the same as that used downstream.
    :return:
    """
    ecs = []
    for k in venn_set:
        this_k = np.array([t for t in k]).astype(bool)
        nmatch = 0
        for grp, grp_idx in subgroup_ind.items():
            if this_k[grp_idx].any():
                nmatch += 1
        if nmatch > 1:
            # add to the expanded core set
            ecs.append(k)
    return ecs


def upset_plot_de(de_data, venn_set, subgroup_ind, set_labels):
    # UpsetR attribute plots

    # set colours for UpsetR plot
    sets_full = {}
    sets_partial = {}
    sets_unique = []

    for k in venn_set:
        this_k = np.array([t for t in k]).astype(bool)
        if this_k.sum() == 1:
            sets_unique.append(k)
        elif this_k.sum() == 2:
            for grp, grp_idx in subgroup_ind.items():
                if this_k[grp_idx].sum() == this_k.sum():
                    sets_partial.setdefault(grp, []).append(k)
        elif this_k.sum() == 3:
            for grp, grp_idx in subgroup_ind.items():
                if this_k[grp_idx].sum() == this_k.sum():
                    sets_full.setdefault(grp, []).append(k)

    set_colours = [
        ('RTK I full', {'sets': sets_full['RTK I'], 'colour': '#0d680f'}),
        ('RTK I partial', {'sets': sets_partial['RTK I'], 'colour': '#6ecc70'}),
        ('RTK II full', {'sets': sets_full['RTK II'], 'colour': '#820505'}),
        ('RTK II partial', {'sets': sets_partial['RTK II'], 'colour': '#d67373'}),
        ('Expanded core', {'sets': expanded_core_sets(venn_set, subgroup_ind), 'colour': '#4C72B0'}),
        ('Unique', {'sets': sets_unique, 'colour': '#f4e842'})
    ]

    # 1. Descending order
    return venn.upset_set_size_plot(
        de_data,
        set_labels,
        set_colours=set_colours,
        min_size=10,
        n_plot=30,
        default_colour='gray'
    )


def upset_plot_dmr(dmr_data, venn_set, subgroup_ind, set_labels):
    res = upset_plot_de(dmr_data, venn_set, subgroup_ind, set_labels)
    res['axes']['main'].set_ylabel('Number of DMRs in set')
    res['axes']['set_size'].set_xlabel('Number of DMRs in single comparison')
    return res


def de_dmr_hash(dat, cols=('gene', 'dmr_cid')):
    """
    Generate a list of hash values from the supplied DataFrame
    :param dat:
    :param cols:
    :return:
    """
    return list(dat.loc[:, cols].itertuples(index=False, name=None))


if __name__ == "__main__":
    outdir = output.unique_output_dir("hgic_de_dmr_two_strategies", reuse_empty=True)
    outdir_s1 = output.unique_output_dir("s1", root_output_dir=outdir, reuse_empty=True)
    outdir_s2 = output.unique_output_dir("s2", root_output_dir=outdir, reuse_empty=True)

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'GLM'
    }

    dmr_params = {
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'fdr': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': mp.cpu_count(),
    }

    pids = ['019', '030', '031', '017', '050', '054']

    subgroups = {
        'RTK I': ['019', '030', '031'],
        'RTK II': ['017', '050', '054'],
    }
    # indicator showing which groups the PIDs belong to
    subgroup_ind = dict([
        (k, pd.Index(pids).isin(v)) for k, v in subgroups.items()
    ])

    external_ref_names_de = ['GSE61794', 'GSE38993']

    external_refs_de = [
        ('GIBCO', 'NSC'),
        ('H9', 'NSC'),
        ('H1', 'NSC'),
    ]

    external_refs_dm = [
        'GIBCO'
    ]

    # plotting parameters
    cmap = 'RdYlGn_r'

    # file location parameters
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    # load methylation
    me_obj, anno = load_methylation(pids)
    me_data = me_obj.data

    # Compute DMR cross-comparison

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    the_hash = dmr_results_hash(pids, dmr_params)
    filename = 'dmr_results.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    loaded = False
    if os.path.isfile(fn):
        dmr_res = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        dmr_res = dmr.compute_cross_dmr(me_data, me_obj.meta, anno, pids, dmr_params)
        # Save DMR results to disk
        dmr_res.to_pickle(fn, include_annotation=False)
        print "Saved DMR results to %s" % fn

    dmr_classes = ['tss', 'island']

    # extract combined significant results for tss and island
    dmr_res_all_cls = dmr_res.results_by_class()
    dmr_res_sign_all_cls = dmr_res.results_significant_by_class()

    dmr_res_full = {}
    dmr_res_sign = {}
    for k1 in dmr_res_sign_all_cls:
        for k2 in dmr_res_sign_all_cls[k1]:
            dmr_res_sign[(k1, k2)] = {}
            dmr_res_full[(k1, k2)] = {}
            for dc in dmr_classes:
                dmr_res_sign[(k1, k2)].update(dmr_res_sign_all_cls[k1][k2][dc])
                dmr_res_full[(k1, k2)].update(dmr_res_all_cls[k1][k2][dc])
            dmr_res_sign[(k1, k2)] = dmr_res_sign[(k1, k2)].keys()

    rnaseq_obj = load_rnaseq(pids, external_ref_names_de)

    # Compute DE cross-comparison

    de_res_full = differential_expression.compute_cross_de(
        rnaseq_obj,
        pids,
        external_references=external_refs_de,
        return_full=True,
        **de_params
    )

    # no need to do this
    # de_res = differential_expression.compute_cross_de(
    #     rnaseq_obj,
    #     pids,
    #     external_references=external_refs_de,
    #     **de_params
    # )

    # have checked and it's identical to this faster option
    de_res = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full.items()])

    # STRATEGY 1: No references, just compare GBM-iNSC for each patient

    # a) DE only, so we can check the numbers
    de_res_full_s1 = dict([
        (pid, de_res_full[(pid, pid)]) for pid in pids
    ])
    de_res_s1 = dict([
        (pid, de_res[(pid, pid)]) for pid in pids
    ])
    de_by_member = [de_res_s1[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*de_by_member)

    # add null set manually from full DE results
    de_genes_all = setops.reduce_union(*venn_set.values())
    k_null = ''.join(['0'] * len(pids))
    venn_set[k_null] = list(de_res_full[(pids[0], pids[0])].index.difference(de_genes_all))
    venn_ct[k_null] = len(venn_set[k_null])

    upset1 = upset_plot_de(
        de_by_member, venn_set, subgroup_ind, pids
    )

    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_de.png"), dpi=200)
    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_de.tiff"), dpi=200)

    # generate wide-form lists and save to Excel file
    data = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, full_data=de_res_full_s1)
    data.to_excel(os.path.join(outdir_s1, 'full_de.xlsx'))

    # expanded core
    ec_sets = expanded_core_sets(venn_set, subgroup_ind)
    data_ec = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, include_sets=ec_sets)
    data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_de.xlsx'))

    # subgroup-specific
    ss_sets = []
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        ss_sets.append(k)
    data_ss = differential_expression.venn_set_to_dataframe(de_res, venn_set, pids, include_sets=ss_sets)
    data_ss.to_excel(os.path.join(outdir, 'subgroup_specific_de.xlsx'))

    # patient unique
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))
    data_pu = differential_expression.venn_set_to_dataframe(de_res, venn_set, pids, include_sets=pu_sets)
    data_pu.to_excel(os.path.join(outdir, 'patient_unique_de.xlsx'))

    # b) DMR only
    dmr_res_sign_s1 = dict([
        (pid, dmr_res_sign[(pid, pid)]) for pid in pids
    ])
    dmr_res_full_s1 = dict([
        (pid, dmr_res_full[(pid, pid)]) for pid in pids
    ])
    dmr_by_member = [dmr_res_sign_s1[pid] for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*dmr_by_member)

    # add null set manually from full DMR results
    dmr_id_all = setops.reduce_union(*venn_set.values())
    k_null = ''.join(['0'] * len(pids))
    venn_set[k_null] = list(set(dmr_res_full_s1[pids[0]].keys()).difference(dmr_id_all))
    venn_ct[k_null] = len(venn_set[k_null])

    # generate an expanded core gene set, defined as DMRs that are DM in both RTK I and RTK II (any number of patients)
    upset1 = upset_plot_dmr(
        dmr_by_member, venn_set, subgroup_ind, pids
    )

    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_dmr.png"), dpi=200)
    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_dmr.tiff"), dpi=200)

    # generate wide-form lists and save to Excel file
    data_for_dmr_table = dict([
        (pid, dmr_res[pid][pid].to_table(include='significant', skip_geneless=False)) for pid in pids
    ])
    data_for_dmr_table_full = dict([
        (pid, dmr_res[pid][pid].to_table(include='all',skip_geneless=False)) for pid in pids
    ])

    data = dmr.venn_set_to_wide_dataframe(data_for_dmr_table, venn_set, pids, full_data=data_for_dmr_table_full)

    # quick sanity check
    for pid in pids:
        xx = data.loc[:, data.columns.str.contains(pid)]
        x1 = xx.loc[xx.loc[:, pid] == 'Y']
        x0 = xx.loc[xx.loc[:, pid] == 'N']
        if x1.loc[:, "%s_padj" % pid].isnull().sum() != 0:
            print "PID %s has failing entries in x1" % pid
        if (x0.loc[:, "%s_padj" % pid] < dmr_params['fdr']).sum() != 0:
            print "PID %s has failing entries in x0" % pid

    data.to_excel(os.path.join(outdir_s1, 'full_dmr.xlsx'))

    # expanded core
    ec_sets = expanded_core_sets(venn_set, subgroup_ind)
    data_ec = dmr.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        include_sets=ec_sets
    )
    data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_dmr.xlsx'))

    # subgroup-specific
    ss_sets = []
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        ss_sets.append(k)
    data_ss = dmr.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        include_sets=ss_sets
    )
    data_ss.to_excel(os.path.join(outdir_s1, 'subgroup_specific_dmr.xlsx'))

    # patient unique
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))
    data_pu = dmr.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        include_sets=pu_sets
    )
    data_pu.to_excel(os.path.join(outdir_s1, 'patient_unique_dmr.xlsx'))

    # c) Layered DE and DMR
    dmr_res_s1 = dict([(pid, dmr_res[pid][pid]) for pid in pids])

    # get the joint table
    joint_de_dmr_s1 = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_s1)

    # filter - only DMRs with TSS or island class
    joint_de_dmr_tss_island_s1 = {}
    for pid in pids:
        idx = joint_de_dmr_s1[pid].dmr_class_island | joint_de_dmr_s1[pid].dmr_class_tss
        joint_de_dmr_tss_island_s1[pid] = joint_de_dmr_s1[pid].loc[idx]
        print "Patient %s. %d DE & DMR rows, %d are gene only" % (pid, idx.size, (~idx).sum())


    de_dmr_by_member = [de_dmr_hash(joint_de_dmr_tss_island_s1[pid]) for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*de_dmr_by_member)

    ## TODO

    # is this necessary? It takes AGES
    # joint_de_dmr_s1_full = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_full_s1, dmr_include='all')

    # add null set manually from full DMR results
    # de_dmr_id_all = setops.reduce_union(*venn_set.values())
    # k_null = ''.join(['0'] * len(pids))
    # venn_set[k_null] = list(set(dmr_res_full_s1[pids[0]].keys()).difference(dmr_id_all))
    # venn_ct[k_null] = len(venn_set[k_null])

    upset1 = upset_plot_dmr(
        de_dmr_by_member, venn_set, subgroup_ind, pids
    )

    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_de_dmr.png"), dpi=200)
    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_de_dmr.tiff"), dpi=200)

    # generate wide-form lists and save to Excel file
    # change the DE DMR table index for something uniquely identifying
    for pid in pids:
        joint_de_dmr_tss_island_s1[pid].index = de_dmr_hash(joint_de_dmr_tss_island_s1[pid])
    data = dmr.venn_set_to_wide_dataframe(
        joint_de_dmr_tss_island_s1,
        venn_set,
        pids,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        direction_col='de_logfc'
    )
    ## FIXME: 50% are not consistent?? Check this.
    data.to_excel(os.path.join(outdir_s1, 'full_de.xlsx'))

    ## TODO: finish these lines
    # expanded core
    # ec_sets = expanded_core_sets(venn_set, subgroup_ind)
    # data_ec = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, include_sets=ec_sets)
    # data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_de.xlsx'))

    # subgroup-specific
    # ss_sets = []
    # for grp in subgroup_ind:
    #     k = ''.join(subgroup_ind[grp].astype(int).astype(str))
    #     ss_sets.append(k)
    # data_ss = differential_expression.venn_set_to_dataframe(de_res, venn_set, pids, include_sets=ss_sets)
    # data_ss.to_excel(os.path.join(outdir, 'subgroup_specific_de.xlsx'))

    # patient unique
    # pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))
    # data_pu = differential_expression.venn_set_to_dataframe(de_res, venn_set, pids, include_sets=pu_sets)
    # data_pu.to_excel(os.path.join(outdir, 'patient_unique_de.xlsx'))



    # Compute DMR cross-comparison correction
    dm_specific_to_all_refs = cross_comparison.compute_cross_comparison_correction(
        dmr_res_sign,
        pids,
        external_refs=external_refs_dm,
        set_type='pair_only'
    )