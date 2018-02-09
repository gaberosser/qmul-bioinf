import multiprocessing as mp
import os
import pandas as pd
import numpy as np
from utils import output, setops, excel
from methylation import dmr, process, loader as methylation_loader
from rnaseq import loader as rnaseq_loader, differential_expression, general
from integrator import rnaseq_methylationarray
from analysis import cross_comparison
from load_data import loader
from plotting import venn


def load_methylation(pids, ref_names=None, norm_method='swan', ref_name_filter=None):
    """
    Load and prepare the Illumina methylation data
    """
    # patient data
    obj = methylation_loader.load_by_patient(pids, norm_method=norm_method)
    anno = methylation_loader.load_illumina_methylationepic_annotation()

    # reference data
    if ref_names is not None:
        ref_obj = methylation_loader.load_reference(ref_names, norm_method=norm_method)
        if ref_name_filter is not None:
            ref_obj.filter_by_sample_name(ref_name_filter, exact=True)
        obj = loader.MultipleBatchLoader([obj, ref_obj])

    me_data = obj.data.dropna()
    me_data = process.m_from_beta(me_data)

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


def paired_dmr(me_data, me_meta, anno, pids, dmr_params):
    # Compute DMR for paired comparisons (only - no other comparisons needed at present)
    dmr_res_obj = dmr.DmrResults(anno=anno)
    dmr_res_obj.identify_clusters(**dmr_params)
    dmr_res = {}

    for pid in pids:
        the_idx1 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == 'GBM')
        the_idx2 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == 'iNSC')
        the_idx = the_idx1 | the_idx2
        the_groups = me_meta.loc[the_idx, 'type'].values
        the_samples = me_meta.index[the_idx].groupby(the_groups).values()
        the_obj = dmr_res_obj.copy()
        the_obj.test_clusters(me_data,
                              samples=the_samples,
                              n_jobs=dmr_params['n_jobs'],
                              min_median_change=dmr_params['delta_m_min'],
                              method=dmr_params['dmr_test_method'],
                              alpha=dmr_params['alpha'],
                              **dmr_params['test_kwargs']
                              )
        dmr_res[pid] = the_obj

    return dmr.DmrResultCollection(**dmr_res)


def dmr_results_hash(pids, dmr_params):
    hash_elements = tuple(sorted(pids)) + (
        dmr_params['d_max'],
        dmr_params['n_min'],
        dmr_params['delta_m_min'],
        dmr_params['dmr_test_method'],
        dmr_params['alpha'],
        dmr_params['norm_method']
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


def annotate_de_dmr_wide_form(data):
    data.insert(0, 'dmr_cluster_id', [t[1] for t in data.index])
    data.insert(0, 'gene', [t[0] for t in data.index])


def expand_dmr_results_table_by_gene(tbl):

    this_keep = tbl.loc[tbl.genes.apply(len) == 1]
    this_keep.insert(0, 'gene_symbol', tbl.genes.apply(lambda t: t[0]))
    this_keep = this_keep.drop('genes', axis=1)
    this_keep.index = zip(this_keep.index, this_keep.gene_symbol)

    to_expand = tbl.loc[tbl.genes.apply(len) > 1]
    expanded = []
    for cl_id, row in to_expand.iterrows():
        the_genes = row.genes
        row = row.drop('genes')
        for g in the_genes:
            new_row = row.copy()
            new_row.loc['gene_symbol'] = g
            expanded.append(new_row)

    expanded = pd.DataFrame(expanded)
    expanded.index = zip(expanded.index, expanded.gene_symbol)

    return pd.concat((this_keep, expanded), axis=0)


if __name__ == "__main__":
    outdir = output.unique_output_dir("hgic_de_dmr_two_strategies", reuse_empty=True)
    outdir_s1 = os.path.join(outdir, 's1')
    outdir_s2 = os.path.join(outdir, 's2')
    if not os.path.isdir(outdir_s1):
        os.makedirs(outdir_s1)
    if not os.path.isdir(outdir_s2):
        os.makedirs(outdir_s2)

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'GLM'
    }

    dmr_params = {
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'alpha': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': mp.cpu_count(),
    }
    norm_method_s1 = 'swan'

    pids = ['018', '019', '031', '017', '050', '054']

    subgroups = {
        'RTK I': ['018', '019', '031'],
        'RTK II': ['017', '050', '054'],
    }
    # indicator showing which groups the PIDs belong to
    subgroup_ind = dict([
        (k, pd.Index(pids).isin(v)) for k, v in subgroups.items()
    ])

    external_ref_names_de = ['GSE61794']
    external_ref_names_dm = ['gse38216']
    external_ref_samples_dm = ['H9 NPC 1', 'H9 NPC 2']

    external_refs_de = [
        ('GIBCO', 'NSC'),
        ('H9', 'NSC'),
    ]
    external_refs_dm = [
        ('GIBCO', 'NSC'),
        ('H9', 'NSC'),
    ]

    # plotting parameters
    cmap = 'RdYlGn_r'

    # file location parameters
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    ##################
    ### Strategy 1 ###
    ##################

    # load methylation
    # For S1, we just need the paired comparison. This is important - any other samples lead to a change in the final
    # probe list (any NA rows are dropped from ALL samples).

    me_obj, anno = load_methylation(pids, norm_method=norm_method_s1)
    me_data = me_obj.data
    me_meta = me_obj.meta

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    the_hash = dmr_results_hash(me_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        dmr_res_s1 = paired_dmr(me_data, me_meta, anno, pids, dmr_params)
        # Save DMR results to disk
        dmr_res_s1.to_pickle(fn, include_annotation=False)
        print "Saved DMR results to %s" % fn

    dmr_classes = ['tss', 'island']

    # extract combined significant results for tss and island
    dmr_res_all_cls = dmr_res_s1.results_by_class()
    dmr_res_sign_all_cls = dmr_res_s1.results_significant_by_class()

    dmr_res_full_s1 = {}
    dmr_res_sign_s1 = {}
    for k1 in dmr_res_sign_all_cls:
        # for k2 in dmr_res_sign_all_cls[k1]:
        dmr_res_sign_s1[k1] = {}
        dmr_res_full_s1[k1] = {}
        for dc in dmr_classes:
            dmr_res_sign_s1[k1].update(dmr_res_sign_all_cls[k1][dc])
            dmr_res_full_s1[k1].update(dmr_res_all_cls[k1][dc])
        # only need the keys for the significant DMR list
        dmr_res_sign_s1[k1] = dmr_res_sign_s1[k1].keys()

    rnaseq_obj = load_rnaseq(pids, external_ref_names_de)

    # Compute DE cross-comparison

    de_res_full = differential_expression.compute_cross_de(
        rnaseq_obj,
        pids,
        external_references=external_refs_de,
        return_full=True,
        **de_params
    )

    # extract the subset of significant results from the full dataframes
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

    # expanded core (genes DE in multiple subgroups)
    ec_sets = expanded_core_sets(venn_set, subgroup_ind)
    data_ec = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, include_sets=ec_sets)
    data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_de.xlsx'))

    # subgroup-specific
    ss_sets = []
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        ss_sets.append(k)
    data_ss = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, include_sets=ss_sets)
    data_ss.to_excel(os.path.join(outdir_s1, 'subgroup_specific_de.xlsx'))

    # patient unique
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))
    data_pu = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, include_sets=pu_sets)
    data_pu.to_excel(os.path.join(outdir_s1, 'patient_unique_de.xlsx'))

    # b) DMR only

    # dmr_res_sign_s1 = dict([
    #     (pid, dmr_res_sign[(pid, pid)]) for pid in pids
    # ])
    # dmr_res_full_s1 = dict([
    #     (pid, dmr_res_full[(pid, pid)]) for pid in pids
    # ])

    dmr_by_member = [dmr_res_sign_s1[pid] for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*dmr_by_member)

    # add null set manually from full DMR results
    dmr_id_all = setops.reduce_union(*venn_set.values())
    k_null = ''.join(['0'] * len(pids))
    venn_set[k_null] = list(set(dmr_res_full_s1[pids[0]].keys()).difference(dmr_id_all))
    venn_ct[k_null] = len(venn_set[k_null])

    upset1 = upset_plot_dmr(
        dmr_by_member, venn_set, subgroup_ind, pids
    )

    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_dmr.png"), dpi=200)
    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_dmr.tiff"), dpi=200)

    # generate wide-form lists and save to Excel file
    data_for_dmr_table = {}
    data_for_dmr_table_full = {}
    for pid in pids:
        this_sign = dmr_res_s1[pid].to_table(include='significant', skip_geneless=True)
        data_for_dmr_table[pid] = expand_dmr_results_table_by_gene(this_sign)

        this_full = dmr_res_s1[pid].to_table(include='all',skip_geneless=True)
        data_for_dmr_table_full[pid] = expand_dmr_results_table_by_gene(this_full)

    # data_for_dmr_table_full = dict([
    #     (pid, dmr_res[pid][pid].to_table(include='all',skip_geneless=False)) for pid in pids
    # ])

    # recalculate venn set
    dmr_by_member = [data_for_dmr_table[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*dmr_by_member)

    # add null set manually from full DMR results
    dmr_id_all = setops.reduce_union(*venn_set.values())
    k_null = ''.join(['0'] * len(pids))
    venn_set[k_null] = list(data_for_dmr_table_full[pids[0]].index.difference(dmr_id_all))
    venn_ct[k_null] = len(venn_set[k_null])

    data = setops.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        cols_to_include=('median_delta', 'padj'),
        consistency_check_col='median_delta',
        consistency_check_method='sign'
    )
    data.insert(0, 'gene_symbol', [t[1] for t in data.index])
    data.insert(0, 'cluster_id', [t[0] for t in data.index])

    # quick (bespoke) sanity check
    for pid in pids:
        xx = data.loc[:, data.columns.str.contains(pid)]
        x1 = xx.loc[xx.loc[:, pid] == 'Y']
        x0 = xx.loc[xx.loc[:, pid] == 'N']
        if x1.loc[:, "%s_padj" % pid].isnull().sum() != 0:
            print "PID %s has failing entries in x1" % pid
        if (x0.loc[:, "%s_padj" % pid] < dmr_params['alpha']).sum() != 0:
            print "PID %s has failing entries in x0" % pid

    data.to_excel(os.path.join(outdir_s1, 'full_dmr.xlsx'))

    # expanded core
    ec_sets = expanded_core_sets(venn_set, subgroup_ind)
    data_ec = setops.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        include_sets=ec_sets,
        cols_to_include=('median_delta', 'padj'),
        consistency_check_col='median_delta',
        consistency_check_method='sign'
    )
    data_ec.insert(0, 'gene_symbol', [t[1] for t in data_ec.index])
    data_ec.insert(0, 'cluster_id', [t[0] for t in data_ec.index])
    data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_dmr.xlsx'))

    # subgroup-specific
    ss_sets = []
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        ss_sets.append(k)
    data_ss = setops.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        include_sets=ss_sets,
        cols_to_include=('median_delta', 'padj'),
        consistency_check_col='median_delta',
        consistency_check_method='sign'
    )
    data_ss.insert(0, 'gene_symbol', [t[1] for t in data_ss.index])
    data_ss.insert(0, 'cluster_id', [t[0] for t in data_ss.index])
    data_ss.to_excel(os.path.join(outdir_s1, 'subgroup_specific_dmr.xlsx'))

    # patient unique
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))
    data_pu = setops.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        include_sets=pu_sets,
        cols_to_include=('median_delta', 'padj'),
        consistency_check_col='median_delta',
        consistency_check_method='sign'
    )
    data_pu.insert(0, 'gene_symbol', [t[1] for t in data_pu.index])
    data_pu.insert(0, 'cluster_id', [t[0] for t in data_pu.index])
    data_pu.to_excel(os.path.join(outdir_s1, 'patient_unique_dmr.xlsx'))

    ####
    #### c) Layered DE and DMR
    ####

    # get the joint table
    # joint_de_dmr_s1 = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_s1)
    joint_de_dmr_s1 = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_s1)
    # filter - only DMRs with TSS or island class
    for pid in pids:
        idx = joint_de_dmr_s1[pid].dmr_class_island | joint_de_dmr_s1[pid].dmr_class_tss
        joint_de_dmr_s1[pid] = joint_de_dmr_s1[pid].loc[idx]
        joint_de_dmr_s1[pid].index = de_dmr_hash(joint_de_dmr_s1[pid])
        print "Patient %s. %d DE & DMR rows, %d are gene only, keeping %d." % (pid, idx.size, (~idx).sum(), idx.sum())

    de_dmr_by_member = [joint_de_dmr_s1[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*de_dmr_by_member)

    # we don't add the null set here - it takes too long

    # generate wide-form lists and save to Excel file
    # first, export everything (concordant and discordant)
    data = dmr.venn_set_to_wide_dataframe(
        joint_de_dmr_s1,
        venn_set,
        pids,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        direction_col='de_logfc'
    )
    # add gene symbol and cluster ID back in
    annotate_de_dmr_wide_form(data)
    data.to_excel(os.path.join(outdir_s1, 'full_de_dmr.xlsx'), index=False)

    # now filter by concordance and recompute venn sets
    joint_de_dmr_concordant_s1 = {}
    for pid in pids:
        a = joint_de_dmr_s1[pid].de_direction
        b = joint_de_dmr_s1[pid].dmr_direction
        idx = (a != b)
        joint_de_dmr_concordant_s1[pid] = joint_de_dmr_s1[pid].loc[idx]
        print "Patient %s has %d combined DE / DMR results, of which %d are concordant" % (
            pid, idx.size, idx.sum()
        )
    de_dmr_by_member_concordant = [joint_de_dmr_concordant_s1[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*de_dmr_by_member_concordant)

    # export only concordant results
    data = dmr.venn_set_to_wide_dataframe(
        joint_de_dmr_concordant_s1,
        venn_set,
        pids,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        direction_col='de_logfc'
    )
    # add gene symbol and cluster ID back in
    annotate_de_dmr_wide_form(data)
    data.to_excel(os.path.join(outdir_s1, 'full_de_dmr_concordant.xlsx'), index=False)

    # expanded core
    ec_sets = expanded_core_sets(venn_set, subgroup_ind)
    data_ec = dmr.venn_set_to_wide_dataframe(
        joint_de_dmr_concordant_s1,
        venn_set,
        pids,
        include_sets=ec_sets,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        direction_col='de_logfc'
    )
    annotate_de_dmr_wide_form(data_ec)
    data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_de_dmr_concordant.xlsx'), index=False)

    # subgroup-specific
    ss_sets = []
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        ss_sets.append(k)
    data_ss = dmr.venn_set_to_wide_dataframe(
        joint_de_dmr_concordant_s1,
        venn_set,
        pids,
        include_sets=ss_sets,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        direction_col='de_logfc'
    )
    annotate_de_dmr_wide_form(data_ss)
    data_ss.to_excel(os.path.join(outdir_s1, 'subgroup_specific_de_dmr_concordant.xlsx'), index=False)

    # patient unique
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))
    data_pu = dmr.venn_set_to_wide_dataframe(
        joint_de_dmr_concordant_s1,
        venn_set,
        pids,
        include_sets=pu_sets,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        direction_col='de_logfc'
    )
    annotate_de_dmr_wide_form(data_pu)
    data_pu.to_excel(os.path.join(outdir_s1, 'patient_unique_de_dmr_concordant.xlsx'), index=False)

    upset1 = upset_plot_de(
        de_dmr_by_member_concordant, venn_set, subgroup_ind, pids
    )
    upset1['axes']['main'].set_ylabel('Number of DM-concordant DE genes in set')
    upset1['axes']['set_size'].set_xlabel('Number of DM-concordant DE genes\nin single comparison')
    upset1['gs'].update(bottom=0.11)

    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_de_dmr.png"), dpi=200)
    upset1['figure'].savefig(os.path.join(outdir_s1, "upset_de_dmr.tiff"), dpi=200)

    ##################
    ### Strategy 2 ###
    ##################

    norm_method_s2 = 'pbc'

    # load methylation with external references
    me_obj_with_ref, anno = load_methylation(
        pids,
        ref_names=external_ref_names_dm,
        ref_name_filter=external_ref_samples_dm,
        norm_method=norm_method_s2
    )
    me_data_with_ref = me_obj_with_ref.data

    # Compute DMR cross-comparison

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s2

    the_hash = dmr_results_hash(me_obj_with_ref.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_cross_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        dmr_res_s2 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        dmr_res_s2 = dmr.compute_cross_dmr(
            me_data_with_ref,
            me_obj_with_ref.meta,
            anno,
            pids,
            dmr_params,
            external_references=external_refs_dm
        )
        # Save DMR results to disk
        dmr_res_s2.to_pickle(fn, include_annotation=False)
        print "Saved DMR results to %s" % fn

    dmr_classes = ['tss', 'island']

    # extract combined significant results for tss and island

    dmr_res_all_cls = dmr_res_s2.results_by_class()
    dmr_res_sign_all_cls = dmr_res_s2.results_significant_by_class()

    dmr_res_full_s2 = {}
    dmr_res_sign_s2 = {}
    for k1 in dmr_res_sign_all_cls:
        for k2 in dmr_res_sign_all_cls[k1]:
            dmr_res_sign_s2[(k1, k2)] = {}
            dmr_res_full_s2[(k1, k2)] = {}
            for dc in dmr_classes:
                dmr_res_sign_s2[(k1, k2)].update(dmr_res_sign_all_cls[k1][k2][dc])
                dmr_res_full_s2[(k1, k2)].update(dmr_res_all_cls[k1][k2][dc])
            dmr_res_sign_s2[(k1, k2)] = dmr_res_sign_s2[(k1, k2)].keys()

    ## a) Pair-only DE

    de_res_s2_idx = dict([
        (k, v.index) for k, v in de_res.items()
    ])

    cc_dict_de = cross_comparison.compute_cross_comparison_correction(
        de_res_s2_idx,
        pids,
        external_refs=[t[0] for t in external_refs_de],
        set_type='pair_only'
    )
    po_specific_to_all_refs_de = sorted(cc_dict_de['specific_to_all_refs'])
    pair_only_de = cc_dict_de['venn_set']

    for_export = rnaseq_obj.data.loc[po_specific_to_all_refs_de]
    general.add_gene_symbols_to_ensembl_data(for_export)
    for_export.to_excel(os.path.join(outdir_s2, 'consistently_in_pair_only_across_all_refs_de.xlsx'))

    # take the intersection of the PO lists, then export to a file
    po_de_export = {}
    po_de_export_full = {}
    for pid in pids:
        po_de_export_full[pid] = de_res_full[(pid, pid)].copy()

        this_row = pair_only_de.loc[pid, [t[0] for t in external_refs_de]]
        this_genes_pre = setops.reduce_intersection(*this_row)

        # can correct as follows
        # this_genes = sorted(this_genes_pre.difference(po_specific_to_all_refs))
        # print "PID %s. Subtracted %d correction genes from the %d PO intersection genes to leave %d PO genes" % (
        #     pid, len(po_specific_to_all_refs), len(this_genes_pre), len(this_genes)
        # )

        # we won't, to avoid cross-comparison issues
        this_genes = this_genes_pre

        po_col = pd.Series('N', index=po_de_export_full[pid].index)
        po_col.loc[this_genes] = 'Y'
        po_de_export_full[pid].insert(po_de_export_full[pid].shape[1], 'pair_only', po_col)

        po_de_export[pid] = de_res[(pid, pid)].loc[this_genes]

    excel.pandas_to_excel(po_de_export, os.path.join(outdir_s2, 'pair_only_de.xlsx'))
    excel.pandas_to_excel(po_de_export_full, os.path.join(outdir_s2, 'pair_only_de_full.xlsx'))

    # export with a different layout, analogous to strategy 1
    venn_set, venn_ct = setops.venn_from_arrays(*[po_de_export[pid].index for pid in pids])
    po_combination_export = differential_expression.venn_set_to_dataframe(po_de_export, venn_set, pids)

    po_combination_export = setops.venn_set_to_wide_dataframe(
        po_de_export,
        venn_set,
        pids,
        cols_to_include=('logFC', 'FDR'),
        consistency_check_col='logFC',
        consistency_check_method="sign"
    )

    po_combination_export.to_excel(os.path.join(outdir_s2, 'pair_only_de_wideform.xlsx'))

    ## b) Pair-only DMR

    # Compute cross-comparison correction
    cc_dict_dmr = cross_comparison.compute_cross_comparison_correction(
        dmr_res_sign_s2,
        pids,
        external_refs=[t[0] for t in external_refs_dm],
        set_type='pair_only'
    )

    po_specific_to_all_refs_dmr = sorted(cc_dict_dmr['specific_to_all_refs'])
    pair_only_dmr = cc_dict_dmr['venn_set']

    #for DMR, just dump some information to a file

    for_export = []
    for i in po_specific_to_all_refs_dmr:
        t = dmr_res_s2.clusters[i].summary_dict
        for_export.append('\n'.join(["%s: %s" % (k, str(v)) for k, v in t.items()]))
    with open(os.path.join(outdir_s2, 'consistently_in_pair_only_across_all_refs_dmr.txt'), 'wb') as f:
        f.write("\n\n***\n\n".join(for_export))

    # take the intersection of the PO lists, then export to a file
    po_dmr_export = {}
    po_dmr_export_full = {}
    for pid in pids:
        po_dmr_export_full[pid] = dmr_res_full_s2[(pid, pid)].copy()

        this_row = pair_only_dmr.loc[pid, [t[0] for t in external_refs_dm]]
        this_genes_pre = setops.reduce_intersection(*this_row)

        # no correction, to avoid cross-comparison issues
        this_genes = this_genes_pre

        po_col = pd.Series('N', index=po_dmr_export_full[pid].index)
        po_col.loc[this_genes] = 'Y'
        po_dmr_export_full[pid].insert(po_dmr_export_full[pid].shape[1], 'pair_only', po_col)

        ## TODO: the results are currently in a dictionary - can we export them to a table (or did we already??)
        # po_dmr_export[pid] = dmr_res_s2[(pid, pid)].loc[this_genes]

    excel.pandas_to_excel(po_de_export, os.path.join(outdir_s2, 'pair_only_dmr.xlsx'))
    excel.pandas_to_excel(po_de_export_full, os.path.join(outdir_s2, 'pair_only_dmr_full.xlsx'))

    # export with a different layout, analogous to strategy 1
    venn_set, venn_ct = setops.venn_from_arrays(*[po_dmr_export[pid].index for pid in pids])
    po_combination_export = setops.venn_set_to_wide_dataframe(
        po_dmr_export,
        venn_set,
        pids,
        cols_to_include=('median_delta', 'padj'),
        consistency_check_col='median_delta',
        consistency_check_method="sign"
    )

    po_combination_export.to_excel(os.path.join(outdir_s2, 'pair_only_dmr_wideform.xlsx'))

    ## c) DE / DMR combined
    # TODO