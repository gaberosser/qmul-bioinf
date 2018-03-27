import multiprocessing as mp
import os
import pandas as pd
import numpy as np
from utils import output, setops, excel
from methylation import dmr, process, loader as methylation_loader
from rnaseq import loader as rnaseq_loader, differential_expression, general, filter
from integrator import rnaseq_methylationarray
from analysis import cross_comparison
from load_data import loader
from plotting import venn
from matplotlib import pyplot as plt
import seaborn as sns


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


def load_rnaseq(pids, ref_names, ref_name_filter='NSC', discard_filter='IPSC', strandedness=None):
    """
    :param strandedness: Iterable of same length as ref_names giving the strandedness of each ref
    """
    if strandedness is None:
        strandedness = ['u'] * len(ref_names)
    else:
        if len(strandedness) != len(ref_names):
            raise ValueError("Supplied strandedness must be a list of the same length as the ref_names.")

    # Load RNA-Seq from STAR
    obj = rnaseq_loader.load_by_patient(pids)

    # load additional references
    ref_objs = []
    for rn, strnd in zip(ref_names, strandedness):
        ref_obj = rnaseq_loader.load_references(rn, strandedness=strnd)
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

    obj.batch_id = obj.batch_id.loc[obj.meta.index]

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


def expand_dmr_results_table_by_gene(tbl, drop_geneless=False):

    this_keep = tbl.loc[tbl.genes.apply(len) == 1]
    this_keep.insert(0, 'gene_symbol', this_keep.genes.apply(lambda t: t[0]))
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

    if drop_geneless:
        return pd.concat((this_keep, expanded), axis=0)
    else:
        this_geneless = tbl.loc[tbl.genes.apply(len) == 0]
        this_geneless.insert(0, 'gene_symbol', None)
        this_geneless = this_geneless.drop('genes', axis=1)
        this_geneless.index = zip(this_geneless.index, this_geneless.gene_symbol)
        return pd.concat((this_keep, expanded, this_geneless), axis=0)


def patient_and_subgroup_specific_set(pids, subgroups):
    """
    For each PID, get the list of sets that contains that patient and is specific to the right subgroup
    :param pids:
    :param subgroups:
    :return: Dict, one entry per PID. Values are the lists of corresponding sets.
    """
    subgroup_ind = dict([
        (k, pd.Index(pids).isin(v)) for k, v in subgroups.items()
    ])
    candidates = list(setops.binary_combinations(len(pids), include_zero=False))
    # for each PID, what subgroup are they in?
    subgroup_membership = {}
    res = {}
    for i, pid in enumerate(pids):
        res[pid] = []
        for sg, arr in subgroup_ind.items():
            if arr[i]:
                if pid in subgroup_membership:
                    raise ValueError("Patient %s matches two or more subgroups: %s %s." % (
                        str(pid),
                        subgroup_membership[pid],
                        sg
                    ))
                subgroup_membership[pid] = sg
        if pid not in subgroup_membership:
            raise ValueError("Patient %s has no subgroup." % str(pid))

        not_sg_idx = np.where(~subgroup_ind[subgroup_membership[pid]])[0]
        for c in candidates:
            c = np.array([t for t in c])
            if any(c[not_sg_idx] == '1'):
                # Not subgroup specific
                continue
            if c[i] == '1':
                res[pid].append(''.join(c))
    return res


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


if __name__ == "__main__":
    outdir = output.unique_output_dir("hgic_de_dmr_two_strategies", reuse_empty=True)
    outdir_s1 = os.path.join(outdir, 's1')
    outdir_s2 = os.path.join(outdir, 's2')
    if not os.path.isdir(outdir_s1):
        os.makedirs(outdir_s1)
    if not os.path.isdir(outdir_s2):
        os.makedirs(outdir_s2)

    min_cpm = 1.

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'QLGLM'
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

    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']

    subgroups = {
        'RTK I': ['018', '019', '030', '031'],
        'RTK II': ['017', '050', '054', '061'],
        'MES': ['026', '052']
    }

    # indicator showing which groups the PIDs belong to
    subgroup_ind = dict([
        (k, pd.Index(pids).isin(v)) for k, v in subgroups.items()
    ])

    external_ref_names_de = ['GSE61794']
    external_ref_strandedness_de = ['u']
    external_ref_names_dm = ['gse38216']
    external_ref_samples_dm = ['H9 NPC 1', 'H9 NPC 2']

    external_refs_de = [
        ('GIBCO', 'NSC'),
        ('H9', 'NSC'),
    ]
    external_refs_de_labels = [t[0] for t in external_refs_de]

    external_refs_dm = [
        ('GIBCO', 'NSC'),
        ('H9', 'NSC'),
    ]
    external_refs_dm_labels = [t[0] for t in external_refs_dm]

    # plotting parameters
    cmap = 'RdYlGn_r'

    # file location parameters
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    # upset plotting colours
    subgroup_set_colours = {
        'RTK I full': '#0d680f',
        'RTK II full': '#820505',
        'MES full': '#7900ad',
        'RTK I partial': '#6ecc70',
        'RTK II partial': '#d67373',
        'MES partial': '#cc88ea',
        'mixed': '#4C72B0',
        'specific': '#f4e842',
    }

    sets_all = setops.full_partial_unique_other_sets_from_groups(pids, subgroups)
    set_colours = []
    for sg in subgroups:
        for x in ['full', 'partial']:
            k = "%s %s" % (sg, x)
            if sg in sets_all[x]:
                set_colours.append(
                    (k, {'sets': sets_all[x][sg], 'colour': subgroup_set_colours[k]})
                )
    set_colours.append(
        ('Expanded core', {'sets': sets_all['mixed'], 'colour': subgroup_set_colours['mixed']})
    )
    set_colours.append(
        ('Specific', {'sets': sets_all['specific'], 'colour': subgroup_set_colours['specific']})
    )

    ##################
    ### Strategy 1 ###
    ##################

    # some quantities relating to set membership
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))

    ss_sets = []
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        ss_sets.append(k)

    ec_sets = expanded_core_sets(setops.binary_combinations(len(pids), include_zero=False), subgroup_ind)

    pss_sets = patient_and_subgroup_specific_set(pids, subgroups)

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

    rnaseq_obj = load_rnaseq(pids, external_ref_names_de, strandedness=external_ref_strandedness_de)

    # remove unneeded samples
    idx = (~rnaseq_obj.meta.index.isin(['DURA061_NSC_N1_P5', 'DURA061_NSC_N6_P4']))
    rnaseq_obj.meta = rnaseq_obj.meta.loc[idx]
    rnaseq_obj.data = rnaseq_obj.data.loc[:, idx]
    rnaseq_obj.batch_id = rnaseq_obj.batch_id.loc[idx]

    #########################################################################
    ### STRATEGY 1: No references, just compare GBM-iNSC for each patient ###
    #########################################################################

    # only keep the syngeneic samples
    dat_s1 = rnaseq_obj.data.loc[
             :,
             rnaseq_obj.batch_id.str.contains('wtchg') & (~rnaseq_obj.data.columns.str.contains('GIBCO'))
    ]
    meta_s1 = rnaseq_obj.meta.loc[dat_s1.columns]
    groups_s1 = pd.Series(index=meta_s1.index)
    comparisons_s1 = []
    for pid in pids:
        groups_s1[groups_s1.index.str.contains('GBM') & groups_s1.index.str.contains(pid)] = "GBM%s" % pid
        groups_s1[groups_s1.index.str.contains('NSC') & groups_s1.index.str.contains(pid)] = "iNSC%s" % pid
        comparisons_s1.append(("GBM%s" % pid, "iNSC%s" % pid))

    de_res_full_s1 = de_grouped_dispersion(
        dat_s1,
        groups_s1,
        comparisons_s1,
        min_cpm=min_cpm,
        return_full=True,
        **de_params
    )
    # rename the keys to simplify
    de_res_full_s1 = dict([(pid, de_res_full_s1[("GBM%s" % pid, "iNSC%s" % pid)]) for pid in pids])
    # extract only significant DE genes
    de_res_s1 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s1.items()])

    ##################
    ### a) DE only ###
    ##################

    de_by_member = [de_res_s1[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*de_by_member)

    # add null set manually from full DE results
    de_genes_all = setops.reduce_union(*venn_set.values())
    k_null = ''.join(['0'] * len(pids))
    venn_set[k_null] = list(dat_s1.index.difference(de_genes_all))
    venn_ct[k_null] = len(venn_set[k_null])

    upset = venn.upset_set_size_plot(
        de_by_member,
        set_labels=pids,
        set_colours=set_colours,
        min_size=10,
        n_plot=30,
        default_colour='gray'
    )

    upset['figure'].savefig(os.path.join(outdir_s1, "upset_de.png"), dpi=200)
    upset['figure'].savefig(os.path.join(outdir_s1, "upset_de.tiff"), dpi=200)

    # generate wide-form lists and save to Excel file
    data = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, full_data=de_res_full_s1)
    data.to_excel(os.path.join(outdir_s1, 'full_de.xlsx'))

    # expanded core (genes DE in multiple subgroups)
    data_ec = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, include_sets=ec_sets)
    data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_de.xlsx'))

    # subgroup-specific
    data_ss = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, include_sets=ss_sets)
    data_ss.to_excel(os.path.join(outdir_s1, 'subgroup_specific_de.xlsx'))

    # patient unique
    data_pu = differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, include_sets=pu_sets)
    data_pu.to_excel(os.path.join(outdir_s1, 'patient_unique_de.xlsx'))

    # patient and subgroup-specific
    data_pss = dict([
        (pid, differential_expression.venn_set_to_dataframe(de_res_s1, venn_set, pids, include_sets=pss_sets[pid]))
        for pid in pids
    ])
    excel.pandas_to_excel(data_pss, os.path.join(outdir_s1, 'patient_and_subgroup_specific_de.xlsx'))

    ###################
    ### b) DMR only ###
    ###################

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

    # patient and subgroup-specific
    data_pss = {}
    for pid in pids:
        this_tbl = setops.venn_set_to_wide_dataframe(
                data_for_dmr_table,
                venn_set,
                pids,
                full_data=data_for_dmr_table_full,
                include_sets=pss_sets[pid],
                cols_to_include=('median_delta', 'padj'),
                consistency_check_col='median_delta',
                consistency_check_method='sign'
            )
        this_tbl.insert(0, 'gene_symbol', [t[1] for t in this_tbl.index])
        this_tbl.insert(0, 'cluster_id', [t[0] for t in this_tbl.index])
        data_pss[pid] = this_tbl
    excel.pandas_to_excel(data_pss, os.path.join(outdir_s1, 'patient_and_subgroup_specific_dmr.xlsx'))

    #############################
    ### c) Layered DE and DMR ###
    #############################

    # get the joint table
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

    # patient and subgroup-specific
    data_pss = {}
    for pid in pids:
        this_tbl =  dmr.venn_set_to_wide_dataframe(
            joint_de_dmr_concordant_s1,
            venn_set,
            pids,
            include_sets=pss_sets[pid],
            cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
            direction_col='de_logfc'
        )
        annotate_de_dmr_wide_form(this_tbl)
        data_pss[pid] = this_tbl

    excel.pandas_to_excel(data_pss, os.path.join(outdir_s1, 'patient_and_subgroup_specific_de_dmr_concordant.xlsx'))

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

    # Compute DE cross-comparison

    dat_s2 = rnaseq_obj.data
    meta_s2 = rnaseq_obj.meta
    groups_s2 = pd.Series(index=meta_s2.index)
    comparisons_s2 = {}
    for er, er_type in external_refs_de:
        groups_s2[groups_s2.index.str.contains(er) & meta_s2.loc[:, 'type'] == er_type] = er_type

    for pid in pids:
        groups_s2[groups_s2.index.str.contains('GBM') & groups_s2.index.str.contains(pid)] = "GBM%s" % pid
        groups_s2[groups_s2.index.str.contains('NSC') & groups_s2.index.str.contains(pid)] = "iNSC%s" % pid

        for pid2 in pids:
            comparisons_s2[("GBM%s" % pid, "iNSC%s" % pid2)] = "GBM%s - iNSC%s" % (pid, pid2)
        for er, er_type in external_refs_de:
            comparisons_s2[("GBM%s" % pid, er_type)] = "GBM%s - %s" % (pid, er_type)


    de_res_full_s2 = de_grouped_dispersion(
        dat_s2,
        groups_s2,
        comparisons_s2,
        min_cpm=min_cpm,
        return_full=True,
        **de_params
    )
    # rename for convenience
    de_res_full_s2 = dict([
        ((k[0].replace('GBM', ''), k[1].replace('iNSC', '')), v) for k, v in de_res_full_s2.items()
    ])

    de_res_s2 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s2.items()])

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

    #######################
    ### a) Pair-only DE ###
    #######################

    de_res_s2_idx = dict([
        (k, v.index) for k, v in de_res_s2.items()
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
        po_de_export_full[pid] = de_res_full_s2[(pid, pid)].copy()

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

        po_de_export[pid] = de_res_s2[(pid, pid)].loc[this_genes]

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

    fig, axs = plt.subplots(nrows=2, ncols=3)
    for pid in pids:
        if pid in subgroups['RTK I']:
            i = 0
            sg = subgroups['RTK I']
        else:
            i = 1
            sg = subgroups['RTK II']
        j = sg.index(pid)
        the_lists = [
            de_res_s2[(pid, r)].index for r in external_refs_de_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_de_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s vs..." % pid)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir_s2, 'number_de_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir_s2, 'number_de_multiple_references.tiff'), dpi=200)

    # plot: how many DE genes are in the pair only comparison when each reference is used?
    # NB apply correction

    # at the same time, get numbers for a bar chart about % overlap
    n_pair_only_intersect = pd.DataFrame(0, index=pids, columns=external_refs_de_labels)

    fig, axs = plt.subplots(nrows=2, ncols=3)
    for pid in pids:
        if pid in subgroups['RTK I']:
            i = 0
            sg = subgroups['RTK I']
        else:
            i = 1
            sg = subgroups['RTK II']
        j = sg.index(pid)
        the_lists = [
            set(pair_only_de.loc[pid, r]) for r in external_refs_de_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_de_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s pair only" % pid)

        for i, r in enumerate(external_refs_de_labels):
            # this will fail for anything other than 3 refs
            n_pair_only_intersect.loc[pid, r] = cts[''.join(['1'] * len(external_refs_de_labels))]
    fig.tight_layout()
    fig.savefig(os.path.join(outdir_s2, 'number_po_de_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir_s2, 'number_po_de_multiple_references.tiff'), dpi=200)

    # plot: overlap between individual references in terms of PO genes shared
    po_counts = pair_only_de.applymap(len)
    pct_pair_only_intersect = n_pair_only_intersect / po_counts.loc[:, external_refs_de_labels] * 100.

    ax = pct_pair_only_intersect.plot.bar(color=['#ff9999', '#99cc99', '#9999ff'], ec='k', legend=False)
    ax.set_xlabel('Patient')
    ax.set_ylabel("% DE genes shared")
    ax.set_ylim([0, 100])
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir_s2, "de_perc_po_gene_correspondence.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir_s2, "de_perc_po_gene_correspondence.tiff"), dpi=200)

    ########################
    ### b) Pair-only DMR ###
    ########################

    # Compute cross-comparison
    cc_dict_dmr = cross_comparison.compute_cross_comparison_correction(
        dmr_res_sign_s2,
        pids,
        external_refs=external_refs_dm_labels,
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

        # identify pair only
        this_row = pair_only_dmr.loc[pid, external_refs_dm_labels]
        this_dmr_pre = setops.reduce_intersection(*this_row)

        # no correction, to avoid cross-comparison issues
        this_dmr = this_dmr_pre

        # get the relevant paired comparison data in tabular form
        this_sign = dmr_res_s2[pid][pid].to_table(include='significant', skip_geneless=False)
        sign_cl_ids = this_sign.index

        # reduce to pair only
        this_sign = this_sign.loc[this_dmr]
        # # filter to include only TSS and Island - should be unnecessary
        # this_sign = this_sign.loc[this_sign.class_island | this_sign.class_tss]

        # expand by gene
        this_sign = expand_dmr_results_table_by_gene(this_sign, drop_geneless=True)

        this_full = dmr_res_s2[pid][pid].to_table(include='all', skip_geneless=False)
        # # filter to include only TSS and Island  - should be unnecessary
        # this_full = this_full.loc[this_full.class_island | this_full.class_tss]
        # annotate pair only
        po_col = pd.Series('N', index=this_full.index)
        po_col.loc[this_dmr] = 'Y'
        this_full.insert(this_full.shape[1], 'pair_only', po_col)
        # annotate significant (otherwise it isn't obvious)
        dm_col = pd.Series('N', index=this_full.index)
        dm_col.loc[sign_cl_ids] = 'Y'
        this_full.insert(this_full.shape[1], 'significant_dmr', dm_col)
        # expand by gene
        this_full = expand_dmr_results_table_by_gene(this_full, drop_geneless=True)

        # for export reasons, we need to convert the class_ boolean fields
        this_sign.loc[:, this_sign.columns.str.contains('class_')] = \
            this_sign.loc[:, this_sign.columns.str.contains('class_')].astype(int)
        this_full.loc[:, this_full.columns.str.contains('class_')] = \
            this_full.loc[:, this_full.columns.str.contains('class_')].astype(int)

        po_dmr_export_full[pid] = this_full
        po_dmr_export[pid] = this_sign

    excel.pandas_to_excel(po_dmr_export, os.path.join(outdir_s2, 'pair_only_dmr.xlsx'))
    excel.pandas_to_excel(po_dmr_export_full, os.path.join(outdir_s2, 'pair_only_dmr_full.xlsx'))

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

    dmr_members = {}
    for pid in pids:
        for r in external_refs_dm_labels:
            this_tbl = dmr_res_s2[pid][r].to_table(include='significant', skip_geneless=False)
            this_tbl = this_tbl.loc[this_tbl.class_island | this_tbl.class_tss]
            dmr_members[(pid, r)] = this_tbl.index

    fig, axs = plt.subplots(nrows=2, ncols=3)
    for pid in pids:
        if pid in subgroups['RTK I']:
            i = 0
            sg = subgroups['RTK I']
        else:
            i = 1
            sg = subgroups['RTK II']
        j = sg.index(pid)
        the_lists = [
            dmr_members[(pid, r)] for r in external_refs_dm_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_dm_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s vs..." % pid)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir_s2, 'number_dmr_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir_s2, 'number_dmr_multiple_references.tiff'), dpi=200)

    # plot: how many DE genes are in the pair only comparison when each reference is used?

    # at the same time, get numbers for a bar chart about % overlap
    n_pair_only_intersect = pd.DataFrame(0, index=pids, columns=external_refs_dm_labels)

    fig, axs = plt.subplots(nrows=2, ncols=3)
    for pid in pids:
        if pid in subgroups['RTK I']:
            i = 0
            sg = subgroups['RTK I']
        else:
            i = 1
            sg = subgroups['RTK II']
        j = sg.index(pid)
        the_lists = [
            set(pair_only_dmr.loc[pid, r]) for r in external_refs_dm_labels
            ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_dm_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s pair only" % pid)

        for i, r in enumerate(external_refs_dm_labels):
            n_pair_only_intersect.loc[pid, r] = cts[''.join(['1'] * len(external_refs_dm_labels))]

    fig.tight_layout()
    fig.savefig(os.path.join(outdir_s2, 'number_po_dmr_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir_s2, 'number_po_dmr_multiple_references.tiff'), dpi=200)

    # plot: overlap between individual references in terms of PO genes shared
    po_counts = pair_only_dmr.applymap(len)
    pct_pair_only_intersect = n_pair_only_intersect / po_counts.loc[:, external_refs_dm_labels] * 100.

    ax = pct_pair_only_intersect.plot.bar(color=['#ff9999', '#99cc99', '#9999ff'], ec='k', legend=False)
    ax.set_xlabel('Patient')
    ax.set_ylabel("% DMRs shared")
    ax.set_ylim([0, 100])
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir_s2, "dmr_perc_po_gene_correspondence.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir_s2, "dmr_perc_po_gene_correspondence.tiff"), dpi=200)

    ############################
    ### c) DE / DMR combined ###
    ############################

    joint_de_dmr_s2 = {}
    pair_joint = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_s1)
    for pid in pids:
        joint_de_dmr_s2[(pid, pid)] = pair_joint[pid]

    for t in zip(external_refs_de, external_refs_dm):
        de_ref = t[0][0]
        dm_ref = t[1][0]
        if de_ref != dm_ref:
            print "Warning: de_ref (%s) != dm ref (%s). Is this intended?" % (de_ref, dm_ref)
        the_dmr = dmr.DmrResultCollection(**dict([
            (pid, dmr_res_s2[pid][dm_ref]) for pid in pids
        ]))
        the_de = dict([
            (pid, de_res_s2[(pid, de_ref)]) for pid in pids
        ])
        ref_joint = rnaseq_methylationarray.compute_joint_de_dmr(the_dmr, the_de)
        for pid in pids:
            joint_de_dmr_s2[(pid, de_ref)] = ref_joint[pid]

    # filter - only DMRs with TSS or island class
    for k in joint_de_dmr_s2:
        idx = joint_de_dmr_s2[k].dmr_class_island | joint_de_dmr_s2[k].dmr_class_tss
        joint_de_dmr_s2[k] = joint_de_dmr_s2[k].loc[idx]
        print "Comparison %s. %d DE & DMR rows, %d are gene only, keeping %d." % (k, idx.size, (~idx).sum(), idx.sum())

    # filter - concordant DE/DMRs
    for k in joint_de_dmr_s2:
        idx = joint_de_dmr_s2[k].de_direction != joint_de_dmr_s2[k].dmr_direction
        joint_de_dmr_s2[k] = joint_de_dmr_s2[k].loc[idx]
        print "Comparison %s. %d DE & DMR rows, %d are concordant." % (k, idx.size, idx.sum())

    # get pair only results
    # do this using only the gene symbol

    de_dmr_by_gene = dict([
        (k, set(v.gene.values)) for k, v in joint_de_dmr_s2.items()
    ])

    pair_only_de_dmr = pd.DataFrame(index=pids, columns=external_refs_de_labels)
    for i in pair_only_de_dmr.index:
        p = de_dmr_by_gene[(i, i)]
        for j in pair_only_de_dmr.columns:
            r = de_dmr_by_gene[(i, j)]
            x, _ = setops.venn_from_arrays(p, r)
            pair_only_de_dmr.loc[i, j] = x['10']

    # take the intersection of the PO lists, then export to a file
    po_de_dmr_export = {}

    for pid in pids:
        this_row = pair_only_de_dmr.loc[pid, [t[0] for t in external_refs_de]]
        this_genes = setops.reduce_intersection(*this_row)

        # now we need to look up those genes in the relevant joint results table
        this_joint = joint_de_dmr_s2[(pid, pid)]
        po_de_dmr_export[pid] = this_joint.loc[this_joint.gene.isin(this_genes)]

    excel.pandas_to_excel(po_de_dmr_export, os.path.join(outdir_s2, "pair_only_de_dmr.xlsx"))

    # export with a different layout, analogous to strategy 1
    venn_set, venn_ct = setops.venn_from_arrays(*[po_de_dmr_export[pid].index for pid in pids])
    po_combination_export = setops.venn_set_to_wide_dataframe(
        po_de_dmr_export,
        venn_set,
        pids,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        consistency_check_col='de_logfc',
        consistency_check_method="sign"
    )

    po_combination_export.to_excel(os.path.join(outdir_s2, 'pair_only_de_dmr_wideform.xlsx'))

    # plot: number of DM-concordant DE genes in each reference comparison
    fig, axs = plt.subplots(nrows=2, ncols=3)
    for pid in pids:
        if pid in subgroups['RTK I']:
            i = 0
            sg = subgroups['RTK I']
        else:
            i = 1
            sg = subgroups['RTK II']
        j = sg.index(pid)
        the_lists = [
            de_dmr_by_gene[(pid, r)] for r in external_refs_de_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_de_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s vs..." % pid)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir_s2, 'number_de_dmr_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir_s2, 'number_de_dmr_multiple_references.tiff'), dpi=200)

    # plot: how many DM-concordant DE genes are in the pair only comparison when each reference is used?

    # at the same time, get numbers for a bar chart about % overlap
    n_pair_only_intersect = pd.DataFrame(0, index=pids, columns=external_refs_de_labels)

    fig, axs = plt.subplots(nrows=2, ncols=3)
    for pid in pids:
        if pid in subgroups['RTK I']:
            i = 0
            sg = subgroups['RTK I']
        else:
            i = 1
            sg = subgroups['RTK II']
        j = sg.index(pid)
        the_lists = [
            set(pair_only_de_dmr.loc[pid, r]) for r in external_refs_de_labels
            ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_de_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s pair only" % pid)

        for i, r in enumerate(external_refs_de_labels):
            n_pair_only_intersect.loc[pid, r] = cts[''.join(['1'] * len(external_refs_de_labels))]

    fig.tight_layout()
    fig.savefig(os.path.join(outdir_s2, 'number_po_de_dmr_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir_s2, 'number_po_de_dmr_multiple_references.tiff'), dpi=200)

    # plot: overlap between individual references in terms of PO genes shared
    po_counts = pair_only_de_dmr.applymap(len)
    pct_pair_only_intersect = n_pair_only_intersect / po_counts.loc[:, external_refs_de_labels] * 100.

    ax = pct_pair_only_intersect.plot.bar(color=['#ff9999', '#99cc99', '#9999ff'], ec='k', legend=False)
    ax.set_xlabel('Patient')
    ax.set_ylabel("% DE genes shared")
    ax.set_ylim([0, 100])
    ax.figure.tight_layout()
    ax.figure.savefig(os.path.join(outdir_s2, "de_dmr_perc_po_gene_correspondence.png"), dpi=200)
    ax.figure.savefig(os.path.join(outdir_s2, "de_dmr_perc_po_gene_correspondence.tiff"), dpi=200)
