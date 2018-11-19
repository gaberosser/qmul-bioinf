import multiprocessing as mp
import os
import pandas as pd
import numpy as np
import pickle
from utils import output, setops, excel, log
from methylation import dmr, process, loader as methylation_loader, annotation_gene_to_ensembl
from rnaseq import loader as rnaseq_loader, differential_expression, general, filter
from integrator import rnaseq_methylationarray
from analysis import cross_comparison
from load_data import loader
from plotting import venn
from matplotlib import pyplot as plt, text, patches
import seaborn as sns
import consts


logger = log.get_console_logger()


def load_methylation(pids, ref_names=None, norm_method='swan', ref_name_filter=None, patient_samples=None):
    """
    Load and prepare the Illumina methylation data
    """
    # patient data
    obj = methylation_loader.load_by_patient(pids, norm_method=norm_method, samples=patient_samples)
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
    # dmr.add_merged_probe_classes(anno)
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
    """
    Compute DMRs for paired GBM-iNSC comparisons (defined in that order) for all patients
    :param me_data:
    :param me_meta:
    :param anno:
    :param pids:
    :param dmr_params:
    :return:
    """
    dmr_res_obj = dmr.DmrResults(anno=anno)
    dmr_res_obj.identify_clusters(**dmr_params)
    dmr_res = {}

    for pid in pids:
        the_idx1 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == 'GBM')
        the_idx2 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == 'iNSC')
        # control comparison order
        the_samples = [
            me_meta.index[the_idx1],
            me_meta.index[the_idx2],
        ]

        # the_idx = the_idx1 | the_idx2
        # the_groups = me_meta.loc[the_idx, 'type'].values
        # the_samples = me_meta.index[the_idx].groupby(the_groups).values()

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


def de_results_hash(pids, de_params):
    hash_elements = tuple(sorted(pids)) + tuple([de_params[k] for k in sorted(de_params.keys())])
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


def expand_dmr_results_table_by_gene_and_cluster(df, drop_geneless=False):
    df.insert(0, 'cluster_id', df.index)

    if not drop_geneless:
        # we need to protect the geneless results by adding dummy values
        g_count = df.genes.apply(len)
        df.loc[g_count == 0, 'genes'] = [(('!', '!'),)] * (g_count == 0).sum()

    # open out the the (gene, relation) lists contained in the `genes` column
    g = pd.DataFrame(
        df.genes.apply(pd.Series).stack().reset_index(level=-1, drop=True),
        columns=['gene_relation']
    )
    g.insert(0, 'relation', g.gene_relation.apply(lambda x: x[1]))
    g.insert(0, 'gene', g.gene_relation.apply(lambda x: x[0]))
    g = g.drop('gene_relation', axis=1).reset_index()

    # pivot table around the relation
    tbl = g.pivot_table(index=['cluster_id', 'gene'], columns='relation', fill_value=0, aggfunc='size').reset_index()

    # extract other attributes from original dataframe
    attrs = df.drop(['cluster_id', 'genes'], axis=1).loc[tbl.cluster_id]

    # combine them
    res = pd.concat((tbl.set_index('cluster_id'), attrs), axis=1).reset_index()

    # if maintaining geneless results, deprotect these now
    if not drop_geneless:
        new_gene_col = res.loc[:, 'gene'].replace('!', np.nan)  # don't use None as the replacement variable!
        res.loc[:, 'gene'] = new_gene_col
        # drop the pivot column, but only if it actually exists
        if '!' in res.columns:
            res.drop('!', axis=1, inplace=True)
        # quick sanity check
        cols = sorted(set(g.relation.unique()).difference(['!']))
        ix = res.gene.isnull()
        if not (res.loc[ix, cols].sum(axis=1) == 0).all():
            logger.error("Some 'geneless' clusters still have gene relations defined.")
            raise ValueError("Some 'geneless' clusters still have gene relations defined.")

    res.index = zip(res.cluster_id, res.gene)
    return res


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


def partial_subgroup_specific(pids, subgroup_ind):
    """
    For each subgroup, get the list of sets that corresponds to >= 2 members of that subgroup and no members of any
    other.
    :param pids: List of PIDs
    :param subgroup_ind: Dict, keyed by subgroup name. Values are np.array of booleans (i.e. boolean indexers). Order
    of indexers must match pids.
    :return: Dict, keyed by subgroup name. Each value is a list of sets.
    """
    ss_part_sets = {}
    candidates = list(setops.binary_combinations_sum_gte(len(pids), 2))
    for grp in subgroup_ind:
        ss_part_sets[grp] = [t for t in candidates if np.array([x for x in t]).astype(int)[~subgroup_ind[grp]].sum() == 0]
    return ss_part_sets


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


def de_export_to_ipa(de_wideform, pids):
    cols = []
    for p in pids:
        cols.append("%s logFC" % p)
        cols.append("%s FDR" % p)
    ipa_export = pd.DataFrame(index=de_wideform.index, columns=cols)
    for pid in pids:
        ix = de_wideform.loc[:, pid] == 'Y'
        ipa_export.loc[ix, "%s logFC" % pid] = de_wideform.loc[ix, "%s_logFC" % pid]
        ipa_export.loc[ix, "%s FDR" % pid] = de_wideform.loc[ix, "%s_FDR" % pid]
    return ipa_export


def export_dmr_data_for_ipa(results_tbl):
    # output data for IPA
    all_genes = set()
    for v in results_tbl.values():
        all_genes.update(v.gene.dropna())
    # convert to Ensembl and drop any that cannot be found
    genes_for_ipa = annotation_gene_to_ensembl.gene_to_ens(all_genes).dropna().sort_values()
    # data import is made easier if the columns contain letters as well as numbers
    ipa_export = pd.DataFrame(index=genes_for_ipa.index, columns=["Patient %s" % p for p in results_tbl.keys()])

    for pid, this in results_tbl.items():
        ix = this.gene.isin(genes_for_ipa.index)
        this = this.loc[ix]
        this = this.loc[~this.padj.isnull()]
        this = this.groupby('gene').padj.min()
        ipa_export.loc[this.index, "Patient %s" % pid] = this
    # easier to leave non-significant values as NaN (missing)
    # alternative: set p value to 1.0
    # ipa_export = ipa_export.fillna(1.)

    # replace the index with the Ensembl IDs
    ipa_export = ipa_export.set_index(genes_for_ipa)
    return ipa_export


if __name__ == "__main__":
    outdir = output.unique_output_dir("hgic_de_dmr_two_strategies", reuse_empty=True)
    outdir_s1 = os.path.join(outdir, 's1')
    outdir_s2 = os.path.join(outdir, 's2')
    if not os.path.isdir(outdir_s1):
        os.makedirs(outdir_s1)
    if not os.path.isdir(outdir_s2):
        os.makedirs(outdir_s2)

    outdir_s1_ipa = os.path.join(outdir_s1, 'ipa')
    outdir_s2_ipa = os.path.join(outdir_s2, 'ipa')
    if not os.path.isdir(outdir_s1_ipa):
        os.makedirs(outdir_s1_ipa)
    if not os.path.isdir(outdir_s2_ipa):
        os.makedirs(outdir_s2_ipa)

    min_cpm = 1.
    # Set this to True to include reference methylation data
    # This will limit the number of available probes (to 450K)
    include_external_dm_refs = False

    de_params = consts.DE_PARAMS
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    norm_method_s1 = 'swan'

    pids = consts.PIDS

    subgroups = consts.SUBGROUPS

    subgroups_lookup = {}
    for grp, arr in subgroups.items():
        subgroups_lookup.update(dict([
            (t, grp) for t in arr
        ]))

    # indicator showing which groups the PIDs belong to
    subgroup_ind = dict([
        (k, pd.Index(pids).isin(v)) for k, v in subgroups.items()
    ])

    external_ref_names_de = ['GSE61794']
    external_ref_strandedness_de = ['u']

    if include_external_dm_refs:
        external_ref_names_dm = ['gse38216']
        external_ref_samples_dm = ['H9 NPC 1', 'H9 NPC 2']
    else:
        external_ref_names_dm = None
        external_ref_samples_dm = None

    external_refs_de = [
        ('GIBCO', 'NSC'),
        ('H9', 'NSC'),
    ]
    external_refs_de_labels = [t[0] for t in external_refs_de]

    if include_external_dm_refs:
        external_refs_dm = [
            ('GIBCO', 'NSC'),
            ('H9', 'NSC'),
        ]
    else:
        external_refs_dm = [
            ('GIBCO', 'NSC'),
        ]

    external_refs_dm_labels = [t[0] for t in external_refs_dm]

    # plotting parameters
    cmap = 'RdYlGn_r'

    # file location parameters
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')
    DE_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'de')

    if not os.path.isdir(DE_LOAD_DIR):
        logger.info("Created DE results save path %s.", DE_LOAD_DIR)
        os.makedirs(DE_LOAD_DIR)

    if not os.path.isdir(DMR_LOAD_DIR):
        logger.info("Created DMR results save path %s.", DMR_LOAD_DIR)
        os.makedirs(DMR_LOAD_DIR)

    # upset plotting colours
    subgroup_set_colours = {
        'RTK I full': '#0d680f',
        'RTK II full': '#820505',
        'MES full': '#7900ad',
        'RTK I partial': '#6ecc70',
        'RTK II partial': '#d67373',
        'MES partial': '#cc88ea',
        'Expanded core': '#4C72B0',
        'Specific': '#f4e842',
    }

    sets_all = setops.full_partial_unique_other_sets_from_groups(pids, subgroups)

    ##################
    ### Strategy 1 ###
    ##################

    # some quantities relating to set membership
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))

    ss_sets = []
    for grp in subgroup_ind:
        k = ''.join(subgroup_ind[grp].astype(int).astype(str))
        ss_sets.append(k)

    ss_part_sets = partial_subgroup_specific(pids, subgroup_ind)

    ec_sets = expanded_core_sets(setops.binary_combinations(len(pids), include_zero=False), subgroup_ind)

    pss_sets = patient_and_subgroup_specific_set(pids, subgroups)

    # load methylation
    # For S1, we just need the paired comparison. This is important - any other samples lead to a change in the final
    # probe list (any NA rows are dropped from ALL samples).

    me_obj, anno = load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
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
        logger.info("Reading S1 DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        dmr_res_s1 = paired_dmr(me_data, me_meta, anno, pids, dmr_params)
        # Save DMR results to disk
        dmr_res_s1.to_pickle(fn, include_annotation=False)
        print "Saved DMR results to %s" % fn

    # extract results
    dmr_res_full_s1 = dmr_res_s1.results
    dmr_res_sign_s1 = dmr_res_s1.results_significant

    rnaseq_obj = load_rnaseq(
        pids,
        external_ref_names_de,
        strandedness=external_ref_strandedness_de,
    )

    rnaseq_obj.filter_by_sample_name(consts.ALL_RNASEQ_SAMPLES)

    #########################################################################
    ### STRATEGY 1: No references, just compare GBM-iNSC for each patient ###
    #########################################################################

    # only keep the required syngeneic samples for this analysis
    dat_s1 = rnaseq_obj.data.loc[
                :, rnaseq_obj.meta.index.isin(consts.S1_RNASEQ_SAMPLES)
             ]
    meta_s1 = rnaseq_obj.meta.loc[dat_s1.columns]

    the_hash = de_results_hash(meta_s1.index.tolist(), de_params)
    filename = 'de_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DE results from %s", fn)
        with open(fn, 'rb') as f:
            de_res_full_s1 = pickle.load(f)
    else:
        groups_s1 = pd.Series(index=meta_s1.index)
        comparisons_s1 = {}
        for pid in pids:
            groups_s1[groups_s1.index.str.contains('GBM') & groups_s1.index.str.contains(pid)] = "GBM%s" % pid
            groups_s1[groups_s1.index.str.contains('NSC') & groups_s1.index.str.contains(pid)] = "iNSC%s" % pid

            comparisons_s1[("GBM%s" % pid, "iNSC%s" % pid)] = "GBM%s - iNSC%s" % (pid, pid)

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

        with open(fn, 'wb') as f:
            pickle.dump(de_res_full_s1, f)

        logger.info("Saved S1 DE results to %s", fn)

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

    upset = venn.upset_plot_with_groups(
        de_by_member,
        pids,
        subgroup_ind,
        subgroup_set_colours,
        min_size=10,
        n_plot=30,
        default_colour='gray'
    )

    upset['figure'].savefig(os.path.join(outdir_s1, "upset_de.png"), dpi=200)
    upset['figure'].savefig(os.path.join(outdir_s1, "upset_de.tiff"), dpi=200)

    # generate wide-form lists and save to Excel file
    data = setops.venn_set_to_wide_dataframe(
        de_res_s1,
        venn_set,
        pids,
        full_data=de_res_full_s1,
        cols_to_include=['logFC', 'FDR'],
        consistency_check_col='logFC',
        consistency_check_method='sign'
    )
    # add gene symbols back in
    general.add_gene_symbols_to_ensembl_data(data)
    data.to_excel(os.path.join(outdir_s1, 'full_de.xlsx'))
    de_export_to_ipa(data, pids).to_excel(os.path.join(outdir_s1_ipa, "full_de_for_ipa.xlsx"))

    # expanded core (genes DE in multiple subgroups)
    data_ec = setops.venn_set_to_wide_dataframe(
        de_res_s1,
        venn_set,
        pids,
        include_sets=ec_sets,
        full_data=de_res_full_s1,
        cols_to_include=['logFC', 'FDR'],
        static_cols_to_include=['Gene Symbol'],
        consistency_check_col='logFC',
        consistency_check_method='sign'
    )
    data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_de.xlsx'))
    de_export_to_ipa(data_ec, pids).to_excel(os.path.join(outdir_s1_ipa, "expanded_core_de_for_ipa.xlsx"))

    # subgroup-specific (full - must be in all patients)
    # split: one subgroup per tab
    data_ss = {}
    for sg in subgroups:
        k = ''.join(subgroup_ind[sg].astype(int).astype(str))
        data_ss[sg] = setops.venn_set_to_wide_dataframe(
            de_res_s1,
            venn_set,
            pids,
            include_sets=[k],
            full_data=de_res_full_s1,
            cols_to_include=['logFC', 'FDR'],
            static_cols_to_include=['Gene Symbol'],
            consistency_check_col='logFC',
            consistency_check_method='sign'
        )
    excel.pandas_to_excel(data_ss, os.path.join(outdir_s1, 'full_subgroup_specific_de.xlsx'))

    # subgroup-specific (partial: at least 2 members)
    # split: one subgroup per tab
    data_ss_part = {}
    for sg in subgroups:
        data_ss_part[sg] = setops.venn_set_to_wide_dataframe(
            de_res_s1,
            venn_set,
            pids,
            include_sets=ss_part_sets[sg],
            full_data=de_res_full_s1,
            cols_to_include=['logFC', 'FDR'],
            static_cols_to_include=['Gene Symbol'],
            consistency_check_col='logFC',
            consistency_check_method='sign'
        )
    excel.pandas_to_excel(data_ss_part, os.path.join(outdir_s1, 'partial_subgroup_specific_de.xlsx'))

    # patient specific
    data_pu = setops.venn_set_to_wide_dataframe(
        de_res_s1,
        venn_set,
        pids,
        include_sets=pu_sets,
        full_data=de_res_full_s1,
        cols_to_include=['logFC', 'FDR'],
        static_cols_to_include=['Gene Symbol'],
        consistency_check_col='logFC',
        consistency_check_method='sign'
    )
    data_pu.to_excel(os.path.join(outdir_s1, 'patient_specific_de.xlsx'))
    de_export_to_ipa(data_pu, pids).to_excel(os.path.join(outdir_s1_ipa, "patient_specific_de_for_ipa.xlsx"))

    # patient or subgroup-specific
    data_pss = {}
    export_pss = pd.DataFrame(index=de_res_full_s1.values()[0].index)
    for pid in pids:
        data_pss[pid] = setops.venn_set_to_wide_dataframe(
            de_res_s1,
            venn_set,
            pids,
            include_sets=pss_sets[pid],
            full_data=de_res_full_s1,
            cols_to_include=['logFC', 'FDR'],
            static_cols_to_include=['Gene Symbol'],
            consistency_check_col='logFC',
            consistency_check_method='sign'
        )
        this_export = de_export_to_ipa(data_pss[pid], pids)
        export_pss.insert(0, "%s FDR" % pid, this_export.loc[:, "%s FDR" % pid])
        export_pss.insert(0, "%s logFC" % pid, this_export.loc[:, "%s logFC" % pid])
    excel.pandas_to_excel(data_pss, os.path.join(outdir_s1, 'patient_or_subgroup_specific_de.xlsx'))
    export_pss.to_excel(os.path.join(outdir_s1_ipa, "patient_or_subgroup_specific_de_for_ipa.xlsx"))

    ###################
    ### b) DMR only ###
    ###################

    dmr_by_member = [dmr_res_sign_s1[pid].keys() for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*dmr_by_member)

    # add null set manually from full DMR results
    dmr_id_all = setops.reduce_union(*venn_set.values())
    k_null = ''.join(['0'] * len(pids))
    venn_set[k_null] = list(set(dmr_res_full_s1[pids[0]].keys()).difference(dmr_id_all))
    venn_ct[k_null] = len(venn_set[k_null])

    upset = venn.upset_plot_with_groups(
        dmr_by_member,
        pids,
        subgroup_ind,
        subgroup_set_colours,
        min_size=10,
        n_plot=30,
        default_colour='gray'
    )
    upset['axes']['main'].set_ylabel('Number of DMRs in set')
    upset['axes']['set_size'].set_xlabel('Number of DMRs in single comparison')

    upset['figure'].savefig(os.path.join(outdir_s1, "upset_dmr.png"), dpi=200)
    upset['figure'].savefig(os.path.join(outdir_s1, "upset_dmr.tiff"), dpi=200)

    # generate wide-form lists and save to Excel file
    data_for_dmr_table = {}
    data_for_dmr_table_full = {}
    for pid in pids:
        this_sign = dmr_res_s1[pid].to_table(include='significant', skip_geneless=True)
        data_for_dmr_table[pid] = expand_dmr_results_table_by_gene_and_cluster(this_sign, drop_geneless=True)

        this_full = dmr_res_s1[pid].to_table(include='all',skip_geneless=True)
        data_for_dmr_table_full[pid] = expand_dmr_results_table_by_gene_and_cluster(this_full, drop_geneless=True)

    # output data for IPA
    export_dmr_data_for_ipa(data_for_dmr_table_full).to_excel(os.path.join(outdir_s1_ipa, "full_dmr_genes_for_ipa.xlsx"))

    # recalculate venn set
    dmr_by_member = [data_for_dmr_table[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*dmr_by_member)

    # add null set manually from full DMR results
    dmr_id_all = setops.reduce_union(*venn_set.values())
    k_null = ''.join(['0'] * len(pids))
    venn_set[k_null] = list(data_for_dmr_table_full[pids[0]].index.difference(dmr_id_all))
    venn_ct[k_null] = len(venn_set[k_null])

    # get all gene relation choices
    gene_relation_choices = set()
    for pp in dmr_res_s1.clusters.values():
        gene_relation_choices.update([t[1] for t in list(pp.genes)])

    data = setops.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        cols_to_include=('median_delta', 'padj'),
        static_cols_to_include=['cluster_id', 'gene'] + sorted(gene_relation_choices),
        consistency_check_col='median_delta',
        consistency_check_method='sign'
    )

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
        static_cols_to_include=['cluster_id', 'gene'] + sorted(gene_relation_choices),
        consistency_check_col='median_delta',
        consistency_check_method='sign'
    )

    data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_dmr.xlsx'))

    # subgroup-specific
    data_ss = setops.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        include_sets=ss_sets,
        cols_to_include=('median_delta', 'padj'),
        static_cols_to_include=['cluster_id', 'gene'] + sorted(gene_relation_choices),
        consistency_check_col='median_delta',
        consistency_check_method='sign'
    )

    data_ss.to_excel(os.path.join(outdir_s1, 'subgroup_specific_dmr.xlsx'))

    # patient unique
    data_pu = setops.venn_set_to_wide_dataframe(
        data_for_dmr_table,
        venn_set,
        pids,
        full_data=data_for_dmr_table_full,
        include_sets=pu_sets,
        cols_to_include=('median_delta', 'padj'),
        static_cols_to_include=['cluster_id', 'gene'] + sorted(gene_relation_choices),
        consistency_check_col='median_delta',
        consistency_check_method='sign'
    )

    data_pu.to_excel(os.path.join(outdir_s1, 'patient_specific_dmr.xlsx'))

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
            static_cols_to_include=['cluster_id', 'gene'] + sorted(gene_relation_choices),
            consistency_check_col='median_delta',
            consistency_check_method='sign'
        )

        data_pss[pid] = this_tbl
    excel.pandas_to_excel(data_pss, os.path.join(outdir_s1, 'patient_or_subgroup_specific_dmr.xlsx'))

    #############################
    ### c) Layered DE and DMR ###
    #############################

    # This has been moved to two_strategies_combine_de_dmr

    ##################
    ### STRATEGY 2 ###
    ##################

    norm_method_s2 = 'pbc'

    # load methylation with external references
    me_obj_with_ref, anno = load_methylation(
        pids,
        ref_names=external_ref_names_dm,
        ref_name_filter=external_ref_samples_dm,
        norm_method=norm_method_s2,
        patient_samples=consts.ALL_METHYL_SAMPLES  # technically contains H9 too, but these will be ignored
    )
    me_data_with_ref = me_obj_with_ref.data

    # Compute DE cross-comparison or load if results already exist

    dat_s2 = rnaseq_obj.data
    meta_s2 = rnaseq_obj.meta

    the_hash = de_results_hash(meta_s2.index.tolist(), de_params)
    filename = 'de_results_cross_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S2 DE results from %s", fn)
        with open(fn, 'rb') as f:
            de_res_full_s2 = pickle.load(f)
    else:
        groups_s2 = pd.Series(index=meta_s2.index)
        comparisons_s2 = {}
        for er, er_type in external_refs_de:
            groups_s2[groups_s2.index.str.contains(er) & (meta_s2.loc[:, 'type'] == er_type)] = er

        for pid in pids:
            groups_s2[groups_s2.index.str.contains('GBM') & groups_s2.index.str.contains(pid)] = "GBM%s" % pid
            groups_s2[groups_s2.index.str.contains('NSC') & groups_s2.index.str.contains(pid)] = "iNSC%s" % pid

            for pid2 in pids:
                comparisons_s2[("GBM%s" % pid, "iNSC%s" % pid2)] = "GBM%s - iNSC%s" % (pid, pid2)
            for er, er_type in external_refs_de:
                comparisons_s2[("GBM%s" % pid, er)] = "GBM%s - %s" % (pid, er)

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

        with open(fn, 'wb') as f:
            pickle.dump(de_res_full_s2, f)

        logger.info("Saved S2 DE results to %s", fn)

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

    dmr_res_full_s2 = {}
    dmr_res_sign_s2 = {}

    for k1 in dmr_res_s2.keys():
        for k2 in dmr_res_s2[k1].keys():
            dmr_res_sign_s2[(k1, k2)] = dmr_res_s2[k1][k2].results_significant
            dmr_res_full_s2[(k1, k2)] = dmr_res_s2[k1][k2].results

    #######################
    ### a) Pair-only DE ###
    #######################

    # export full syngeneic and reference comparisons (not cross comparisons) to Excel
    for_export = {}
    for pid in pids:
        this_full = {}
        this_export = {}
        export_order = []
        this_export["%s_syngeneic" % pid] = de_res_s2[(pid, pid)]
        this_full["%s_syngeneic" % pid] = de_res_full_s2[(pid, pid)]
        export_order.append("%s_syngeneic" % pid)
        for r in external_refs_de_labels:
            this_export["%s_%s" % (pid, r)] = de_res_s2[(pid, r)]
            this_full["%s_%s" % (pid, r)] = de_res_full_s2[(pid, r)]
            export_order.append("%s_%s" % (pid, r))

        venn_set, _ = setops.venn_from_arrays(*[this_export[k].index for k in export_order])

        for_export[pid] = differential_expression.venn_set_to_dataframe(
            this_export,
            venn_set,
            export_order,
            full_data=this_full,
            add_null_set=True
        )
    excel.pandas_to_excel(for_export, os.path.join(outdir_s2, "full_de_by_patient.xlsx"))

    # wideform version of this (i.e. 30 blocks)
    # we can't use the Venn approach here, but we don't need to
    all_genes = sorted(setops.reduce_union(*[t.index for t in de_res_full_s2.values()]))
    for_export = pd.DataFrame(index=all_genes)
    member_cols = []
    general.add_gene_symbols_to_ensembl_data(for_export)
    for pid in pids:
        for r in [pid] + external_refs_de_labels:
            this_sign = de_res_s2[(pid, r)]
            this_full = de_res_full_s2[(pid, r)]
            this_yn = pd.Series('N', index=all_genes)
            this_yn.loc[this_sign.index] = 'Y'
            k = "%s_%s" % (pid, r if r != pid else 'syngeneic')
            member_cols.append(k)
            for_export.insert(
                for_export.shape[1],
                k,
                this_yn
            )
            for_export.insert(
                for_export.shape[1],
                "%s_logFC" % k,
                this_full.reindex(all_genes)['logFC']
            )
            for_export.insert(
                for_export.shape[1],
                "%s_FDR" % k,
                this_full.reindex(all_genes)['FDR']
            )

    # consistency column: is this useful with this many samples?
    all_yn = for_export.loc[:, member_cols]
    all_logFC = for_export.loc[:, ["%s_logFC" % t for t in member_cols]].values
    all_logFC[(all_yn == 'N').values] = np.nan
    row_sum_abs = np.abs(np.nansum(np.sign(all_logFC), axis=1))
    row_nz = (~np.isnan(np.sign(all_logFC))).sum(axis=1)

    consist = pd.Series('N', index=all_genes)
    consist.loc[row_sum_abs == row_nz] = 'Y'

    for_export.insert(
        for_export.shape[1],
        'consistent',
        consist
    )
    for_export.to_excel(os.path.join(outdir_s2, "full_de.xlsx"))

    # export for IPA
    # NB: we still need the syngeneic version, because it isn't *quite* the same as S1 (due to lumped dispersion)
    ipa_export = for_export.drop(['Gene Symbol', 'consistent'], axis=1)

    # since the maximum number of observations is 20, split this over 3 files
    ix = ipa_export.columns.str.contains('_syngeneic') & ipa_export.columns.str.contains(r'(_logFC)|(_FDR)')
    ipa_export.loc[:, ix].to_excel(os.path.join(outdir_s2_ipa, "full_de_syngeneic.xlsx"))
    ix = ipa_export.columns.str.contains('_H9') & ipa_export.columns.str.contains(r'(_logFC)|(_FDR)')
    ipa_export.loc[:, ix].to_excel(os.path.join(outdir_s2_ipa, "full_de_h9.xlsx"))
    ix = ipa_export.columns.str.contains('_GIBCO') & ipa_export.columns.str.contains(r'(_logFC)|(_FDR)')
    ipa_export.loc[:, ix].to_excel(os.path.join(outdir_s2_ipa, "full_de_gibco.xlsx"))

    # Find DE genes that are "pair only" in ALL cross comparisons
    # These are not syngeneic-specific, but are "our iNSC specific", suggesting they may be an artefact?

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
    # po_de_export_full = {}
    for pid in pids:
        # po_de_export_full[pid] = de_res_full_s2[(pid, pid)].copy()

        this_row = pair_only_de.loc[pid, external_refs_de_labels]
        this_genes_pre = setops.reduce_intersection(*this_row)

        # can correct as follows
        # this_genes = sorted(this_genes_pre.difference(po_specific_to_all_refs))
        # print "PID %s. Subtracted %d correction genes from the %d PO intersection genes to leave %d PO genes" % (
        #     pid, len(po_specific_to_all_refs), len(this_genes_pre), len(this_genes)
        # )

        # we won't, to avoid cross-comparison issues
        this_genes = this_genes_pre

        po_col = pd.Series('N', index=de_res_full_s2[(pid, pid)].index)
        po_col.loc[this_genes] = 'Y'
        po_de_export[pid] = de_res_s2[(pid, pid)].loc[this_genes]

    # export
    excel.pandas_to_excel(po_de_export, os.path.join(outdir_s2, "pair_only_de.xlsx"))

    # export with a different layout, analogous to strategy 1
    # this is helpful for comparing BETWEEN patients
    venn_set, venn_ct = setops.venn_from_arrays(*[po_de_export[pid].index for pid in pids])
    po_combination_export = differential_expression.venn_set_to_dataframe(po_de_export, venn_set, pids)

    po_combination_export = setops.venn_set_to_wide_dataframe(
        po_de_export,
        venn_set,
        pids,
        cols_to_include=('logFC', 'FDR'),
        static_cols_to_include=['Gene Symbol'],
        consistency_check_col='logFC',
        consistency_check_method="sign"
    )
    po_combination_export.to_excel(os.path.join(outdir_s2, 'pair_only_de_inter_patient.xlsx'))

    # can also use this representation to export to IPA
    de_export_to_ipa(po_combination_export, pids).to_excel(os.path.join(outdir_s2_ipa, "pair_only_de_for_ipa.xlsx"))


    # array of Venn plots: GBM vs X (# DE genes)
    # only possible if number of references <= 3 (otherwise too many sets)
    if len(external_refs_de_labels) < 4:
        nrows = len(subgroups)
        ncols = max([len(t) for t in subgroups.values()])
        set_labels = ['Syngeneic iNSC'] + external_refs_de_labels
        set_colours = ['r', 'g', 'b']
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
        # assume we'll have one axis leftover!
        # TODO: if not, we'll need to add one or move the legend outside the main figure?
        cax = axs[-1, -1]
        ax_set = set(axs.flat)
        for pid in pids:
            sg = subgroups_lookup[pid]
            sg_members = subgroups[sg]
            i = subgroups.keys().index(sg)
            j = sg_members.index(pid)
            the_lists = [
                de_res_s2[(pid, r)].index for r in [pid] + external_refs_de_labels
            ]
            venn_sets, cts = setops.venn_from_arrays(*the_lists)
            venn.venn_diagram(*the_lists, set_labels=None, ax=axs[i, j])
            axs[i, j].set_title("GBM%s" % pid, y=0.95)
            ax_set.remove(axs[i, j])
            # resize text
            h_txt = [t for t in axs[i, j].get_children() if isinstance(t, text.Text)]
            plt.setp(h_txt, fontsize=12)
        for ax in ax_set:
            ax.axis('off')

        # legend
        leg_objs = [
            patches.Rectangle([0, 0], 1, 1, color=c, alpha=0.4, label=l) for c, l in zip(set_colours, set_labels)
        ]
        cax.legend(handles=leg_objs, loc='center')

        fig.subplots_adjust(left=0., right=1., bottom=0., top=.93, wspace=0., hspace=0.05)

        fig.savefig(os.path.join(outdir_s2, 'venn_number_de_all_nsc.png'), dpi=200)
        fig.savefig(os.path.join(outdir_s2, 'venn_number_de_all_nsc.tiff'), dpi=200)

    # array of Venn plots: GBM vs Ref X (# DE genes)

    nrows = len(subgroups)
    ncols = max([len(t) for t in subgroups.values()])
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
    ax_set = set(axs.flat)
    for pid in pids:
        sg = subgroups_lookup[pid]
        sg_members = subgroups[sg]
        i = subgroups.keys().index(sg)
        j = sg_members.index(pid)
        the_lists = [
            de_res_s2[(pid, r)].index for r in external_refs_de_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_de_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s vs..." % pid)
        ax_set.remove(axs[i, j])
    fig.subplots_adjust(left=0.02, right=1., bottom=0.05, top=.95, wspace=0.1, hspace=0.3)
    for ax in ax_set:
        ax.set_visible(False)
    fig.savefig(os.path.join(outdir_s2, 'number_de_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir_s2, 'number_de_multiple_references.tiff'), dpi=200)

    # array of Venn plots: GBM vs Ref X (# pair only DE genes)

    # at the same time, get numbers for a bar chart about % overlap
    n_pair_only_intersect = pd.DataFrame(0, index=pids, columns=external_refs_de_labels)

    nrows = len(subgroups)
    ncols = max([len(t) for t in subgroups.values()])
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
    ax_set = set(axs.flat)
    for pid in pids:
        sg = subgroups_lookup[pid]
        sg_members = subgroups[sg]
        i = subgroups.keys().index(sg)
        j = sg_members.index(pid)
        the_lists = [
            set(pair_only_de.loc[pid, r]) for r in external_refs_de_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_de_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s pair only" % pid)
        ax_set.remove(axs[i, j])

        for i, r in enumerate(external_refs_de_labels):
            n_pair_only_intersect.loc[pid, r] = cts[''.join(['1'] * len(external_refs_de_labels))]

    fig.tight_layout()
    for ax in ax_set:
        ax.set_visible(False)
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

    # export full syngeneic and reference comparisons (not cross comparisons) to Excel
    for_export = {}
    for_ipa = {'syngeneic': {}}
    for r in external_refs_dm_labels:
        for_ipa[r] = {}

    for pid in pids:
        this_full = {}
        this_export = {}
        export_order = []

        this_export["%s_syngeneic" % pid] = dmr_res_s2[pid][pid].to_table(include='significant', skip_geneless=False)
        this_export["%s_syngeneic" % pid] = expand_dmr_results_table_by_gene_and_cluster(
            this_export["%s_syngeneic" % pid],
            drop_geneless=False
        )
        for_ipa['syngeneic']["%s_syngeneic" % pid] = this_export["%s_syngeneic" % pid]

        this_full["%s_syngeneic" % pid] = dmr_res_s2[pid][pid].to_table(include='all', skip_geneless=False)
        this_full["%s_syngeneic" % pid] = expand_dmr_results_table_by_gene_and_cluster(
            this_full["%s_syngeneic" % pid],
            drop_geneless=False
        )

        export_order.append("%s_syngeneic" % pid)
        for r in external_refs_dm_labels:
            this_export["%s_%s" % (pid, r)] = dmr_res_s2[pid][r].to_table(include='significant', skip_geneless=False)
            this_export["%s_%s" % (pid, r)] = expand_dmr_results_table_by_gene_and_cluster(
                this_export["%s_%s" % (pid, r)],
                drop_geneless=False
            )
            for_ipa[r]["%s_%s" % (pid, r)] = this_export["%s_%s" % (pid, r)]

            this_full["%s_%s" % (pid, r)] = dmr_res_s2[pid][r].to_table(include='all', skip_geneless=False)
            this_full["%s_%s" % (pid, r)] = expand_dmr_results_table_by_gene_and_cluster(
                this_full["%s_%s" % (pid, r)],
                drop_geneless=False
            )
            export_order.append("%s_%s" % (pid, r))

        venn_set, _ = setops.venn_from_arrays(*[this_export[k].index for k in export_order])

        for_export[pid] = setops.venn_set_to_wide_dataframe(
            this_export,
            venn_set,
            export_order,
            full_data=this_full,
            cols_to_include=('median_delta', 'padj'),
            static_cols_to_include=['cluster_id', 'gene'] + sorted(gene_relation_choices),
            consistency_check_col='median_delta',
            consistency_check_method="sign"
        )
    excel.pandas_to_excel(for_export, os.path.join(outdir_s2, "full_dmr.xlsx"))

    # Export data for IPA analysis TODO: check
    # split into 3 x 10 to allow uploading (max: 20)
    for k, v in for_ipa.items():
        export_dmr_data_for_ipa(v).to_excel(os.path.join(outdir_s2_ipa, "full_dmr_genes_for_ipa_%s.xlsx" % k.lower()))

    # Compute cross-comparison using cluster IDs (only)

    dmr_res_s2_idx = dict([
        (k, v.keys()) for k, v in dmr_res_sign_s2.items()
    ])

    cc_dict_dmr = cross_comparison.compute_cross_comparison_correction(
        dmr_res_s2_idx,
        pids,
        external_refs=external_refs_dm_labels,
        set_type='pair_only'
    )

    po_specific_to_all_refs_dmr = sorted(cc_dict_dmr['specific_to_all_refs'])
    pair_only_dmr = cc_dict_dmr['venn_set']

    # export list of reference-specific pair only cluster IDs, together with DMR information
    if len(po_specific_to_all_refs_dmr) > 0:
        this_export = {}
        for pid in pids:
            this_tbl = expand_dmr_results_table_by_gene_and_cluster(
                dmr_res_s2[pid][pid].to_table(include='all', skip_geneless=False),
                drop_geneless=False
            )
            this_tbl = this_tbl.loc[this_tbl.cluster_id.isin(po_specific_to_all_refs_dmr)]
            this_export[pid] = this_tbl

        excel.pandas_to_excel(this_export, os.path.join(outdir_s2, "consistently_in_pair_only_across_all_refs_dmr.xlsx"))

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

        # expand by gene
        this_sign = expand_dmr_results_table_by_gene_and_cluster(this_sign, drop_geneless=False)

        this_full = dmr_res_s2[pid][pid].to_table(include='all', skip_geneless=False)

        # annotate pair only
        po_col = pd.Series('N', index=this_full.index)
        po_col.loc[this_dmr] = 'Y'
        this_full.insert(this_full.shape[1], 'pair_only', po_col)

        # annotate significant (otherwise it isn't obvious)
        dm_col = pd.Series('N', index=this_full.index)
        dm_col.loc[sign_cl_ids] = 'Y'
        this_full.insert(this_full.shape[1], 'significant_dmr', dm_col)

        # expand by gene
        this_full = expand_dmr_results_table_by_gene_and_cluster(this_full, drop_geneless=False)

        po_dmr_export_full[pid] = this_full
        po_dmr_export[pid] = this_sign

    excel.pandas_to_excel(po_dmr_export, os.path.join(outdir_s2, 'pair_only_dmr.xlsx'))
    excel.pandas_to_excel(po_dmr_export_full, os.path.join(outdir_s2, 'pair_only_dmr_full.xlsx'))

    # export with a different layout, analogous to strategy 1

    # get all gene relation choices
    gene_relation_choices = set()
    for pp in dmr_res_s2.clusters.values():
        gene_relation_choices.update([t[1] for t in list(pp.genes)])

    venn_set, venn_ct = setops.venn_from_arrays(*[po_dmr_export[pid].index for pid in pids])
    po_combination_export = setops.venn_set_to_wide_dataframe(
        po_dmr_export,
        venn_set,
        pids,
        cols_to_include=('median_delta', 'padj'),
        static_cols_to_include=['cluster_id', 'gene'] + sorted(gene_relation_choices),
        consistency_check_col='median_delta',
        consistency_check_method="sign"
    )

    po_combination_export.to_excel(os.path.join(outdir_s2, 'pair_only_dmr_wideform.xlsx'))

    # array of Venn plots: GBM vs X (# DMRs)
    # only possible if number of references <= 3 (otherwise too many sets)
    if len(external_refs_dm_labels) < 4:
        nrows = len(subgroups)
        ncols = max([len(t) for t in subgroups.values()])
        set_labels = ['Syngeneic iNSC'] + external_refs_dm_labels
        set_colours = ['r', 'g', 'b']
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
        # assume we'll have one axis leftover!
        # TODO: if not, we'll need to add one or move the legend outside the main figure?
        cax = axs[-1, -1]
        ax_set = set(axs.flat)
        for pid in pids:
            sg = subgroups_lookup[pid]
            sg_members = subgroups[sg]
            i = subgroups.keys().index(sg)
            j = sg_members.index(pid)
            the_lists = [
                dmr_res_s2_idx[(pid, r)] for r in [pid] + external_refs_dm_labels
            ]
            venn_sets, cts = setops.venn_from_arrays(*the_lists)
            venn.venn_diagram(*the_lists, set_labels=None, ax=axs[i, j])
            axs[i, j].set_title("GBM%s" % pid, y=0.95)
            ax_set.remove(axs[i, j])
            # resize text
            h_txt = [t for t in axs[i, j].get_children() if isinstance(t, text.Text)]
            plt.setp(h_txt, fontsize=12)
        for ax in ax_set:
            ax.axis('off')

        # legend
        leg_objs = [
            patches.Rectangle([0, 0], 1, 1, color=c, alpha=0.4, label=l) for c, l in zip(set_colours, set_labels)
        ]
        cax.legend(handles=leg_objs, loc='center')

        fig.subplots_adjust(left=0., right=1., bottom=0., top=.93, wspace=0., hspace=0.05)

        fig.savefig(os.path.join(outdir_s2, 'venn_number_dmr_all_nsc.png'), dpi=200)
        fig.savefig(os.path.join(outdir_s2, 'venn_number_dmr_all_nsc.tiff'), dpi=200)

    dmr_members = {}
    for pid in pids:
        for r in external_refs_dm_labels:
            this_tbl = dmr_res_s2[pid][r].to_table(include='significant', skip_geneless=False)
            dmr_members[(pid, r)] = this_tbl.index

    # Venn plot showing the number of DMRs when each external ref is used as a comparator
    # This doesn't make any sense unless the number of external references is > 1
    if len(external_refs_dm_labels) > 1:
        nrows = len(subgroups)
        ncols = max([len(t) for t in subgroups.values()])
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
        ax_set = set(axs.flat)
        for pid in pids:
            sg = subgroups_lookup[pid]
            sg_members = subgroups[sg]
            i = subgroups.keys().index(sg)
            j = sg_members.index(pid)
            the_lists = [
                dmr_members[(pid, r)] for r in external_refs_dm_labels
            ]
            venn_sets, cts = setops.venn_from_arrays(*the_lists)
            venn.venn2(cts, set_labels=external_refs_dm_labels, ax=axs[i, j])
            axs[i, j].set_title("GBM%s vs..." % pid)
            ax_set.remove(axs[i, j])

        fig.subplots_adjust(left=0.03, right=0.98, bottom=.02, top=0.95, hspace=0.2, wspace=0.15)
        for ax in ax_set:
            ax.set_visible(False)
        fig.savefig(os.path.join(outdir_s2, 'number_dmr_multiple_references.png'), dpi=200)
        fig.savefig(os.path.join(outdir_s2, 'number_dmr_multiple_references.tiff'), dpi=200)

    # plot: how many DE genes are in the pair only comparison when each reference is used?
    # this only makes sense when the number of external references is > 1
    if len(external_refs_dm_labels) > 1:
        # at the same time, get numbers for a bar chart about % overlap
        n_pair_only_intersect = pd.DataFrame(0, index=pids, columns=external_refs_dm_labels)

        fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
        ax_set = set(axs.flat)
        for pid in pids:
            sg = subgroups_lookup[pid]
            sg_members = subgroups[sg]
            i = subgroups.keys().index(sg)
            j = sg_members.index(pid)
            the_lists = [
                set(pair_only_dmr.loc[pid, r]) for r in external_refs_dm_labels
            ]
            venn_sets, cts = setops.venn_from_arrays(*the_lists)
            venn.venn2(cts, set_labels=external_refs_dm_labels, ax=axs[i, j])
            axs[i, j].set_title("GBM%s pair only" % pid)
            ax_set.remove(axs[i, j])

            for i, r in enumerate(external_refs_dm_labels):
                n_pair_only_intersect.loc[pid, r] = cts[''.join(['1'] * len(external_refs_dm_labels))]

        fig.tight_layout()
        for ax in ax_set:
            ax.set_visible(False)
        fig.savefig(os.path.join(outdir_s2, 'number_po_dmr_multiple_references.png'), dpi=200)
        fig.savefig(os.path.join(outdir_s2, 'number_po_dmr_multiple_references.tiff'), dpi=200)

    # plot: overlap between individual references in terms of PO genes shared
    # this only makes sense when the number of external references is > 1
    if len(external_refs_dm_labels) > 1:
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

    nrows = len(subgroups)
    ncols = max([len(t) for t in subgroups.values()])
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
    ax_set = set(axs.flat)
    for pid in pids:
        sg = subgroups_lookup[pid]
        sg_members = subgroups[sg]
        i = subgroups.keys().index(sg)
        j = sg_members.index(pid)
        the_lists = [
            de_dmr_by_gene[(pid, r)] for r in external_refs_de_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_de_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s vs..." % pid)
        ax_set.remove(axs[i, j])
    fig.tight_layout()
    for ax in ax_set:
        ax.set_visible(False)
    fig.savefig(os.path.join(outdir_s2, 'number_de_dmr_multiple_references.png'), dpi=200)
    fig.savefig(os.path.join(outdir_s2, 'number_de_dmr_multiple_references.tiff'), dpi=200)

    # plot: how many DM-concordant DE genes are in the pair only comparison when each reference is used?

    # at the same time, get numbers for a bar chart about % overlap
    n_pair_only_intersect = pd.DataFrame(0, index=pids, columns=external_refs_de_labels)

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
    ax_set = set(axs.flat)
    for pid in pids:
        sg = subgroups_lookup[pid]
        sg_members = subgroups[sg]
        i = subgroups.keys().index(sg)
        j = sg_members.index(pid)
        the_lists = [
            set(pair_only_de_dmr.loc[pid, r]) for r in external_refs_de_labels
        ]
        venn_sets, cts = setops.venn_from_arrays(*the_lists)
        venn.venn2(cts, set_labels=external_refs_de_labels, ax=axs[i, j])
        axs[i, j].set_title("GBM%s pair only" % pid)
        ax_set.remove(axs[i, j])

        for i, r in enumerate(external_refs_de_labels):
            n_pair_only_intersect.loc[pid, r] = cts[''.join(['1'] * len(external_refs_de_labels))]

    fig.tight_layout()
    for ax in ax_set:
        ax.set_visible(False)
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
