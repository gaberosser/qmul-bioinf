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


if __name__ == "__main__":
    outdir = output.unique_output_dir(reuse_empty=True)
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
        raise AttributeError("To run this script, we require pre-computed DE results residing in %s" % DE_LOAD_DIR)

    if not os.path.isdir(DMR_LOAD_DIR):
        raise AttributeError("To run this script, we require pre-computed DMR results residing in %s" % DMR_LOAD_DIR)

    #########################################################################
    ### STRATEGY 1: No references, just compare GBM-iNSC for each patient ###
    #########################################################################

    # some quantities relating to set membership

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
        raise AttributeError("Unable to load pre-computed DMR results, expected at %s" % fn)

    # extract results
    dmr_res_full_s1 = dmr_res_s1.results
    dmr_res_sign_s1 = dmr_res_s1.results_significant

    rnaseq_obj = load_rnaseq(
        pids,
        external_ref_names_de,
        strandedness=external_ref_strandedness_de,
    )

    rnaseq_obj.filter_by_sample_name(consts.ALL_RNASEQ_SAMPLES)

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
        raise AttributeError("Unable to load pre-computed DE results, expected at %s" % fn)

    de_res_s1 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s1.items()])

    # get the joint table
    joint_de_dmr_s1 = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_s1)

    ## TODO: update from here

    # filter - only DMRs with TSS or island class
    # for pid in pids:
    #     idx = joint_de_dmr_s1[pid].dmr_class_island | joint_de_dmr_s1[pid].dmr_class_tss
    #     joint_de_dmr_s1[pid] = joint_de_dmr_s1[pid].loc[idx]
    #     joint_de_dmr_s1[pid].index = de_dmr_hash(joint_de_dmr_s1[pid])
    #     print "Patient %s. %d DE & DMR rows, %d are gene only, keeping %d." % (pid, idx.size, (~idx).sum(), idx.sum())

    de_dmr_by_member = [joint_de_dmr_s1[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*de_dmr_by_member)

    # we don't add the null set here - it takes too long

    # generate wide-form lists and save to Excel file
    # first, export everything (concordant and discordant)
    data = setops.venn_set_to_wide_dataframe(
        joint_de_dmr_s1,
        venn_set,
        pids,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        consistency_check_col='de_logfc',
        consistency_check_method='sign'
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
    data = setops.venn_set_to_wide_dataframe(
        joint_de_dmr_concordant_s1,
        venn_set,
        pids,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        consistency_check_col='de_logfc',
        consistency_check_method='sign'
    )
    # add gene symbol and cluster ID back in
    annotate_de_dmr_wide_form(data)
    data.to_excel(os.path.join(outdir_s1, 'full_de_dmr_concordant.xlsx'), index=False)

    # expanded core
    data_ec = setops.venn_set_to_wide_dataframe(
        joint_de_dmr_concordant_s1,
        venn_set,
        pids,
        include_sets=ec_sets,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        consistency_check_col='de_logfc',
        consistency_check_method='sign'
    )
    annotate_de_dmr_wide_form(data_ec)
    data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_de_dmr_concordant.xlsx'), index=False)

    # subgroup-specific
    data_ss = setops.venn_set_to_wide_dataframe(
        joint_de_dmr_concordant_s1,
        venn_set,
        pids,
        include_sets=ss_sets,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        consistency_check_col='de_logfc',
        consistency_check_method='sign'
    )
    annotate_de_dmr_wide_form(data_ss)
    data_ss.to_excel(os.path.join(outdir_s1, 'subgroup_specific_de_dmr_concordant.xlsx'), index=False)

    # patient unique
    data_pu = setops.venn_set_to_wide_dataframe(
        joint_de_dmr_concordant_s1,
        venn_set,
        pids,
        include_sets=pu_sets,
        cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
        consistency_check_col='de_logfc',
        consistency_check_method='sign'
    )
    annotate_de_dmr_wide_form(data_pu)
    data_pu.to_excel(os.path.join(outdir_s1, 'patient_unique_de_dmr_concordant.xlsx'), index=False)

    # patient and subgroup-specific
    data_pss = {}
    for pid in pids:
        this_tbl =  setops.venn_set_to_wide_dataframe(
            joint_de_dmr_concordant_s1,
            venn_set,
            pids,
            include_sets=pss_sets[pid],
            cols_to_include=('de_logfc', 'de_padj', 'dmr_median_delta', 'dmr_padj'),
            consistency_check_col='de_logfc',
            consistency_check_method='sign'
        )
        annotate_de_dmr_wide_form(this_tbl)
        data_pss[pid] = this_tbl

    excel.pandas_to_excel(data_pss, os.path.join(outdir_s1, 'patient_and_subgroup_specific_de_dmr_concordant.xlsx'))

    upset = venn.upset_plot_with_groups(
        de_dmr_by_member_concordant,
        pids,
        subgroup_ind,
        subgroup_set_colours,
        min_size=5,
        n_plot=30,
        default_colour='gray'
    )
    upset['axes']['main'].set_ylabel('Number of DM-concordant DE genes in set')
    upset['axes']['set_size'].set_xlabel('Number of DM-concordant DE genes\nin single comparison')
    upset['gs'].update(bottom=0.11)

    upset['figure'].savefig(os.path.join(outdir_s1, "upset_de_dmr.png"), dpi=200)
    upset['figure'].savefig(os.path.join(outdir_s1, "upset_de_dmr.tiff"), dpi=200)
