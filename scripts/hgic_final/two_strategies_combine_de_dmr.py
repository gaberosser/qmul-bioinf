import multiprocessing as mp
import os
import pandas as pd
import numpy as np
import pickle
import re
from scipy import stats
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


def annotate_de_dmr_wide_form(data):
    """
    Annotate (in-place) wideform DataFrame containing combined DE/DMR results.
    :param data:
    :return:
    """
    data.insert(0, 'dmr_cluster_id', [t[0] for t in data.index])
    data.insert(0, 'gene', [t[1] for t in data.index])



def fisher_test_concordance(de_dmr_dat, de_col='de_logFC', dmr_col='dmr_median_delta'):
    x = de_dmr_dat[de_col]
    y = de_dmr_dat[dmr_col]
    cont = np.array([
        [((x < 0) & (y < 0)).sum(), ((x < 0) & (y > 0)).sum()],
        [((x > 0) & (y < 0)).sum(), ((x > 0) & (y > 0)).sum()],
    ])
    prior_odds, pval = stats.fisher_exact(cont)
    return cont, prior_odds, pval


def export_to_ipa(
        data_wideform,
        pids,
        effect_size_col='{pid}_de_logFC',
        pval_col='{pid}_de_FDR'
):
    cols = []
    for p in pids:
        cols.append("%s logFC" % p)
        cols.append("%s FDR" % p)
    ipa_export = pd.DataFrame(index=data_wideform['gene'].dropna().unique(), columns=cols)
    for pid in pids:
        ix = data_wideform.loc[:, pid] == 'Y'
        g_ix = data_wideform.loc[ix, 'gene'].dropna().unique()
        # grab effect size and pval columns
        z = data_wideform.loc[ix, effect_size_col.format(pid=pid)].astype(float)
        p = data_wideform.loc[ix, pval_col.format(pid=pid)].astype(float)
        # reduce to unique gene list by a mean groupby() aggregation
        # in practice, values for identical genes are also identical, so the mean operation is purely for convenience
        # check this first, just to be sure
        if not (z.groupby(data_wideform.loc[ix, 'gene']).apply(lambda x: len(set(x))) == 1).all():
            raise ValueError("Expected effect size for DE/DMR results to be identical for all matching genes, but this was contravened.")
        if not (p.groupby(data_wideform.loc[ix, 'gene']).apply(lambda x: len(set(x))) == 1).all():
            raise ValueError("Expected pval for DE/DMR results to be identical for all matching genes, but this was contravened.")
        ipa_export.loc[g_ix, "%s logFC" % pid] = z.groupby(data_wideform.loc[ix, 'gene']).mean()
        ipa_export.loc[g_ix, "%s FDR" % pid] = p.groupby(data_wideform.loc[ix, 'gene']).mean()
    return ipa_export


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
        external_ref_names_dm = ['GSE38216']
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

    external_refs_de_dm_labels = list(setops.reduce_intersection(external_refs_dm_labels, external_refs_de_labels))


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
    me_meta.insert(0, 'patient_id', me_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

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

    de_dmr_by_member = [joint_de_dmr_s1[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*de_dmr_by_member)

    # we don't add the null set here - it takes too long

    # generate wide-form lists and save to Excel file
    # first, export everything (concordant and discordant)
    dyn_col_include = ['de_logFC', 'de_FDR', 'de_direction', 'dmr_median_delta', 'dmr_padj', 'dmr_direction']
    static_col_include = [t for t in joint_de_dmr_s1[pids[0]].columns if re.match('(de|dmr)_', t) and t not in dyn_col_include]

    data_s1 = setops.venn_set_to_wide_dataframe(
        joint_de_dmr_s1,
        venn_set,
        pids,
        cols_to_include=dyn_col_include,
        static_cols_to_include=static_col_include,
        consistency_check_col='de_logFC',
        consistency_check_method='sign'
    )
    # add gene symbol and cluster ID back in
    annotate_de_dmr_wide_form(data_s1)
    data_s1.to_excel(os.path.join(outdir_s1, 'de_dmr_significant.xlsx'), index=False)

    # now filter by concordance and recompute venn sets
    joint_de_dmr_concordant_s1 = {}
    for pid in pids:
        a = joint_de_dmr_s1[pid].de_direction
        b = joint_de_dmr_s1[pid].dmr_direction
        idx = ((a == 'U') & (b == 'Hypo')) | ((a == 'D') & (b == 'Hyper'))
        joint_de_dmr_concordant_s1[pid] = joint_de_dmr_s1[pid].loc[idx]
        print "Patient %s has %d combined DE / DMR results, of which %d are concordant (%.2f %%)" % (
            pid, idx.size, idx.sum(), 100. * idx.sum() / float(idx.size)
        )
    de_dmr_by_member_concordant = [joint_de_dmr_concordant_s1[pid].index for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*de_dmr_by_member_concordant)

    # export only concordant results
    data_concordant_s1 = setops.venn_set_to_wide_dataframe(
        joint_de_dmr_concordant_s1,
        venn_set,
        pids,
        cols_to_include=dyn_col_include,
        static_cols_to_include=static_col_include,
        consistency_check_col='de_logFC',
        consistency_check_method='sign'
    )
    # add gene symbol and cluster ID back in
    annotate_de_dmr_wide_form(data_concordant_s1)
    data_concordant_s1.to_excel(os.path.join(outdir_s1, 'de_dmr_significant_concordant.xlsx'), index=False)

    # expanded core
    # data_ec = setops.venn_set_to_wide_dataframe(
    #     joint_de_dmr_concordant_s1,
    #     venn_set,
    #     pids,
    #     include_sets=ec_sets,
    #     cols_to_include=col_include,
    #     consistency_check_col='de_logFC',
    #     consistency_check_method='sign'
    # )
    # annotate_de_dmr_wide_form(data_ec)
    # data_ec.to_excel(os.path.join(outdir_s1, 'expanded_core_de_dmr_concordant.xlsx'), index=False)

    # subgroup-specific
    # data_ss = setops.venn_set_to_wide_dataframe(
    #     joint_de_dmr_concordant_s1,
    #     venn_set,
    #     pids,
    #     include_sets=ss_sets,
    #     cols_to_include=col_include,
    #     consistency_check_col='de_logFC',
    #     consistency_check_method='sign'
    # )
    # annotate_de_dmr_wide_form(data_ss)
    # data_ss.to_excel(os.path.join(outdir_s1, 'subgroup_specific_de_dmr_concordant.xlsx'), index=False)

    # patient unique
    # data_pu = setops.venn_set_to_wide_dataframe(
    #     joint_de_dmr_concordant_s1,
    #     venn_set,
    #     pids,
    #     include_sets=pu_sets,
    #     cols_to_include=col_include,
    #     consistency_check_col='de_logFC',
    #     consistency_check_method='sign'
    # )
    # annotate_de_dmr_wide_form(data_pu)
    # data_pu.to_excel(os.path.join(outdir_s1, 'patient_unique_de_dmr_concordant.xlsx'), index=False)

    # patient and subgroup-specific
    # data_pss = {}
    # for pid in pids:
    #     this_tbl =  setops.venn_set_to_wide_dataframe(
    #         joint_de_dmr_concordant_s1,
    #         venn_set,
    #         pids,
    #         include_sets=pss_sets[pid],
    #         cols_to_include=col_include,
    #         consistency_check_col='de_logFC',
    #         consistency_check_method='sign'
    #     )
    #     annotate_de_dmr_wide_form(this_tbl)
    #     data_pss[pid] = this_tbl
    #
    # excel.pandas_to_excel(data_pss, os.path.join(outdir_s1, 'patient_and_subgroup_specific_de_dmr_concordant.xlsx'))

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

    upset['figure'].savefig(os.path.join(outdir_s1, "upset_de_dmr_concordant.png"), dpi=200)
    upset['figure'].savefig(os.path.join(outdir_s1, "upset_de_dmr_concordant.tiff"), dpi=200)

    # export to IPA for pathway analysis
    export_to_ipa(data_s1, pids).to_excel(os.path.join(outdir_s1_ipa, "de_dmr_significant_for_ipa.xlsx"))
    export_to_ipa(data_concordant_s1, pids).to_excel(os.path.join(outdir_s1_ipa, "concordant_de_dmr_significant_for_ipa.xlsx"))

    # out of interest, how does this compare to the alternative approach where we require patient-specific DE and DMR
    # separately, then look for concordance/matching?

    # start with concordant results (data_s1_concordant)
    venn_set, venn_ct = setops.venn_from_arrays(*de_dmr_by_member_concordant)

    # repeat with the other approach
    vs_dm, vc_dm = setops.venn_from_arrays(*[dmr_res_s1[pid].results_significant.keys() for pid in pids])
    vs_de, vc_de = setops.venn_from_arrays(*[de_res_s1[pid]['Gene Symbol'].dropna() for pid in pids])
    ps_dm = dict([
        (pid, vs_dm[v]) for pid, v in setops.specific_sets(pids).items()
    ])
    ps_dm_genes = {}
    for pid, cids in ps_dm.items():
        ps_dm_genes[pid] = sorted(setops.reduce_union(*[[t[0] for t in dmr_res_s1.clusters[c].genes] for c in cids]))
    ps_de = dict([
        (pid, vs_de[v]) for pid, v in setops.specific_sets(pids).items()
    ])

    ps_de_dm = {}
    for pid in pids:
        this_genes = set(ps_de[pid]).intersection(ps_dm_genes[pid])
        this_clusters = []
        for c in ps_dm[pid]:
            this = [t[0] for t in dmr_res_s1.clusters[c].genes]
            if len(this_genes.intersection(this)):
                this_clusters.extend([(c, g) for g in this_genes.intersection(this)])
        ps_de_dm[pid] = joint_de_dmr_concordant_s1[pid].reindex(this_clusters).dropna(axis=0, how='all')

    from plotting.genomics import MethylationExpressionLocusPlotter
    import collections

    dmr_comparison_groups = collections.OrderedDict([(pid, {}) for pid in consts.PIDS])
    gg = me_data.columns.groupby(zip(me_meta.patient_id, me_meta.type))
    for (pid, typ), samples in gg.items():
        dmr_comparison_groups[pid][typ] = samples

    plt_obj = MethylationExpressionLocusPlotter()
    plt_obj.set_mvalues(me_data)
    plt_obj.set_dmr_res(dmr_res_s1, dmr_comparison_groups)
    plt_obj.set_de_res(de_res_s1)

    # shortlist
    shortlist = sorted(setops.reduce_union(*[[u[1] for u in t.index] for t in ps_de_dm.values()]))

    import datetime
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    # Create the PdfPages object to which we will save the pages:
    # The with statement makes sure that the PdfPages object is closed properly at
    # the end of the block, even if an Exception occurs.
    with PdfPages(os.path.join(outdir, "de_dmr_shortlist_mex_plots.pdf")) as pdf:
        for g in shortlist:
            the_obj = plt_obj.plot_gene(g)
            the_obj.dm_axs[pids[0]].set_title(g)
            the_obj.gs.update(top=0.95)
            pdf.savefig(the_obj.fig)
            plt.close(the_obj.fig)

        # pdf.attach_note("plot of sin(x)")  # you can add a pdf note to
        # attach metadata to a page

        # We can also set the file's metadata via the PdfPages object:
        d = pdf.infodict()
        d['Title'] = 'Patient-specific DE/DMR, shortlist'
        d['Author'] = 'Gabriel Rosser'
        d['Subject'] = ''
        d['CreationDate'] = datetime.datetime.now()
        d['ModDate'] = datetime.datetime.today()


    raise StopIteration

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
        raise AttributeError("Unable to load pre-computed DE results, expected at %s" % fn)

    de_res_s2 = {}
    for pid in pids:
        t = de_res_full_s2[(pid, pid)]
        de_res_s2.setdefault('syngeneic', {})[pid] = t.loc[t.FDR < de_params['fdr']]
        for k in external_refs_de_labels:
            t = de_res_full_s2[(pid, k)]
            de_res_s2.setdefault(k, {})[pid] = t.loc[t.FDR < de_params['fdr']]

    # Load DMR cross-comparison

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s2

    the_hash = dmr_results_hash(me_obj_with_ref.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_cross_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S2 DMR results from %s", fn)
        dmr_res_s2 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise AttributeError("Unable to load pre-computed DMR results, expected at %s" % fn)

    dmr_res_full_s2 = {}
    dmr_res_sign_s2 = {}

    for k1 in dmr_res_s2.keys():
        dmr_res_full_s2.setdefault('syngeneic', {})[k1] = dmr_res_s2[k1][k1]
        dmr_res_sign_s2.setdefault('syngeneic', {})[k1] = dmr_res_s2[k1][k1]
        for k2 in external_refs_dm_labels:
            dmr_res_full_s2.setdefault(k2, {})[k1] = dmr_res_s2[k1][k2]
            dmr_res_sign_s2.setdefault(k2, {})[k1] = dmr_res_s2[k1][k2]

    # get the joint table (full)
    joint_de_dmr_s2 = {}
    for k in ['syngeneic'] + external_refs_de_dm_labels:
        this_joint = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_sign_s2[k], de_res_s2[k])
        for k2, v in this_joint.items():
            joint_de_dmr_s2["%s_%s" % (k2, k)] = v

    # generate wide-form lists and save to Excel file
    # first, export everything (concordant and discordant)

    # we can't use the Venn approach here, but we don't need to
    ix = sorted(setops.reduce_union(*[t.index for t in joint_de_dmr_s2.values()]))
    blocks_for_export = []
    static_block = pd.DataFrame(index=ix, columns=static_col_include)

    for k, v in joint_de_dmr_s2.items():
        this_yn = pd.Series('N', index=ix)
        this_yn.loc[v.index] = 'Y'
        this_block = v.reindex(ix)[dyn_col_include]
        this_block.columns = ["%s_%s" % (k, t) for t in this_block.columns]
        this_block.insert(0, k, this_yn)
        blocks_for_export.append(this_block)

        # fill in static attributes where found
        static_block.loc[v.index] = v[static_col_include]

    blocks_for_export = [static_block] + blocks_for_export

    data_s2 = pd.concat(blocks_for_export, axis=1)
    annotate_de_dmr_wide_form(data_s2)

    data_s2.to_excel(os.path.join(outdir_s2, "de_dmr_significant.xlsx"))

    # get the joint table (concordant only)

    joint_de_dmr_concordant_s2 = {}
    for k, v in joint_de_dmr_s2.items():
        a = v.de_direction
        b = v.dmr_direction
        idx = ((a == 'U') & (b == 'Hypo')) | ((a == 'D') & (b == 'Hyper'))
        joint_de_dmr_concordant_s2[k] = v.loc[idx]
        print "Comparison %s has %d combined DE / DMR results, of which %d are concordant (%.2f %%)" % (
            k, idx.size, idx.sum(), 100. * idx.sum() / float(idx.size)
        )

    # generate wide-form lists and save to Excel file
    # now, export only concordant results
    # we can't use the Venn approach here, but we don't need to
    ix = sorted(setops.reduce_union(*[t.index for t in joint_de_dmr_concordant_s2.values()]))
    blocks_for_export = []
    static_block = pd.DataFrame(index=ix, columns=static_col_include)

    for k, v in joint_de_dmr_concordant_s2.items():
        this_yn = pd.Series('N', index=ix)
        this_yn.loc[v.index] = 'Y'
        this_block = v.reindex(ix)[dyn_col_include]
        this_block.columns = ["%s_%s" % (k, t) for t in this_block.columns]
        this_block.insert(0, k, this_yn)
        blocks_for_export.append(this_block)

        # fill in static attributes where found
        static_block.loc[v.index] = v[static_col_include]

    blocks_for_export = [static_block] + blocks_for_export
    data_concordant_s2 = pd.concat(blocks_for_export, axis=1)
    annotate_de_dmr_wide_form(data_concordant_s2)

    data_concordant_s2.to_excel(os.path.join(outdir_s2, "de_dmr_significant_concordant.xlsx"))

    # export to IPA for pathway analysis
    # do this in separate blocks (by comparator) to fit in batch upload size restriction
    for k in ['syngeneic'] + external_refs_de_dm_labels:
        cols = data_s2.columns[data_s2.columns.str.contains(k)
        ]
        for_export = data_s2[['gene', 'dmr_cluster_id', 'dmr_chr'] + cols.tolist()]
        for_export.columns = [t.replace("_%s" % k, "") for t in for_export.columns]
        export_to_ipa(for_export, pids).to_excel(os.path.join(outdir_s2_ipa, "de_dmr_significant_for_ipa_%s.xlsx" % k))

        for_export = data_concordant_s2[['gene', 'dmr_cluster_id', 'dmr_chr'] + cols.tolist()]
        for_export.columns = [t.replace("_%s" % k, "") for t in for_export.columns]
        export_to_ipa(for_export, pids).to_excel(os.path.join(outdir_s2_ipa, "de_dmr_significant_concordant_for_ipa_%s.xlsx" % k))

    # Venn diagram array to quantify numbers
    # for this purpose, reduce to the number of genes (as one gene can appear more than once)
    # array of Venn plots: GBM vs X (# DE+DM genes)
    # only possible if number of references <= 3 (otherwise too many sets)
    if len(external_refs_de_dm_labels) < 4:
        nrows = len(subgroups)
        ncols = max([len(t) for t in subgroups.values()])
        set_labels = ['Syngeneic iNSC'] + external_refs_de_dm_labels
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
                data_concordant_s2['gene'].loc[data_concordant_s2["%s_%s" % (pid, k)] == 'Y']
                for k in ['syngeneic'] + external_refs_de_dm_labels
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

        fig.savefig(os.path.join(outdir_s2, 'venn_number_de_dmr_all_nsc.png'), dpi=200)
        fig.savefig(os.path.join(outdir_s2, 'venn_number_de_dmr_all_nsc.tiff'), dpi=200)
