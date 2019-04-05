from plotting import bar, common, pie, polar, venn
from methylation import loader, dmr, process, annotation_gene_to_ensembl
import pandas as pd
from stats import nht
from utils import output, setops, genomics, log, ipa, dictionary
from cytoscape import cyto
import multiprocessing as mp
import os
import collections
import operator
import pickle
import numpy as np
from scipy import stats, cluster
import matplotlib
from matplotlib import pyplot as plt, patches, gridspec
from matplotlib.colors import Normalize
from matplotlib import cm
from sklearn.neighbors import KernelDensity
import seaborn as sns

import references
from scripts.hgic_final import \
    two_strategies_grouped_dispersion as tsgd, \
    two_strategies_combine_de_dmr as tscd, \
    analyse_dmrs_s1_direction_distribution as addd, \
    consts
from scripts.dmr_direction_bias_story import \
    same_process_applied_to_de as same_de
from rnaseq import loader as rnaseq_loader, filter
from integrator import rnaseq_methylationarray
from scripts.methylation import dmr_values_to_bigwig

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
logger = log.get_console_logger()

"""
Here we seek to identify DEs corresponding to DMRs that distinguish the hypo- and hypermethylated groups (known as 
group-specific DMRs).

Created this new script because the previous one (originally on DMRs only) was getting bloated. 
Restrict ourselves to DE/DMRs here.

"""

def get_genes_relations(cluster_ids, clusters, relation_map=None, relation_priority=None, relation_filter=None):
    if relation_priority is None:
        relation_priority = [
            'TSS200',
            'TSS1500',
            '1stExon',
            "3'UTR",
            "5'UTR",
            'ExonBnd',
            'Body'
        ]
    if relation_map is None:
        relation_map = dict([(k, k) for k in relation_priority])

    if relation_filter is not None:
        if not hasattr(relation_filter, '__iter__'):
            relation_filter = [relation_filter]
        relation_filter = set(relation_filter)

    rels = []
    genes = set()

    for t in cluster_ids:
        pc = clusters[t]
        if len(pc.genes) > 0:
            gs, rs = zip(*clusters[t].genes)
            if relation_filter is not None:
                # filter gene and relation
                # if this leaves no results, skip this DMR
                tmp = [(g, r) for g, r in zip(gs, rs) if r in relation_filter]
                if len(tmp) == 0:
                    continue
                else:
                    gs, rs = zip(*tmp)
            for rp in relation_priority:
                if rp in rs:
                    rels.append(relation_map[rp])
                    break
            genes.update(gs)
        else:
            # if we aren't filtering by relation, add the (arbitrary) label 'intergene' here to express a lack of
            # gene.
            if relation_filter is None:
                rels.append('Intergene')
    return sorted(genes), rels


def pct_concordant(joint_res):
    a = (np.sign(joint_res.de_logFC) != np.sign(joint_res.dmr_median_delta)).sum()
    b = float(joint_res.shape[0])
    return a / b * 100


def tabulate_de_counts_by_direction(de_res, pids=consts.PIDS, **gene_lists):
    """
    Given the lists of genes associated with the hypo and hyper group, tabulate the numbers of matching genes in
    the supplied DE results for each patient. Split counts by DE logFC direction.
    :param de_res:
    :param genes_hypo:
    :param genes_hyper:
    :param pids:
    :return: Raw counts table, Table expressing coutns as a % of the total in that direction
    """
    cols = reduce(
        lambda x, y: x + y,
        [["%s up" % k, "%s down" % k] for k in gene_lists]
    )

    # table of DE counts (in DMR-linked context)
    de_count_table = pd.DataFrame(
        0,
        index=pids,
        columns=cols + ['Total up', 'Total down']
    )

    for pid in pids:
        de_count_table.loc[pid, 'Total up'] = (de_res[pid]['logFC'] > 0).sum()
        de_count_table.loc[pid, 'Total down'] = (de_res[pid]['logFC'] < 0).sum()
        for k, g_arr in gene_lists.items():
            ix = de_res[pid].index.intersection(references.gene_symbol_to_ensembl(g_arr).dropna().values)
            de_count_table.loc[pid, '%s up' % k] = (de_res[pid].loc[ix, 'logFC'] > 0).sum()
            de_count_table.loc[pid, '%s down' % k] = (de_res[pid].loc[ix, 'logFC'] < 0).sum()

    # express this as a pct of the total up/down
    de_count_table_pct = pd.DataFrame(
        index=pids,
        columns=cols
    )
    for k in gene_lists:
        for k2 in ['up', 'down']:
            de_count_table_pct.loc[pids, "%s %s" % (k, k2)] = \
                de_count_table.loc[pids, "%s %s" % (k, k2)] / de_count_table.loc[pids, 'Total %s' % k2].astype(float) * 100.

    return de_count_table, de_count_table_pct


def get_de_dmr_groups(
        joint_de_dmr,
        clusters,
        groups,
        pids=consts.PIDS,
        relation_filter=None
):
    """
    Get group-specific DE/DMRs. These are defined as DEs that are consistent with the DMRs in a given selection of
    patients (from one to many) that are NOT shared across groups.
    :param joint_de_dmr:
    :param clusters:
    :param groups: Dictionary, keyed by group name. Values are iterables giving patient IDs in each group.
    :param pids:
    :param relation_filter:
    :return:
    """
    venn_sets_by_group = setops.full_partial_unique_other_sets_from_groups(pids, groups)

    if relation_filter is not None:
        if not hasattr(relation_filter, '__iter__'):
            relation_filter = [relation_filter]

    de_dmr_groups = {}
    de_dmr_de_logfc = {}
    de_dmr_de_fdr = {}
    de_dmr_dmr_delta = {}

    if relation_filter is None:
        de_dmr_by_member = [joint_de_dmr[pid].index for pid in pids]
    else:
        de_dmr_by_member = []
        for pid in pids:
            this_members = []
            for t in joint_de_dmr[pid].index:
                gene_rel_options = [(t[1], rel) for rel in relation_filter]
                if len(set(clusters[t[0]].genes).intersection(gene_rel_options)) > 0:
                    this_members.append(t)
            de_dmr_by_member.append(this_members)
    venn_set, venn_count = setops.venn_from_arrays(*de_dmr_by_member)

    for grp in groups:
        this_sets = venn_sets_by_group['full'][grp] + venn_sets_by_group['partial'][grp]
        this_de_dmrs = sorted(setops.reduce_union(*[venn_set[k] for k in this_sets]))

        if relation_filter is not None:
            new_de_dmrs = []
            for t in this_de_dmrs:
                # look for any intersection here
                gene_rel_options = [(t[1], rel) for rel in relation_filter]
                if len(set(clusters[t[0]].genes).intersection(gene_rel_options)) > 0:
                    new_de_dmrs.append(t)
            this_de_dmrs = new_de_dmrs

        de_dmr_groups[grp] = this_de_dmrs

        # get separate lists of DE genes and DMR IDs
        # DMRs is straightforward
        de_dmr_dmr_delta[grp] = pd.DataFrame(
            index=sorted(set([t[0] for t in this_de_dmrs])),
            columns=pids + ['consistent'],
        )
        # DEs is trickier: some genes have mapped twice because I was so diligent in curating the original lists!
        this_de_genes = sorted(set([t[1] for t in this_de_dmrs]))
        this_de_ens = annotation_gene_to_ensembl.gene_to_ens(this_de_genes)
        this_de_ens = this_de_ens[~this_de_ens.duplicated()]
        this_de_genes = this_de_ens.index

        de_dmr_de_logfc[grp] = pd.DataFrame(
            index=this_de_genes.tolist(),
            columns=pids + ['consistent'],
        )
        de_dmr_de_fdr[grp] = pd.DataFrame(
            index=this_de_genes.tolist(),
            columns=pids + ['consistent'],
        )

        # fill them in
        for k in this_sets:
            this_vs = [t for t in venn_set[k] if t[1] in this_de_genes]
            this_pids = [pids[i] for i, t in enumerate(k) if t == '1']
            for pid in this_pids:
                de_dmr_dmr_delta[grp].loc[[t[0] for t in this_vs], pid] = joint_de_dmr[pid].loc[
                    this_vs, 'dmr_median_delta'].values
                de_dmr_de_logfc[grp].loc[[t[1] for t in this_vs], pid] = joint_de_dmr[pid].loc[
                    this_vs, 'de_logFC'].values
                de_dmr_de_fdr[grp].loc[[t[1] for t in this_vs], pid] = joint_de_dmr[pid].loc[
                    this_vs, 'de_FDR'].values

        for k, row in de_dmr_dmr_delta[grp].iterrows():
            tmp_dm = np.sign(row.dropna().astype(float))
            row['consistent'] = (tmp_dm == tmp_dm.iloc[0]).all()

        for k, row in de_dmr_de_logfc[grp].iterrows():
            tmp_de = np.sign(row.dropna().astype(float))
            row['consistent'] = (tmp_de == tmp_de.iloc[0]).all()
            de_dmr_de_fdr[grp].loc[k, 'consistent'] = row['consistent']

    return {
        'dmr_median_delta_m': de_dmr_dmr_delta,
        'de_logFC': de_dmr_de_logfc,
        'de_FDR': de_dmr_de_fdr,
        'de_dmr_groups': de_dmr_groups
    }


def export_de_dmr_groups_for_ipa(de_fdr, de_logfc, groups, fn_out=None, pids=consts.PIDS):
    """

    :param de_fdr: Output of `get_de_dmr_groups`
    :param de_logfc: Output of `get_de_dmr_groups`
    :param groups:
    :param fn_out: If supplied, the IPA results (in Excel format) will be written to this path
    :param pids:
    :return:
    """
    # export these for IPA analysis
    df_for_ipa = pd.DataFrame(
        index=sorted(setops.reduce_union(*[t.index for t in de_logfc.values()])),
        columns=reduce(operator.add, [["%s_logFC" % pid, "%s_FDR" % pid] for pid in pids])
    )
    for grp in groups:
        for pid in groups[grp]:
            this_logfc = de_logfc[grp][pid].dropna()
            this_fdr = de_fdr[grp][pid].dropna()
            df_for_ipa.loc[this_logfc.index, "%s_logFC" % pid] = this_logfc
            df_for_ipa.loc[this_fdr.index, "%s_FDR" % pid] = this_fdr
    if fn_out is not None:
        df_for_ipa.to_excel(fn_out)
    return df_for_ipa


def plot_venn_de_directions(
    logfc,
    set_colours_dict,
    ax=None,
    set_labels=('Hypo', 'Hyper')
):
    if ax is None:
        fig = plt.figure(figsize=(5., 3.3))
        ax = fig.add_subplot(111)

    vv, vs, vc = venn.venn_diagram(
        *[logfc[k].index[logfc[k]['consistent'].astype(bool)] for k in set_labels],
        set_labels=set_labels,
        set_colors=[set_colours_dict[t] for t in set_labels],
        ax=ax
    )
    ax.figure.tight_layout()

    # modify labels based on direction
    this_members = collections.OrderedDict()
    for k in set_labels:
        this_res = logfc[k]
        this_ix = this_res['consistent'].astype(bool)
        this_res = np.sign(this_res.loc[this_ix].astype(float).mean(axis=1))
        this_members["%s up" % k] = this_res.index[this_res > 0].difference(vs['11'])
        this_members["%s down" % k] = this_res.index[this_res < 0].difference(vs['11'])
        # get the corresponding label
        lbl = vv.get_label_by_id(setops.specific_sets(set_labels)[k])
        lbl.set_text(
            lbl.get_text()
            + '\n'
            + r'$%d\uparrow$' % len(this_members["%s up" % k])
            + '\n'
            + r'$%d\downarrow$' % (len(this_members["%s down" % k])
        ))

    return ax


if __name__ == '__main__':
    # set a minimum pval for pathways to be used
    alpha = 0.01
    plogalpha = -np.log10(alpha)
    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.05
    plogalpha_relevant = -np.log10(alpha_relevant)

    relations_tss = ['TSS1500', 'TSS200']

    IPA_PATHWAY_DIR = os.path.join(
        HGIC_LOCAL_DIR,
        'current/dmr_direction_bias_story/ipa/dmr_genes/pathways'
    )

    outdir = output.unique_output_dir()
    de_res_fn = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq', 'full_de_syngeneic_only.xlsx')
    pids = consts.PIDS
    generate_bw = False

    DE_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'de')
    de_params = consts.DE_PARAMS

    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    # should we add 'chr' prefix to bigwig output?
    chr_prefix = True

    subgroups = consts.SUBGROUPS
    subgroups_ind = setops.groups_to_ind(pids, subgroups)

    # load gene expression values
    rna_obj = rnaseq_loader.load_by_patient(pids, include_control=False, source='salmon')
    rna_obj.filter_samples(rna_obj.meta.index.isin(consts.S1_RNASEQ_SAMPLES))
    rna_tpm = rna_obj.data
    rna_meta = rna_obj.meta
    rna_meta.insert(0, 'patient_id', rna_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    # load DE results
    the_hash = tscd.de_results_hash(rna_obj.meta.index.tolist(), de_params)
    filename = 'de_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DE results from %s", fn)
        with open(fn, 'rb') as f:
            de_res_full_s1 = pickle.load(f)
    else:
        raise AttributeError("Unable to load pre-computed DE results, expected at %s" % fn)

    de_res_s1 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s1.items()])

    # load methylation data
    me_obj, anno = tsgd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
    me_data = me_obj.data
    me_meta = me_obj.meta
    me_meta.insert(0, 'patient_id', me_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    # load DMR results
    the_hash = tsgd.dmr_results_hash(me_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        logger.info("Unable to locate pre-existing results. Computing from scratch (this can take a while).")
        dmr_res_s1 = tsgd.paired_dmr(me_data, me_meta, anno, pids, dmr_params)
        # Save DMR results to disk
        dmr_res_s1.to_pickle(fn, include_annotation=False)
        logger.info("Saved DMR results to %s", fn)

    # extract full (all significant) results
    dmr_res_all = dmr_res_s1.results_significant

    # combine DE/DMR
    joint_de_dmr_s1 = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_s1)

    ## Identify group-specific DMRs

    groups = {
        'Hypo': ['019', '030', '031', '017'],
        'Hyper': ['018', '050', '054', '061', '026', '052']
    }
    group_ind = setops.groups_to_ind(pids, groups)
    groups_inv = dictionary.complement_dictionary_of_iterables(groups, squeeze=True)

    set_colours_dict = {
        'Hypo': 'g',
        'Hyper': 'r',
        'Discordant': 'b'
    }

    dmr_by_member = [dmr_res_all[pid].keys() for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*dmr_by_member)

    venn_sets_by_group = setops.full_partial_unique_other_sets_from_groups(pids, groups)
    dmr_groups = {}
    for grp in groups:
        # generate bar chart showing number / pct in each direction (DM)
        this_sets = venn_sets_by_group['full'][grp] + venn_sets_by_group['partial'][grp]
        this_dmrs = sorted(setops.reduce_union(*[venn_set[k] for k in this_sets]))
        dmr_groups[grp] = this_dmrs

    # Rather than just looking at genes corresponding to group-specific DMRs, we make the requirements more
    # stringent. For each Venn set (e.g. 018, 054, 052 - hyper group), we require DE genes in the same patients.
    # Simplest approach is to use the joint_de_dmr dataframes, which have already been combined.
    # all relations

    tmp = get_de_dmr_groups(joint_de_dmr_s1, dmr_res_s1.clusters, groups)
    de_dmrs_all = tmp['de_dmr_groups']
    de_dmr_de_fdr_all = tmp['de_FDR']
    de_dmr_de_logfc_all = tmp['de_logFC']
    de_dmr_ipa_res_all = export_de_dmr_groups_for_ipa(
        de_dmr_de_fdr_all,
        de_dmr_de_logfc_all,
        groups,
        os.path.join(outdir, "group_specific_de_for_ipa.xlsx")
    )

    # TSS relation only

    tmp = get_de_dmr_groups(joint_de_dmr_s1, dmr_res_s1.clusters, groups, relation_filter=['TSS200', 'TSS1500'])
    de_dmrs_tss = tmp['de_dmr_groups']
    de_dmr_de_fdr_tss = tmp['de_FDR']
    de_dmr_de_logfc_tss = tmp['de_logFC']
    de_dmr_ipa_res_tss = export_de_dmr_groups_for_ipa(
        de_dmr_de_fdr_tss,
        de_dmr_de_logfc_tss,
        groups,
        os.path.join(outdir, "group_specific_de_for_ipa_tss.xlsx")
    )

    # Venn diagrams
    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    plot_venn_de_directions(de_dmr_de_logfc_all, set_colours_dict, ax=ax)
    fig.savefig(os.path.join(outdir, "genes_from_group_spec_dmrs_all.png"), dpi=200)

    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    plot_venn_de_directions(de_dmr_de_logfc_tss, set_colours_dict, ax=ax)
    fig.savefig(os.path.join(outdir, "genes_from_group_spec_dmrs_tss.png"), dpi=200)

    # assess concordance between DM and DE direction
    # start with scatterplot

    # now run a null model and compare