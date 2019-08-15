import collections
import multiprocessing as mp
import operator
import os
import pickle

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt, gridspec

from cytoscape import cyto
from integrator import rnaseq_methylationarray
from methylation import dmr, annotation_gene_to_ensembl
from plotting import bar, venn
from rnaseq import loader as rnaseq_loader
from scripts.dmr_direction_bias_story import \
    same_process_applied_to_de as same_de
from scripts.hgic_final import \
    two_strategies_grouped_dispersion as tsgd, \
    two_strategies_combine_de_dmr as tscd, \
    analyse_dmrs_s1_direction_distribution as addd, \
    consts
from scripts.methylation import dmr_values_to_bigwig
from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, INTERMEDIATE_DIR
from utils import output, setops, genomics, log, ipa, dictionary, reference_genomes

logger = log.get_console_logger()

"""
Here we seek to identify DMRs that distinguish the hypo- and hypermethylated groups.

This is initially carried out directly using the GIC-iNSC results.

- We query these DMRs to determine whether they have different location distributions.
- We look for shared but discordant DMRs between the two groups
- We run IPA on the genes associated with the DMRs and visualise the pathways (and upstream regulators?)

Possible directions:
- Run a direct DM comparison to identify DMRs without the need for a comparator
- Look for shared but discordant between the two groups
- Link to DE and look for a change in the concordance relative to 'all DMRs' (try this with TSS only)
"""

def pathway_involvement_heatmap_by_p(
        p_dat,
        pathway_order,
        orientation='vertical',
        vmin=None,
        vmax=None,
):

    n = p_dat.columns.levels[0].size + 1

    if orientation == 'vertical':
        figsize = (7., 7.)
        gs_kw = dict(
            left=0.7,
            right=0.92,
            top=0.998,
            bottom=0.15,
            wspace=0.1,
            hspace=0.1,
            width_ratios=[1] * (n - 1) + [.5],
        )
        ncols = n
        nrows = 5
    else:
        figsize = (8., 7.)
        ncols = 5
        nrows = n
        raise NotImplementedError("TODO: define gs_kw.")

    fig = plt.figure(figsize=figsize)

    # set up axis grid
    gs = gridspec.GridSpec(
        ncols=ncols,
        nrows=nrows,
        **gs_kw
    )
    if orientation == 'vertical':
        cax = fig.add_subplot(gs[1:-1, -1])
    else:
        cax = fig.add_subplot(gs[-1, 1:-1])

    axs = []

    for i, col in enumerate(p_dat.columns.levels[0]):
        if orientation == 'vertical':
            ax = fig.add_subplot(gs[:, i])
        else:
            ax = fig.add_subplot(gs[i, :])
        axs.append(ax)

        this_dat = p_dat.loc[pathway_order, col]

        if orientation != 'vertical':
            this_dat = this_dat.transpose()

        sns.heatmap(
            this_dat,
            mask=this_dat.isnull(),
            cmap='YlOrRd',
            linewidths=.2,
            linecolor='w',
            vmin=vmin,
            vmax=vmax,
            yticklabels=False,
            cbar=i == 0,
            cbar_ax=cax,
            cbar_kws={"orientation": orientation},
            ax=ax,
        )
        if orientation == 'vertical':
            cax.set_ylabel('$-\log(p)$')
            ax.xaxis.set_ticklabels(p_dat[col].columns, rotation=90.)
            ax.set_xlabel(col)
        else:
            cax.set_xlabel('$-\log(p)$')
            ax.yaxis.set_ticklabels(p_dat[col].columns, rotation=0.)
            ax.set_ylabel(col)

        ax.set_facecolor('0.6')

    # colorbar outline
    cbar = axs[0].collections[0].colorbar
    cbar.outline.set_linewidth(1.)
    cbar.outline.set_edgecolor('k')

    if orientation == 'vertical':
        axs[0].set_yticks(0.5 + np.arange(len(pathway_order)))
        axs[0].set_yticklabels(pathway_order[::-1], rotation=0, fontsize=7)
    else:
        # TODO: check this is correct
        axs[-1].set_xticks(0.5 + np.arange(len(pathway_order)))
        axs[-1].set_xticklabels(pathway_order, rotation=90, fontsize=7)

    return {
        'figure': fig,
        'axs': axs,
        'cax': cax,
        'cbar': cbar,
    }


def plot_clustermap_tpm_levels(
        rna_tpm,
        gois,
        goi_hypo,
        goi_hyper,
        min_tpm=1.,
        log=True,
        c_hyper=consts.METHYLATION_DIRECTION_COLOURS['hyper'],
        c_hypo=consts.METHYLATION_DIRECTION_COLOURS['hypo'],
        **kwargs
):
    g_hyper = reference_genomes.gene_symbol_to_ensembl(goi_hyper).dropna()
    g_hypo = reference_genomes.gene_symbol_to_ensembl(goi_hypo).dropna()

    all_genes = reference_genomes.gene_symbol_to_ensembl(gois).dropna()
    all_genes = all_genes.loc[(~all_genes.index.duplicated()) & (~all_genes.duplicated())]

    row_colors = pd.DataFrame('0.5', index=all_genes.values, columns=['Group'])
    row_colors.loc[row_colors.index.intersection(g_hypo)] = c_hypo
    row_colors.loc[row_colors.index.intersection(g_hyper)] = c_hyper
    col_colors = pd.DataFrame(c_hyper, index=rna_tpm.columns,
                              columns=['Group', 'Cell type'])
    col_colors.loc[rna_meta.patient_id.isin(groups['Hypo']), 'Group'] = c_hypo
    col_colors.loc[rna_meta.type == 'GBM', 'Cell type'] = '0.2'
    col_colors.loc[rna_meta.type == 'iNSC', 'Cell type'] = '0.8'

    if len(col_colors['Cell type'].unique()) < 2:
        col_colors.drop('Cell type', axis=1, inplace=True)

    dat_for_plot = rna_tpm.reindex(all_genes).dropna(how='all', axis=0)
    # filter genes that aren't expressed
    dat_for_plot = dat_for_plot.loc[(dat_for_plot > min_tpm).sum(axis=1) > 2]

    if log:
        dat_for_plot = np.log10(dat_for_plot + 0.01)

    cg = sns.clustermap(
        dat_for_plot,
        yticklabels=False,
        row_colors=row_colors,
        col_colors=col_colors,
        figsize=(8.5, 7.2),
        **kwargs
    )
    plt.setp(cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90)
    cg.ax_heatmap.yaxis.label.set_visible(False)
    cg.gs.update(bottom=0.28, top=0.98, left=0.05)
    return cg


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
            ix = de_res[pid].index.intersection(reference_genomes.gene_symbol_to_ensembl(g_arr).dropna().values)
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
        pids=consts.PIDS,
        relation_filter=None
):
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


if __name__ == '__main__':
    """
    Use the DMR bias results to lookup into DE results (and vice versa?)
    """
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

    DE_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'de')
    de_params = consts.DE_PARAMS

    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()
    DMR_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'dmr')

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

    # load chromosome lengths
    chroms = [str(t) for t in range(1, 23)]
    fa_fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human/ensembl/GRCh38.release87/fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    )
    cl = genomics.feature_lengths_from_fasta(fa_fn, features=chroms)

    chrom_length = collections.OrderedDict()
    for k in chroms:
        if chr_prefix:
            new_k = "chr%s" % k
        else:
            new_k = k
        chrom_length[new_k] = cl[k]

    ## Identify group-specific DMRs

    groups = {
        'Hypo': ['019', '030', '031', '017'],
        'Hyper': ['018', '050', '054', '061', '026', '052']
    }
    group_ind = setops.groups_to_ind(pids, groups)
    groups_inv = dictionary.complement_dictionary_of_iterables(groups, squeeze=True)

    dmr_by_member = [dmr_res_all[pid].keys() for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*dmr_by_member)

    venn_sets_by_group = setops.full_partial_unique_other_sets_from_groups(pids, groups)
    dmr_groups = {}
    genes_from_dmr_groups = {}
    for grp in groups:
        # generate bar chart showing number / pct in each direction (DM)
        this_sets = venn_sets_by_group['full'][grp] + venn_sets_by_group['partial'][grp]
        this_dmrs = sorted(setops.reduce_union(*[venn_set[k] for k in this_sets]))
        dmr_groups[grp] = this_dmrs

        if generate_bw:
            for pid in groups[grp]:
                for_bw = dict([(t, dmr_res_all[pid][t]) for t in this_dmrs if t in dmr_res_all[pid]])
                # where are these DMRs? write bigwig
                fn = os.path.join(outdir, "%s_%s_dmrs.bw" % (grp, pid))
                dmr_values_to_bigwig.write_bigwig(
                    for_bw,
                    dmr_res_s1.clusters,
                    chrom_length,
                    fn,
                    chr_prefix=chr_prefix
                )
        # combine with DE to assess bias (if any)
        this_genes = sorted(
            set(
                [u[0] for u in setops.reduce_union(*[dmr_res_s1.clusters[t].genes for t in this_dmrs])]
            )
        )
        genes_from_dmr_groups[grp] = this_genes

    ## Identify DMRs shared between the groups but discordant

    # Go through all 'mixed' DMRs (across groups) and assemble delta M values
    shared_dmrs = collections.defaultdict(lambda: collections.defaultdict(dict))
    for k in venn_sets_by_group['mixed']:
        this_members = setops.key_to_members(k, pids)
        this_cluster_ids = venn_set[k]
        # get effect size / direction
        for pid in this_members:
            this_grp = [g for g in groups if pid in groups[g]][0]
            for c in venn_set[k]:
                shared_dmrs[c][this_grp][pid] = dmr_res_all[pid][c]['median_change']

    # Reiterate and search for DMRs with discordant directions
    shared_dmrs_discordant = {}
    discord_within_group = {}
    for c, d in shared_dmrs.items():
        this_dmr = {}
        keep = False
        for g in groups:
            tmp = np.sign(d[g].values())
            if (tmp != tmp[0]).any():
                discord_within_group[c] = d
                # don't bother continuing with this DMR
                break
            else:
                # all signs the same in this group
                this_dmr[g] = tmp[0]
        if len(this_dmr) >= 2:
            # is there a discrepancy?
            tmp = np.array(this_dmr.values())
            keep = (tmp != tmp[0]).any()

        if keep:
            shared_dmrs_discordant[c] = d

    # some of these clusters are only discordant because one of the groups has only a single patient
    # let's require at least 2 patients in each group
    shared_dmrs_discordant_gte2 = dict([
        (k, v) for k, v in shared_dmrs_discordant.items() if (np.array([len(t) for t in v.values()]) > 1).all()
    ])

    # for completeness, bar chart showing DMR direction in discordant DMRs
    dmr_res_discordant = collections.defaultdict(dict)

    for cid, res in shared_dmrs_discordant.items():
        this_pids = setops.reduce_union(*res.values())
        for pid in this_pids:
            dmr_res_discordant[pid][cid] = dmr_res_s1[pid].results[cid]

    plt_dict = addd.dm_direction_bar_plot(dmr_res_discordant, keys=pids)
    plt_dict['fig'].savefig(os.path.join(outdir, "discordant_dmr_directions.png"), dpi=200)

    dmr_res_discordant_gte2 = collections.defaultdict(dict)

    for cid, res in shared_dmrs_discordant_gte2.items():
        this_pids = setops.reduce_union(*res.values())
        for pid in this_pids:
            dmr_res_discordant_gte2[pid][cid] = dmr_res_s1[pid].results[cid]

    plt_dict = addd.dm_direction_bar_plot(dmr_res_discordant_gte2, keys=pids)
    plt_dict['fig'].savefig(os.path.join(outdir, "discordant_dmr_directions_gte2.png"), dpi=200)


    genes_discordant, rels_discordant = get_genes_relations(shared_dmrs_discordant, dmr_res_s1.clusters)
    genes_discordant_gte2, rels_discordant_gte2 = get_genes_relations(shared_dmrs_discordant_gte2, dmr_res_s1.clusters)
    genes_hypo, rels_hypo = get_genes_relations(dmr_groups['Hypo'], dmr_res_s1.clusters)
    genes_hyper, rels_hyper = get_genes_relations(dmr_groups['Hyper'], dmr_res_s1.clusters)

    # export for IPA
    ix = sorted(setops.reduce_union(genes_hyper, genes_hypo, genes_discordant, genes_discordant_gte2))
    df_for_ipa = pd.DataFrame(1., columns=['Hypo', 'Hyper', 'Shared discordant', 'Shared discordant GTE2'], index=ix)
    df_for_ipa.loc[genes_hypo, 'Hypo'] = 0.01
    df_for_ipa.loc[genes_hyper, 'Hyper'] = 0.01
    df_for_ipa.loc[genes_discordant, 'Shared discordant'] = 0.01
    df_for_ipa.loc[genes_discordant_gte2, 'Shared discordant GTE2'] = 0.01
    df_for_ipa.to_excel(os.path.join(outdir, "group_specific_gene_lists_for_ipa.xlsx"))

    # repeat but only keeping TSS-related genes
    genes_hypo_tss, rels_hypo_tss = get_genes_relations(dmr_groups['Hypo'], dmr_res_s1.clusters, relation_filter=relations_tss)
    genes_hyper_tss, rels_hyper_tss = get_genes_relations(dmr_groups['Hyper'], dmr_res_s1.clusters, relation_filter=relations_tss)

    ix = sorted(setops.reduce_union(genes_hyper_tss, genes_hypo_tss))
    df_for_ipa = pd.DataFrame(1., columns=['Hypo', 'Hyper'], index=ix)
    df_for_ipa.loc[genes_hypo_tss, 'Hypo'] = 0.01
    df_for_ipa.loc[genes_hyper_tss, 'Hyper'] = 0.01
    df_for_ipa.to_excel(os.path.join(outdir, "group_specific_gene_lists_for_ipa_tss_only.xlsx"))

    # discordant AND TSS only
    genes_discordant_tss, rels_discordant_tss = get_genes_relations(shared_dmrs_discordant, dmr_res_s1.clusters, relation_filter=relations_tss)
    genes_discordant_tss_gte2, rels_discordant_tss_gte2 = get_genes_relations(
        shared_dmrs_discordant_gte2, dmr_res_s1.clusters, relation_filter=relations_tss)

    # plot heatmap and generate Cytoscape session for the IPA results

    ipa_res = ipa.load_raw_reports(IPA_PATHWAY_DIR, "{0}_{1}.txt", ['hyper', 'hypo'], ['all_relations', 'tss'])
    pathways_significant = setops.reduce_union(*[
        v.loc[(v['-logp'] >= plogalpha)].index for v in ipa_res.values()
    ])

    # heatmap
    from scripts.hgic_final import ipa_results_s1_s2 as irss

    dd = irss.generate_plotting_structures(ipa_res, pathways_significant, plogalpha)
    de_all_p = dd['p']
    de_all_z = dd['z']
    de_all_in = dd['in']

    # plot 1) P values, ordered by sum of -log(P)
    p_order = de_all_p.sum(axis=1).sort_values(ascending=False).index

    new_cols = [
        de_all_p.columns.levels[0].str.capitalize().tolist(),
        ['All', 'TSS'],
    ]
    de_all_p.columns = pd.MultiIndex.from_product(new_cols, names=['Comparison', 'Inclusion'])

    plot_dict = pathway_involvement_heatmap_by_p(
        de_all_p,
        p_order,
    )
    plot_dict['figure'].savefig(os.path.join(outdir, "ipa_heatmap.png"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "ipa_heatmap.tiff"), dpi=200)

    # cytoscape
    min_edge_count = 4
    cyto_colours = consts.METHYLATION_DIRECTION_COLOURS
    try:
        cy_session = cyto.CytoscapeSession()
    except cyto.cyrest.requests.ConnectionError:
        logger.exception("Unable to launch Cytoscape session: is it running?")
    else:
        cyto_nets = {}

        for k in ['all_relations', 'tss']:
            this_dat = dict([
                (k2, ipa_res[(k2, k)].loc[ipa_res[(k2, k)]['-logp'] >= plogalpha]) for k2 in ['hyper', 'hypo']
            ])
            gg = ipa.nx_graph_from_ipa_multiple(this_dat, name='IPA DMR groups %s' % k, min_edge_count=min_edge_count)
            this_net = cy_session.add_networkx_graph(gg, name=gg.name)
            cyto_nets[gg.name] = this_net

            # formatting
            this_net.passthrough_node_label('name')
            this_net.passthrough_node_size_linear('n_genes')
            this_net.passthrough_edge_width_linear('n_genes', xmin=min_edge_count, ymin=0.4, ymax=5)
            this_net.set_node_border_width(0.)
            this_net.set_edge_colour('#b7b7b7')
            this_net.set_node_fill_colour('#ffffff')
            this_net.set_node_transparency(255)

            this_net.node_pie_charts(['hyper', 'hypo'], colours=[cyto_colours[t] for t in ['hyper', 'hypo']])

        layout_kwargs = dict(
            EdgeAttribute='n_genes',
            defaultSpringLength=50,
            defaultSpringCoefficient=50,
        )
        cy_session.apply_layout(**layout_kwargs)
        # in practice, best layout can be achieved by manual movement

        cy_session.cy_cmd.session.save(os.path.join(outdir, "ipa_cytoscape.cys"))

    # distribution of cluster locations in the two groups (TSS, exon, etc.)
    relation_priority = [
        'TSS200',
        'TSS1500',
        '1stExon',
        "3'UTR",
        "5'UTR",
        # 'ExonBnd',
        'Body'
    ]

    rel_counts = pd.DataFrame(index=relation_priority)
    rel_counts.insert(
        0,
        'Background',
        pd.Index(get_genes_relations(dmr_res_s1.clusters.keys(), dmr_res_s1.clusters)[1]).value_counts().loc[rel_counts.index]
    )
    rel_counts.insert(
        0,
        'Hypo',
        pd.Index(rels_hypo).value_counts().reindex(rel_counts.index)
    )
    rel_counts.insert(
        0,
        'Hyper',
        pd.Index(rels_hyper).value_counts().reindex(rel_counts.index)
    )
    rel_counts.insert(0, 'Discordant', pd.Index(rels_discordant).value_counts().reindex(rel_counts.index))
    rel_counts.insert(0, 'Discordant GTE2', pd.Index(rels_discordant_gte2).value_counts().reindex(rel_counts.index))
    rel_counts = rel_counts.fillna(0).astype(int)

    fig = plt.figure(figsize=(5.5, 3.3))
    ax = fig.add_subplot(111)
    sns.heatmap(rel_counts.divide(rel_counts.sum(), axis=1) * 100., cmap='Reds')
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    cax = [t for t in fig.get_axes() if t is not ax][0]
    cax.set_title('% of DMRs')

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "relation_to_gene_pct_dmrs.png"), dpi=200)

    # how do these compare to the group-specific DMR genes?
    set_colours_dict = {
        'Hypo': 'g',
        'Hyper': 'r',
        'Discordant': 'b'
    }

    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    set_labels = genes_from_dmr_groups.keys()
    venn.venn_diagram(
        *genes_from_dmr_groups.values(),
        set_labels=set_labels,
        set_colors=[set_colours_dict[t] for t in set_labels],
        ax=ax
    )
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "genes_from_group_spec_dmrs.png"), dpi=200)

    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    set_labels = ['Hypo', 'Hyper']
    venn.venn_diagram(
        genes_hypo_tss, genes_hyper_tss,
        set_labels=set_labels,
        set_colors=[set_colours_dict[t] for t in set_labels],
        ax=ax
    )
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "genes_from_group_spec_dmrs_tss_only.png"), dpi=200)

    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    set_labels = (genes_from_dmr_groups.keys() + ['Discordant'])
    venn.venn_diagram(
        *(genes_from_dmr_groups.values() + [list(genes_discordant)]),
        set_labels=set_labels,
        set_colors=[set_colours_dict[t] for t in set_labels],
        ax=ax
    )
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "genes_from_group_spec_and_discordant_dmrs.png"), dpi=200)

    # how do these compare to the group-specific DMR genes (GTE2)?
    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    venn.venn_diagram(
        *(genes_from_dmr_groups.values() + [list(genes_discordant_gte2)]),
        set_labels=set_labels,
        set_colors=[set_colours_dict[t] for t in set_labels],
        ax=ax
    )
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "genes_from_group_spec_and_discordant_dmrs_gte2.png"), dpi=200)

    # discordant and TSS only
    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    set_labels = ['Hypo', 'Hyper', 'Discordant']
    venn.venn_diagram(
        genes_hypo_tss, genes_hyper_tss, list(genes_discordant_tss),
        set_labels=set_labels,
        set_colors=[set_colours_dict[t] for t in set_labels],
        ax=ax
    )
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "genes_from_group_spec_and_discordant_dmrs_tss_only.png"), dpi=200)

    # discordant and TSS only (GTE2)
    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    set_labels = ['Hypo', 'Hyper', 'Discordant']
    venn.venn_diagram(
        genes_hypo_tss, genes_hyper_tss, list(genes_discordant_tss_gte2),
        set_labels=set_labels,
        set_colors=[set_colours_dict[t] for t in set_labels],
        ax=ax
    )
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "genes_from_group_spec_and_discordant_dmrs_tss_only_gte2.png"), dpi=200)

    ## April 2019: New DE/DM matching process for IPA
    # Idea here: rather than just looking at genes corresponding to group-specific DMRs, we make the requirements more
    # stringent. For each Venn set (e.g. 018, 054, 052 - hyper group), we require DE genes in the same patients.
    # Simplest approach is to use the joint_de_dmr dataframes, which have already been combined.

    # all relations

    tmp = get_de_dmr_groups(joint_de_dmr_s1, dmr_res_s1.clusters)
    de_dmr_de_fdr = tmp['de_FDR']
    de_dmr_de_logfc = tmp['de_logFC']

    # export these for IPA analysis
    df_for_ipa = pd.DataFrame(
        index=sorted(setops.reduce_union(*[t.index for t in de_dmr_de_logfc.values()])),
        columns=reduce(operator.add, [["%s_logFC" % pid, "%s_FDR" % pid] for pid in pids])
    )
    for grp in groups:
        for pid in groups[grp]:
            this_logfc = de_dmr_de_logfc[grp][pid].dropna()
            this_fdr = de_dmr_de_fdr[grp][pid].dropna()
            df_for_ipa.loc[this_logfc.index, "%s_logFC" % pid] = this_logfc
            df_for_ipa.loc[this_fdr.index, "%s_FDR" % pid] = this_fdr
    df_for_ipa.to_excel(os.path.join(outdir, "group_specific_de_for_ipa.xlsx"))

    # TSS relation only

    tmp = get_de_dmr_groups(joint_de_dmr_s1, dmr_res_s1.clusters, relation_filter=['TSS200', 'TSS1500'])
    de_dmr_de_fdr = tmp['de_FDR']
    de_dmr_de_logfc = tmp['de_logFC']

    # export these for IPA analysis
    df_for_ipa = pd.DataFrame(
        index=sorted(setops.reduce_union(*[t.index for t in de_dmr_de_logfc.values()])),
        columns=reduce(operator.add, [["%s_logFC" % pid, "%s_FDR" % pid] for pid in pids])
    )
    for grp in groups:
        for pid in groups[grp]:
            this_logfc = de_dmr_de_logfc[grp][pid].dropna()
            this_fdr = de_dmr_de_fdr[grp][pid].dropna()
            df_for_ipa.loc[this_logfc.index, "%s_logFC" % pid] = this_logfc
            df_for_ipa.loc[this_fdr.index, "%s_FDR" % pid] = this_fdr
    df_for_ipa.to_excel(os.path.join(outdir, "group_specific_de_for_ipa_tss.xlsx"))


    # permutation tests for these: is the lack of overlap actually significant?
    # null: pick the same number of distinct DMRs without replacement from the pool
    def permutation_test_dmr_gene_overlap(n_hypo, n_hyper, relation_filter=None, n_iter=1000):
    # relation_filter = None

    # n_iter = 1000
        venn_counts = []
        for i in range(n_iter):
            picked_hypo = np.random.choice(dmr_res_s1.clusters.keys(), n_hypo, replace=False)
            picked_hyper = np.random.choice(list(set(dmr_res_s1.clusters.keys()).difference(picked_hypo)), n_hyper, replace=False)
            if relation_filter is None:
                picked_hyper_genes = setops.reduce_union(*[[x[0] for x in dmr_res_s1.clusters[t].genes] for t in picked_hyper])
                picked_hypo_genes = setops.reduce_union(*[[x[0] for x in dmr_res_s1.clusters[t].genes] for t in picked_hypo])
            else:
                picked_hyper_genes = setops.reduce_union(
                    *[[x[0] for x in dmr_res_s1.clusters[t].genes if x[1] in relation_filter] for t in picked_hyper]
                )
                picked_hypo_genes = setops.reduce_union(
                    *[[x[0] for x in dmr_res_s1.clusters[t].genes if x[1] in relation_filter] for t in picked_hypo]
                )
            _, vc = setops.venn_from_arrays(picked_hypo_genes, picked_hyper_genes)
            venn_counts.append(vc)
        return venn_counts

    n_hypo = len(dmr_groups['Hypo'])
    n_hyper = len(dmr_groups['Hyper'])
    vc_all = permutation_test_dmr_gene_overlap(n_hypo, n_hyper)
    # reduce to single metric
    vc_metric = np.array([
        max(
            x['11'] / float(x['10']) * 100,
            x['11'] / float(x['01']) * 100,
        ) for x in vc_all
    ])
    _, vc = setops.venn_from_arrays(*genes_from_dmr_groups.values())
    our_vc_metric = max(
        vc['11'] / float(vc['10']) * 100.,
        vc['11'] / float(vc['01']) * 100.,
    )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(vc_metric, 20, alpha=0.6, label='Null')
    ax.axvline(our_vc_metric, c='k', label='Our value')
    ax.set_xlabel('Maximum % overlap')
    ax.set_ylabel('Frequency')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "gene_overlap_vs_null_all.png"), dpi=200)

    # again with TSS only
    dmr_groups_hypo = [t for t in dmr_groups['Hypo'] if len([x for x in dmr_res_s1.clusters[t].genes if x[1] in {'TSS200', 'TSS1500'}])]
    dmr_groups_hyper = [t for t in dmr_groups['Hyper'] if len([x for x in dmr_res_s1.clusters[t].genes if x[1] in {'TSS200', 'TSS1500'}])]

    n_hypo = len(dmr_groups_hypo)
    n_hyper = len(dmr_groups_hyper)
    vc_all = permutation_test_dmr_gene_overlap(n_hypo, n_hyper)
    # reduce to single metric
    vc_metric = np.array([
        max(
            x['11'] / float(x['10']) * 100,
            x['11'] / float(x['01']) * 100,
        ) for x in vc_all
    ])
    _, vc = setops.venn_from_arrays(genes_hyper_tss, genes_hypo_tss)
    our_vc_metric = max(
        vc['11'] / float(vc['10']) * 100.,
        vc['11'] / float(vc['01']) * 100.,
    )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(vc_metric, 20, alpha=0.6, label='Null')
    ax.axvline(our_vc_metric, c='k', label='Our value')
    ax.set_xlabel('Maximum % overlap')
    ax.set_ylabel('Frequency')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "gene_overlap_vs_null_tss.png"), dpi=200)

    # Plot gene expr levels (in GIC) of genes that correspond to the group-specific DMRs

    # add 'group' to RNA meta for ease of ordering
    rna_meta_grp = pd.Series('hypo', index=rna_meta.index)
    rna_meta_grp.loc[rna_meta.patient_id.isin(groups['Hyper'])] = 'hyper'
    rna_meta.insert(1, 'dmr_group', rna_meta_grp.loc[rna_meta.index])

    # all GS DMRs
    cg = plot_clustermap_tpm_levels(
        rna_tpm[rna_meta.sort_values(by=['type', 'dmr_group']).index],
        genes_hyper + genes_hypo,
        genes_hypo,
        genes_hyper,
        col_cluster=False,
        z_score=0
    )
    cg.savefig(os.path.join(outdir, "group_spec_genes_all_rels_clustermap.png"), dpi=200)

    # TSS-related only
    cg = plot_clustermap_tpm_levels(
        rna_tpm[rna_meta.sort_values(by=['type', 'dmr_group']).index],
        genes_hyper_tss + genes_hypo_tss,
        genes_hypo,
        genes_hyper,
        col_cluster=False,
        z_score=0
    )
    cg.savefig(os.path.join(outdir, "group_spec_genes_tss_clustermap.png"), dpi=200)

    # discordant
    cg = plot_clustermap_tpm_levels(
        rna_tpm[rna_meta.sort_values(by=['type', 'dmr_group']).index],
        genes_discordant,
        genes_hypo,
        genes_hyper,
        col_cluster=False,
        z_score=0
    )
    cg.savefig(os.path.join(outdir, "discordant_genes_all_rels_clustermap.png"), dpi=200)

    # discordant GTE2
    cg = plot_clustermap_tpm_levels(
        rna_tpm[rna_meta.sort_values(by=['type', 'dmr_group']).index],
        genes_discordant_gte2,
        genes_hypo,
        genes_hyper,
        col_cluster=False,
        z_score=0
    )
    cg.savefig(os.path.join(outdir, "discordant_gte2_genes_all_rels_clustermap.png"), dpi=200)

    # discordant GTE2 TSS only
    cg = plot_clustermap_tpm_levels(
        rna_tpm[rna_meta.sort_values(by=['type', 'dmr_group']).index],
        genes_discordant_tss_gte2,
        genes_hypo,
        genes_hyper,
        col_cluster=False,
        z_score=0
    )
    cg.savefig(os.path.join(outdir, "discordant_gte2_genes_tss_clustermap.png"), dpi=200)

    # extract the DE/DMR results that are only associated with the group-specific DMRs
    tss_cols = ['dmr_TSS1500', 'dmr_TSS200']
    group_specific_de_dmr = dict([
        (
            pid,
            joint_de_dmr_s1[pid].loc[joint_de_dmr_s1[pid].cluster_id.isin(dmr_groups[groups_inv[pid]]), ['de_logFC', 'dmr_median_delta']]
        ) for pid in pids])

    group_specific_de_dmr_tss = dict([
        (
            pid,
            joint_de_dmr_s1[pid].loc[
                joint_de_dmr_s1[pid].cluster_id.isin(dmr_groups[groups_inv[pid]]) & joint_de_dmr_s1[pid][tss_cols].sum(axis=1).astype(bool),
                ['de_logFC', 'dmr_median_delta']
            ]
        ) for pid in pids
    ])

    # table of DE counts (in DMR-linked context)
    de_count_table, de_count_table_pct = tabulate_de_counts_by_direction(
        de_res_s1,
        Hypo=genes_hypo,
        Hyper=genes_hyper
    )
    colours = [consts.METHYLATION_DIRECTION_COLOURS[groups_inv[pid].lower()] for pid in pids]
    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(5., 5.5), sharex=True, sharey=False)
    for i, col in enumerate(de_count_table_pct.columns):
        ax = axs[i]
        ax.bar(range(len(pids)), de_count_table_pct[col], color=colours, edgecolor='k', linewidth=1.)
        ax.set_ylabel("%% %s" % col)
        ax.set_xticks(range(len(pids)))
        ax.set_xticklabels(pids, rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "de_gene_count_by_direction_and_group.png"), dpi=200)

    # again with TSS
    de_count_table_tss, de_count_table_pct_tss = tabulate_de_counts_by_direction(
        de_res_s1,
        Hypo=genes_hypo_tss,
        Hyper=genes_hyper_tss
    )
    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(5., 5.5), sharex=True, sharey=False)
    for i, col in enumerate(de_count_table_pct_tss.columns):
        ax = axs[i]
        ax.bar(range(len(pids)), de_count_table_pct_tss[col], color=colours, edgecolor='k', linewidth=1.)
        ax.set_ylabel("%% %s" % col)
        ax.set_xticks(range(len(pids)))
        ax.set_xticklabels(pids, rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "de_gene_count_by_direction_and_group_tss.png"), dpi=200)

    # again with discordant
    de_count_table_disc, de_count_table_pct_disc = tabulate_de_counts_by_direction(
        de_res_s1,
        Discordant=genes_discordant,
    )

    # plot clustermaps of gene expr / DE logFC for these genes
    all_genes = setops.reduce_union(
        *[[t[1] for t in group_specific_de_dmr[pid].index] for pid in pids]
    )

    tmp = pd.DataFrame(0., index=sorted(all_genes), columns=pids)
    for pid in pids:
        this_res = de_res_s1[pid].loc[de_res_s1[pid]['Gene Symbol'].isin(all_genes)]
        tmp.loc[this_res['Gene Symbol'], pid] = this_res.logFC

    cg = plot_clustermap_tpm_levels(
        rna_tpm[rna_meta.loc[rna_meta.type == 'GBM'].sort_values(by=['dmr_group', 'patient_id']).index],
        all_genes,
        genes_hypo,
        genes_hyper,
        z_score=0,
        metric='correlation',
        method='average'
    )

    ## run through DE genes (per patient) and look at direction distribution
    # look at the direction distribution of genes that correspond to a DMR (full and specific lists)

    tss_cols = ['dmr_TSS1500', 'dmr_TSS200']

    # full list, all relations
    de_linked = dict([
        (
            pid,
            de_res_s1[pid].loc[de_res_s1[pid]['Gene Symbol'].isin(joint_de_dmr_s1[pid].gene)]
        ) for pid in pids
    ])

    de_by_direction = same_de.count_de_by_direction(de_linked)

    plt_dict = same_de.bar_plot(de_linked, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "de_linked_syngeneic_all_rels_directions.png"), dpi=200)

    # full list, TSS only
    de_linked = {}
    for pid in pids:
        ix = joint_de_dmr_s1[pid][tss_cols].sum(axis=1).astype(bool)
        this_genes = joint_de_dmr_s1[pid].loc[ix, 'gene']
        de_linked[pid] = de_res_s1[pid].loc[de_res_s1[pid]['Gene Symbol'].isin(this_genes)]

    de_by_direction = same_de.count_de_by_direction(de_linked)

    plt_dict = same_de.bar_plot(de_linked, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "de_linked_syngeneic_tss_directions.png"), dpi=200)

    # patient-specific DMRs linked to genes, all relations
    spec_ix = setops.specific_features(*[dmr_res_all[pid].keys() for pid in pids])
    dm_specific = dict([
        (
            pid,
            dict([
                (
                    k,
                    dmr_res_all[pid][k]
                ) for k in s
            ])
        ) for pid, s in zip(pids, spec_ix)
    ])

    # manually link these
    dm_specific_genes = {}
    for pid in pids:
        cl_ids = dm_specific[pid].keys()
        dm_specific_genes[pid] = setops.reduce_union(*[[t[0] for t in dmr_res_s1.clusters[c].genes] for c in cl_ids])

    de_linked_spec = dict([
        (
            pid,
            de_res_s1[pid].loc[de_res_s1[pid]['Gene Symbol'].isin(dm_specific_genes[pid])]
        ) for pid in pids
    ])

    de_by_direction_spec = same_de.count_de_by_direction(de_linked_spec)

    plt_dict = same_de.bar_plot(de_linked_spec, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "de_specific_linked_syngeneic_all_rels_directions.png"), dpi=200)

    # patient-specific DMRs linked to genes, TSS only
    de_linked_spec_tss = {}
    for pid in pids:
        ix = joint_de_dmr_s1[pid][tss_cols].sum(axis=1).astype(bool)
        this_genes = joint_de_dmr_s1[pid].loc[ix, 'gene']
        de_linked_spec_tss[pid] = de_linked_spec[pid].loc[de_linked_spec[pid]['Gene Symbol'].isin(this_genes)]

    de_by_direction_spec_tss = same_de.count_de_by_direction(de_linked_spec_tss)

    plt_dict = same_de.bar_plot(de_linked_spec_tss, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "de_specific_linked_syngeneic_tss_directions.png"), dpi=200)

    groups = {
        'Hypo': ['019', '030', '031', '017'],
        'Hyper': ['018', '050', '054', '061', '026', '052']
    }
    group_ind = setops.groups_to_ind(pids, groups)
    groups_inv = dictionary.complement_dictionary_of_iterables(groups, squeeze=True)
    # upset plotting colours
    subgroup_set_colours = {
        'Hypo full': '#427425',
        'Hyper full': '#b31500',
        'Hypo partial': '#c5e7b1',
        'Hyper partial': '#ffa599',
        'Expanded core': '#4C72B0',
        'Specific': '#f4e842',
    }

    # start with an UpSet to represent this
    dmr_by_member = [dmr_res_all[pid].keys() for pid in pids]
    venn_set, venn_ct = setops.venn_from_arrays(*dmr_by_member)

    # add null set manually from full DMR results
    dmr_id_all = setops.reduce_union(*venn_set.values())
    k_null = ''.join(['0'] * len(pids))
    venn_set[k_null] = list(set(dmr_res_s1.clusters.keys()).difference(dmr_id_all))
    venn_ct[k_null] = len(venn_set[k_null])

    upset = venn.upset_plot_with_groups(
        dmr_by_member,
        pids,
        group_ind,
        subgroup_set_colours,
        min_size=10,
        n_plot=30,
        default_colour='gray'
    )
    upset['axes']['main'].set_ylabel('Number of DMRs in set')
    upset['axes']['set_size'].set_xlabel('Number of DMRs in single comparison')
    upset['figure'].savefig(os.path.join(outdir, "upset_dmr_by_direction_group.png"), dpi=200)

    # export DMR data to bigwig
    venn_sets_by_group = setops.full_partial_unique_other_sets_from_groups(pids, groups)
    for grp in groups:
        # generate bar chart showing number / pct in each direction (DM)
        this_sets = venn_sets_by_group['full'][grp] + venn_sets_by_group['partial'][grp]
        this_dmrs = dmr_groups[grp]
        this_genes = genes_from_dmr_groups[grp]

        this_results = {}
        for pid in groups[grp]:
            this_results[pid] = dict([(t, dmr_res_all[pid][t]) for t in this_dmrs if t in dmr_res_all[pid]])
            # where are these DMRs? write bigwig
            fn = os.path.join(outdir, "%s_%s_dmrs.bw" % (grp, pid))
            dmr_values_to_bigwig.write_bigwig(
                this_results[pid],
                dmr_res_s1.clusters,
                chrom_length,
                fn,
                chr_prefix=chr_prefix
            )

        # bar chart showing DMR direction
        plt_dict = addd.dm_direction_bar_plot(this_results, groups[grp])
        plt_dict['fig'].savefig(os.path.join(outdir, "dmr_direction_by_group_%s.png" % grp), dpi=200)

        # combine with DE to assess bias (if any)
        this_results = {}
        for pid in pids:
            this_results[pid] = de_res_s1[pid].loc[de_res_s1[pid]['Gene Symbol'].isin(this_genes)]
        plt_dict = same_de.bar_plot(this_results, keys=pids)
        plt_dict['fig'].savefig(os.path.join(outdir, "de_direction_by_dm_group_%s.png" % grp), dpi=200)

        # same again, but density
        fig, axs = plt.subplots(nrows=1, ncols=len(groups[grp]), sharex=True, figsize=(11 / 6. * len(groups[grp]), 2.5))
        for i, pid in enumerate(groups[grp]):
            xi, y = addd.direction_kde(this_results[pid].logFC)
            axs[i].fill_between(xi[xi < 0], y[xi < 0], facecolor=consts.METHYLATION_DIRECTION_COLOURS['hypo'], edgecolor='k', alpha=0.6)
            axs[i].fill_between(xi[xi > 0], y[xi > 0], facecolor=consts.METHYLATION_DIRECTION_COLOURS['hyper'], edgecolor='k', alpha=0.6)
            axs[i].axvline(0., color='k', linestyle=':')
            # axs[i].yaxis.set_visible(False)
            axs[i].set_yticks([])
            axs[i].set_title(pid)
        axs[0].set_ylabel('Density')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "de_genes_from_dmr_%s_logfc_kde.png" % grp), dpi=200)

    # grouped bar chart showing dmr direction (in one figure)
    colours = pd.Series(
        {
            'Hyper': consts.METHYLATION_DIRECTION_COLOURS['hyper'],
            'Hypo': consts.METHYLATION_DIRECTION_COLOURS['hypo']
        }
    )
    group_order = pd.Series([len(t) for t  in groups.values()], index=groups.keys()).sort_values(ascending=False).index

    fig = plt.figure(figsize=(3.5, 4))
    gs = plt.GridSpec(nrows=2, ncols=2, width_ratios=[len(groups[grp]) for grp in group_order])
    sharey = None

    for i, grp in enumerate(group_order):
        this_sets = venn_sets_by_group['full'][grp] + venn_sets_by_group['partial'][grp]
        this_dmrs = sorted(setops.reduce_union(*[venn_set[k] for k in this_sets]))
        this_results = {}
        for pid in groups[grp]:
            this_results[pid] = dict([(t, dmr_res_all[pid][t]) for t in this_dmrs if t in dmr_res_all[pid]])

        for_plot = addd.count_dm_by_direction(this_results, pids=groups[grp])
        for_plot_pct = for_plot.divide(for_plot.sum(), axis=1) * 100.

        ax = fig.add_subplot(gs[1, i], sharey=sharey)
        if sharey is None:
            sharey = ax
        if i == 0:
            ax.set_ylabel('Number DMRs', fontsize=12)
        else:
            ax.yaxis.set_visible(False)

        bar.stacked_bar_chart(for_plot, colours=colours, width=0.8, legend=False, ax=ax)
        plt.setp(ax.xaxis.get_ticklabels(), rotation=90, fontsize=12)
        plt.setp(ax.yaxis.get_ticklabels(), fontsize=12)

        ax = fig.add_subplot(gs[0, i])
        bar.stacked_bar_chart(for_plot_pct, colours=colours, width=0.8, legend=False, ax=ax)
        ax.xaxis.set_ticklabels([])
        plt.setp(ax.yaxis.get_ticklabels(), fontsize=12)
        ax.set_ylim([0, 100])
        if i == 0:
            ax.set_ylabel('% DMRs', fontsize=12)
        else:
            ax.yaxis.set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_direction_all_groups.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "dmr_direction_all_groups.tiff"), dpi=200)
