from plotting import bar, common, pie, polar, venn
from methylation import loader, dmr, process
import pandas as pd
from stats import nht
from utils import output, setops, genomics, log, ipa, dictionary
from cytoscape import cyto
import multiprocessing as mp
import os
import collections
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
from rnaseq import loader as rnaseq_loader
from integrator import rnaseq_methylationarray
from scripts.methylation import dmr_values_to_bigwig

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
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
    rna_obj = rnaseq_loader.load_by_patient(pids, include_control=False)
    rna_obj.filter_samples(rna_obj.meta.index.isin(consts.S1_RNASEQ_SAMPLES))
    tmp = rna_obj.data
    rna_cpm = tmp.divide((tmp + 1).sum(axis=0), axis=1) * 1e6

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
    genes_hypo_tss, rels_hypo_tss = get_genes_relations(dmr_groups['Hypo'], dmr_res_s1.clusters, relation_filter=['TSS200', 'TSS1500'])
    genes_hyper_tss, rels_hyper_tss = get_genes_relations(dmr_groups['Hyper'], dmr_res_s1.clusters, relation_filter=['TSS200', 'TSS1500'])

    ix = sorted(setops.reduce_union(genes_hyper_tss, genes_hypo_tss))
    df_for_ipa = pd.DataFrame(1., columns=['Hypo', 'Hyper'], index=ix)
    df_for_ipa.loc[genes_hypo_tss, 'Hypo'] = 0.01
    df_for_ipa.loc[genes_hyper_tss, 'Hyper'] = 0.01
    df_for_ipa.to_excel(os.path.join(outdir, "group_specific_gene_lists_for_ipa_tss_only.xlsx"))

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
    cy_session = cyto.CytoscapeSession()
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

    # look at the direction distribution of genes that correspond to a DMR (full and specific lists)
    joint_de_dmr_s1 = rnaseq_methylationarray.compute_joint_de_dmr(dmr_res_s1, de_res_s1)

    # run through DE genes (per patient) and look at direction distribution
    # full list
    de_linked = dict([
        (
            pid,
            de_res_s1[pid].loc[de_res_s1[pid]['Gene Symbol'].isin(joint_de_dmr_s1[pid].gene)]
        ) for pid in pids
    ])

    de_by_direction = same_de.count_de_by_direction(de_linked)

    plt_dict = same_de.bar_plot(de_linked, pids)
    plt_dict['fig'].tight_layout()
    plt_dict['fig'].savefig(os.path.join(outdir, "de_linked_syngeneic_full_list_directions.png"), dpi=200)

    # patient-specific DMRs linked to genes
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
    plt_dict['fig'].savefig(os.path.join(outdir, "de_specific_linked_syngeneic_full_list_directions.png"), dpi=200)

    # hypo/hyper DMRs in the hypo group linked (SEPARATELY) to DE

    groups = {
        'Hypo': ['019', '030', '031', '017'],
        'Hyper': ['018', '050', '054', '061', '026', '052']
    }
    group_ind = setops.groups_to_ind(pids, groups)
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

    fig = plt.figure(figsize=(5, 5))
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
            ax.set_ylabel('Number DMRs')
        else:
            ax.yaxis.set_visible(False)

        bar.stacked_bar_chart(for_plot, colours=colours, width=0.8, legend=False, ax=ax)
        ax = fig.add_subplot(gs[0, i])
        bar.stacked_bar_chart(for_plot_pct, colours=colours, width=0.8, legend=False, ax=ax)
        if i == 0:
            ax.set_ylabel('% DMRs')
        else:
            ax.yaxis.set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_direction_all_groups.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "dmr_direction_all_groups.tiff"), dpi=200)

    # is there much overlap in the gene sets between the two groups?
    fig = plt.figure(figsize=(5, 3))
    ax = fig.add_subplot(111)
    venn.venn_diagram(genes_from_dmr_groups['Hyper'], genes_from_dmr_groups['Hypo'], set_labels=['Hyper', 'Hypo'], ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "genes_from_dmr_groups_venn.png"), dpi=200)

    # no, but if we look at the intersection genes, are they in different directions (DE) between the two groups?
    groups_inv = dictionary.complement_dictionary_of_iterables(groups, squeeze=True)
    in_both = setops.reduce_intersection(*genes_from_dmr_groups.values())
    in_both_ens = references.gene_symbol_to_ensembl(in_both)

    # some of these will have no DE results
    tmp = {}
    for pid in pids:
        tmp[pid] = de_res_s1[pid].reindex(in_both_ens).logFC.dropna()
    in_both_ens = in_both_ens[in_both_ens.isin(pd.DataFrame(tmp).index)]

    fig = plt.figure(figsize=(10.5, 4.))
    ax = fig.add_subplot(111)
    for pid in pids:
        this_grp = groups_inv[pid].lower()
        this_de = de_res_s1[pid].reindex(in_both_ens).logFC
        x = np.arange(len(in_both_ens)) + 0.1 * (-1 if this_grp == 'hypo' else 1)
        ax.scatter(x, this_de, c=consts.METHYLATION_DIRECTION_COLOURS[this_grp])
    ax.set_xticks(np.arange(len(in_both_ens)))
    ax.set_xticklabels(in_both_ens.index, rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "de_logfc_genes_from_dmr_group_intersection.png"), dpi=200)


    # generate scatter plots of joint DE and DMR for the two DMR groups
    # FIXME: this is hard coded to the number of patients in each group. Sloppy.

    fig, axs = plt.subplots(nrows=2, ncols=6, sharex=True, sharey=True)
    big_ax = common.add_big_ax_to_subplot_fig(fig)

    subplots_filled = set()
    for i, grp in enumerate(groups):
        for j, pid in enumerate(groups[grp]):
            subplots_filled.add((i, j))
            ax = axs[i, j]

            this_joint = joint_de_dmr_s1[pid].loc[joint_de_dmr_s1[pid].cluster_id.isin(dmr_groups[grp])]
            ax.scatter(
                this_joint.de_logFC,
                this_joint.dmr_median_delta,
                c=consts.METHYLATION_DIRECTION_COLOURS[grp.lower()],
                alpha=0.5
            )

            # highlight TSS-linked (only)
            this_joint_tss = this_joint.loc[(this_joint.dmr_TSS1500 | this_joint.dmr_TSS200).astype(bool)]
            ax.scatter(
                this_joint_tss.de_logFC,
                this_joint_tss.dmr_median_delta,
                facecolors='none',
                edgecolor='k',
                linewidth=1.,
                alpha=0.5
            )

            ax.axhline(0., c='k', linestyle='--')
            ax.axvline(0., c='k', linestyle='--')
            ax.set_title(pid)
    big_ax.set_xlabel(r'DE logFC')
    big_ax.set_ylabel(r'DMR median $\Delta M$')
    fig.tight_layout()

    for i in range(2):
        for j in range(6):
            if (i, j) not in subplots_filled:
                axs[i, j].set_visible(False)
    axs[0,0].set_ylim([-7, 7])
    fig.savefig(os.path.join(outdir, "de_dmr_scatter_from_dmr_groups.png"), dpi=200)

    # what is the extent of 'concordance' here and is it greater than a measure of 'background'?
    fig = plt.figure(figsize=(8, 3.))
    ax = fig.add_subplot(111)

    pct_conc = pd.DataFrame(index=pids, columns=['All', 'TSS'])

    for pid in pids:
        grp = groups_inv[pid]
        this_joint = joint_de_dmr_s1[pid].loc[joint_de_dmr_s1[pid].cluster_id.isin(dmr_groups[grp])]
        # report extent of 'non-bijective mapping'
        logger.info("Patient %s. Of the %d DMR-DE results being considered", pid, this_joint.shape[0])
        logger.info("%d of the DMRs appear more than once", (this_joint.cluster_id.value_counts() > 1).sum())
        logger.info("%d of the DE genes appear more than once", (this_joint.gene.value_counts() > 1).sum())

        pct_conc.loc[pid, 'All'] = pct_concordant(this_joint)

        this_joint = this_joint.loc[(this_joint.dmr_TSS1500 | this_joint.dmr_TSS200).astype(bool)]
        # report extent of 'non-bijective mapping'
        logger.info("Patient %s. Of the %d DMR-DE TSS results being considered", pid, this_joint.shape[0])
        logger.info("%d of the DMRs appear more than once", (this_joint.cluster_id.value_counts() > 1).sum())
        logger.info("%d of the DE genes appear more than once", (this_joint.gene.value_counts() > 1).sum())

        pct_conc.loc[pid, 'TSS'] = pct_concordant(this_joint)

    # background estimate: for each PID, pick DMRs at random with delta M signs matching the number in the
    # group-specific DMRs

    n_iter_bg = 1000
    perm_res = {}

    for pid in pids:
        grp = groups_inv[pid]
        # BG clusters - these are guaranteed to have associated DE genes
        bg_res = joint_de_dmr_s1[pid]
        # split by DMR direction
        bg_hypo = bg_res.loc[bg_res.dmr_median_delta < 0]
        bg_hyper = bg_res.loc[bg_res.dmr_median_delta > 0]
        n_bg_hypo = bg_hypo.shape[0]
        n_bg_hyper = bg_hyper.shape[0]
        # how many do we need to draw?
        this_joint = joint_de_dmr_s1[pid].loc[joint_de_dmr_s1[pid].cluster_id.isin(dmr_groups[grp])]
        n_hypo = (this_joint.dmr_median_delta < 0).sum()
        n_hyper = (this_joint.dmr_median_delta > 0).sum()

        logger.info(
            "Patient %s. Drawing %d hypo and %d hyper DMRs from background (%d hypo / %d hyper).",
            pid,
            n_hypo,
            n_hyper,
            bg_hypo.shape[0],
            bg_hyper.shape[0]
        )

        this_perm_res = []
        for i in range(n_iter_bg):


        # multiple random draws
        pass
