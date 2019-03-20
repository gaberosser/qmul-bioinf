from plotting import bar, common, pie, polar, venn
from methylation import loader, dmr, process
import pandas as pd
from stats import nht
from utils import output, setops, genomics, log, dictionary
import multiprocessing as mp
import os
import collections
import pickle
import numpy as np
from scipy import stats, cluster
import matplotlib
from matplotlib import pyplot as plt, patches
from matplotlib.colors import Normalize
from matplotlib import cm
from sklearn.neighbors import KernelDensity
import seaborn as sns

import references
from scripts.hgic_final import \
    two_strategies_grouped_dispersion as tsgd, \
    two_strategies_combine_de_dmr as tscd, \
    consts
from scripts.dmr_direction_bias_story import \
    analyse_dmr_direction_and_distribution as addd, \
    same_process_applied_to_de as same_de
from rnaseq import loader as rnaseq_loader
from integrator import rnaseq_methylationarray
from scripts.methylation import dmr_values_to_bigwig

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
logger = log.get_console_logger()

"""
Here we seek to identify DMRs that distinguish the hypo- and hypermethylated groups.

This is initially carried out directly using the GIC-iNSC results.

We then query these DMRs to determine whether they have different location distributions.

Possible directions:
- Run a direct DM comparison to identify DMRs without the need for a comparator
- Look for shared but discordant between the two groups
- Link to DE and look for a change in the concordance relative to 'all DMRs' (try this with TSS only)
"""

if __name__ == '__main__':
    """
    Use the DMR bias results to lookup into DE results (and vice versa?)
    """
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

    # look up probes / genes
    genes_discordant_gte2 = setops.reduce_union(
        *[zip(*dmr_res_s1.clusters[t].genes)[0] for t in shared_dmrs_discordant_gte2 if dmr_res_s1.clusters[t].genes]
    )
    probes_discordant_gte2 = reduce(lambda x, y: x+y, [dmr_res_s1.clusters[t].pids for t in shared_dmrs_discordant_gte2])

    # how do these compare to the group-specific DMR genes?
    venn.venn_diagram(
        *(genes_from_dmr_groups.values() + [list(genes_discordant_gte2)]),
        set_labels=(genes_from_dmr_groups.keys() + ['Discordant'])
    )

    # distribution of probe/region locations in the two groups

    raise StopIteration

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
    dmr_groups = {}
    genes_from_dmr_groups = {}
    for grp in groups:
        # generate bar chart showing number / pct in each direction (DM)
        this_sets = venn_sets_by_group['full'][grp] + venn_sets_by_group['partial'][grp]
        this_dmrs = sorted(setops.reduce_union(*[venn_set[k] for k in this_sets]))
        dmr_groups[grp] = this_dmrs

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
        this_genes = sorted(
            set(
                [u[0] for u in setops.reduce_union(*[dmr_res_s1.clusters[t].genes for t in this_dmrs])]
            )
        )
        genes_from_dmr_groups[grp] = this_genes

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
        fig.savefig(os.path.join(outdir, "de_genes_from_dmr_groups_logfc_kde.png"), dpi=200)

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
            ax.scatter(this_joint.de_logFC, this_joint.dmr_median_delta, c=consts.METHYLATION_DIRECTION_COLOURS[grp.lower()])
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