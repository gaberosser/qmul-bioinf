import collections
import itertools
import multiprocessing as mp
import operator
import os
import pickle

import numpy as np
import pandas as pd
from matplotlib import colors
from matplotlib import pyplot as plt
from scipy import stats, special

from integrator import rnaseq_methylationarray
from methylation import dmr, annotation_gene_to_ensembl
from plotting import common, venn
from rnaseq import loader as rnaseq_loader
from scripts.dmr_direction_bias_story import \
    same_process_applied_to_de as same_de
from scripts.hgic_final import \
    two_strategies_grouped_dispersion as tsgd, \
    two_strategies_combine_de_dmr as tscd, \
    analyse_dmrs_s1_direction_distribution as addd, \
    consts
from settings import HGIC_LOCAL_DIR, INTERMEDIATE_DIR
from stats import basic
from utils import output, setops, log, dictionary, reference_genomes

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
    set_labels=('Hypo', 'Hyper'),
    fontsize=16
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
        ## FIXME: font family seems to change from old text to new - related to LaTeX rendering?
        lbl.set_text(
            lbl.get_text()
            + '\n'
            + r'$%d\uparrow$' % len(this_members["%s up" % k])
            + '\n'
            + r'$%d\downarrow$' % len(this_members["%s down" % k]),
        )
    plt.setp(common.get_children_recursive(ax, type_filt=plt.Text), fontsize=fontsize)

    return ax


def fishers_downsample(x, y, n, n_iter=100):
    """
    'Downsample' the data by picking a subset of n data points, then run Fisher's test. Repeat n_iter times.
    :param x:
    :param y:
    :param n:
    :param n_iter:
    :return:
    """
    if len(x) != len(y):
        raise ValueError("Length of x and y must be equal")
    if len(x) < n:
        raise ValueError("len(x) must be greater than or equal to n")

    x = np.array(x)
    y = np.array(y)
    # check the maximum number of permutations and warn if it is too low
    exhaust = False
    n_poss = special.comb(len(x), n, exact=True)
    if n_poss < n_iter:
        logger.warn("Requested %d iterations but the maximum number of permutations is %d, so we'll run an "
                    "exhaustive search.", n_iter, n_poss)
        exhaust = True

    if exhaust:
        it = itertools.combinations(range(len(x)), n)
    else:
        it = (np.random.choice(range(len(x)), n, replace=False) for i in range(n_iter))

    res = []
    for perm in it:
        ct = basic.construct_contingency(x[list(perm)], y[list(perm)])
        res.append(stats.fisher_exact(ct))

    return res


def line_plot_pvalues_concordance(
        pvals,
        concords,
        cmap=plt.get_cmap('Reds'),
        vmin=None,
        vmax=None,
):
    """
    Plot summarising pvalue and concordance level simultaneously using position (concordance) and colour (pval).
    Hard coded into vertical orientation (TODO: make this a parameter??)
    :param pvals: DataFrame. Index and columns will be used for ticklabels. Suggest using -log10
    :param concords: DataFrame, must match pvals
    :param pct_to_size_func:
    :param cmap: cmap used to represent
    :param vmin: For pvalue shading. If not supplied, min of data will be used
    :param vmax: For pvalue shading. If not supplied, max of data will be used
    :return:
    """
    if sorted(pvals.index) != sorted(concords.index):
        raise AttributeError("Index of pvals and concords must match")
    if sorted(pvals.columns) != sorted(concords.columns):
        raise AttributeError("Columns of pvals and concords must match")

    concords = concords.loc[pvals.index, pvals.columns]

    if vmin is None:
        vmin = pvals.values.min()
    if vmax is None:
        vmax = pvals.values.max()

    ny, nx = pvals.shape
    markers = common.get_best_marker_map(nx)

    gs = plt.GridSpec(nrows=nx, ncols=1, height_ratios=[1, 19])
    fig = plt.figure(figsize=(1.5 * nx, .5 * ny))
    ax = fig.add_subplot(gs[1])
    ax.invert_yaxis()

    cax = fig.add_subplot(gs[0])

    plogp_norm = colors.Normalize(vmin=vmin, vmax=vmax)
    plogp_sm = plt.cm.ScalarMappable(cmap=cmap, norm=plogp_norm)

    for i, col in enumerate(concords.columns):
        ax.scatter(
            concords[col],
            range(ny),
            c=[plogp_sm.to_rgba(t) for t in pvals[col].values],
            s=60,
            edgecolor='k',
            marker=markers[i]
        )

    ax.set_xlabel('% concordance', fontsize=12)
    ax.set_yticks(range(ny))
    ax.set_yticklabels(pvals.index, fontsize=12)
    ax.set_xlim([50, 100])
    plogp_sm.set_array(pvals)
    fig.colorbar(plogp_sm, cax=cax, orientation='horizontal')
    cax.xaxis.set_label_position('top')
    cax.set_xlabel(r'$-\log_{10}(p)$')

    type_attrs = {
        'class': 'line',
        'linestyle': 'none',
        'markeredgecolor': 'k',
        'markeredgewidth': 1.,
        'markerfacecolor': 'none',
        'markersize': 8
    }

    leg_dict = {}
    for i, col in enumerate(pvals.columns):
        leg_dict[col] = dict(type_attrs)
        leg_dict[col].update({'marker': markers[i]})

    common.add_custom_legend(ax, leg_dict, loc_outside=True, loc_outside_horiz='right')
    gs.update(left=0.2, bottom=0.1, right=0.72, top=0.95, hspace=0.12)

    return {
        'fig': fig,
        'ax': ax,
        'gs': gs,
        'cax': cax
    }


def de_dmr_concordance_scatter(de, dm, groups, de_edge=None, dm_edge=None,
                               group_colours=consts.METHYLATION_DIRECTION_COLOURS):
    """

    :param de:
    :param dm:
    :param groups:
    :param de_edge, dm_edge: If supplied, both must be specified. These provide a subset of the data in de and dm
    that will be highlighted on the plots with a black outline. Typically used to highlight DMR / TSS overlap.
    :param group_colours:
    :return:
    """
    a = de_edge is None
    b = dm_edge is None
    if a != b:
        raise AttributeError("Must specify either both of de_edge and dm_edge or neither.")

    ncol = max([len(t) for t in groups.values()])
    nrow = len(groups)

    fig, axs = plt.subplots(nrows=nrow, ncols=ncol, sharex=True, sharey=True, figsize=(1.14 * ncol, 1.8 * nrow))
    big_ax = common.add_big_ax_to_subplot_fig(fig)

    xmin = xmax = ymin = ymax = 0

    subplots_filled = set()
    for i, grp in enumerate(groups):
        for j, pid in enumerate(groups[grp]):
            subplots_filled.add((i, j))
            ax = axs[i, j]
            x = de[pid]
            y = dm[pid]
            xmin = min(xmin, x.min())
            xmax = max(xmax, x.max())
            ymin = min(ymin, y.min())
            ymax = max(ymax, y.max())

            ax.scatter(
                x,
                y,
                c=group_colours[grp.lower()],
                alpha=0.6
            )
            if not a:
                xm = de_edge[pid]
                ym = dm_edge[pid]
                ax.scatter(
                    xm,
                    ym,
                    c='none',
                    edgecolor='k',
                    linewidth=0.5,
                    alpha=0.6
                )
            ax.axhline(0., c='k', linestyle='--')
            ax.axvline(0., c='k', linestyle='--')
            ax.set_title(pid)
    big_ax.set_xlabel(r'DE logFC')
    big_ax.set_ylabel(r'DMR median $\Delta M$')
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.15, top=0.93, hspace=0.2, wspace=0.05)

    for i in range(2):
        for j in range(6):
            if (i, j) not in subplots_filled:
                axs[i, j].set_visible(False)
    axs[0, 0].set_xlim([np.floor(xmin * 1.05), np.ceil(xmax * 1.05)])
    axs[0, 0].set_ylim([np.floor(ymin), np.ceil(ymax)])

    return fig, axs


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

    ## Identify group-specific DMRs

    groups = {
        'Hypo': ['019', '030', '031', '017'],
        'Hyper': ['018', '050', '054', '061', '026', '052']
    }
    group_ind = setops.groups_to_ind(pids, groups)
    groups_inv = dictionary.complement_dictionary_of_iterables(groups, squeeze=True)

    # UpSet plot of joint DE/DMRs relative to groupings
    subgroup_set_colours = {
        'Hypo full': '#427425',
        'Hyper full': '#b31500',
        'Hypo partial': '#c5e7b1',
        'Hyper partial': '#ffa599',
        'Expanded core': '#4C72B0',
        'Specific': '#f4e842',
    }

    # start with an UpSet to represent this
    # all DE/DM relations
    de_dmr_by_member_all = [joint_de_dmr_s1[pid].index for pid in pids]

    upset = venn.upset_plot_with_groups(
        de_dmr_by_member_all,
        pids,
        group_ind,
        subgroup_set_colours,
        # min_size=10,
        n_plot=30,
        default_colour='gray'
    )
    upset['axes']['main'].set_ylabel('Number of DE/DMRs in set')
    upset['axes']['set_size'].set_xlabel('Number of DE/DMRs in single comparison')
    upset['figure'].savefig(os.path.join(outdir, "upset_de_dmr_by_direction_group_all.png"), dpi=200)

    # only TSS DE/DM relations
    de_dmr_by_member_tss = []
    for pid in pids:
        ix = joint_de_dmr_s1[pid][['dmr_%s' % t for t in relations_tss]].sum(axis=1).astype(bool)
        de_dmr_by_member_tss.append(
            joint_de_dmr_s1[pid].index[ix.values]
        )

    upset = venn.upset_plot_with_groups(
        de_dmr_by_member_tss,
        pids,
        group_ind,
        subgroup_set_colours,
        # min_size=10,
        n_plot=30,
        default_colour='gray'
    )
    upset['axes']['main'].set_ylabel('Number of DE/DMRs in set')
    upset['axes']['set_size'].set_xlabel('Number of DE/DMRs in single comparison')
    upset['figure'].savefig(os.path.join(outdir, "upset_de_dmr_by_direction_group_tss.png"), dpi=200)


    set_colours_dict = {
        'Hypo': 'g',
        'Hyper': 'r',
        'Discordant': 'b'
    }

    # Don't think we need this, but may be useful for a comparison?
    if False:
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
    de_dmr_dmr_median_delta_all = tmp['dmr_median_delta_m']
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
    de_dmr_dmr_median_delta_tss = tmp['dmr_median_delta_m']
    de_dmr_ipa_res_tss = export_de_dmr_groups_for_ipa(
        de_dmr_de_fdr_tss,
        de_dmr_de_logfc_tss,
        groups,
        os.path.join(outdir, "group_specific_de_for_ipa_tss.xlsx")
    )

    # look at the DMR direction dist (only) in the joint results - is it the same as when we consider:
    # - DMR only?
    # - DE only

    for grp in groups:
        # generate bar chart showing number / pct in each direction (DM - all)
        this_de_dmrs = de_dmrs_all[grp]
        this_dmrs = sorted(set([t[0] for t in this_de_dmrs]))

        # bar chart showing DMR direction
        this_results = {}
        for pid in groups[grp]:
            this_results[pid] = dict([(k, {'median_change': t}) for k, t in de_dmr_dmr_median_delta_all[grp][pid].dropna().iteritems()])

        plt_dict = addd.dm_direction_bar_plot(this_results, groups[grp], figsize=(len(groups[grp]) - .5, 4.5))
        plt_dict['fig'].savefig(os.path.join(outdir, "dmr_direction_by_group_%s_all.png" % grp.lower()), dpi=200)

        # generate bar chart showing number / pct in each direction (DM - TSS)
        this_de_dmrs = de_dmrs_tss[grp]
        this_dmrs = sorted(set([t[0] for t in this_de_dmrs]))

        # bar chart showing DMR direction
        this_results = {}
        for pid in groups[grp]:
            this_results[pid] = dict([(k, {'median_change': t}) for k, t in de_dmr_dmr_median_delta_tss[grp][pid].dropna().iteritems()])

        plt_dict = addd.dm_direction_bar_plot(this_results, groups[grp], figsize=(len(groups[grp]) - .5, 4.5))
        plt_dict['fig'].savefig(os.path.join(outdir, "dmr_direction_by_group_%s_tss.png" % grp.lower()), dpi=200)

        # bar chart showing DE direction
        for_plot = {}
        for pid in groups[grp]:
            for_plot[pid] = de_dmr_de_logfc_all[grp][[pid]].dropna()
            for_plot[pid].columns = ['logFC']
        plt_dict = same_de.bar_plot(for_plot, keys=groups[grp], figsize=(len(groups[grp]) - .5, 4.5))
        plt_dict['fig'].savefig(os.path.join(outdir, "de_direction_by_group_%s_all.png" % grp.lower()), dpi=200)

        for_plot = {}
        for pid in groups[grp]:
            for_plot[pid] = de_dmr_de_logfc_tss[grp][[pid]].dropna()
            for_plot[pid].columns = ['logFC']
        plt_dict = same_de.bar_plot(for_plot, keys=groups[grp], figsize=(len(groups[grp]) - .5, 4.5))
        plt_dict['fig'].savefig(os.path.join(outdir, "de_direction_by_group_%s_tss.png" % grp.lower()), dpi=200)

    # Venn diagrams of DE
    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    plot_venn_de_directions(de_dmr_de_logfc_all, set_colours_dict, ax=ax)
    fig.savefig(os.path.join(outdir, "de_from_group_spec_de_dmrs_all.png"), dpi=200)

    fig = plt.figure(figsize=(5., 3.3))
    ax = fig.add_subplot(111)
    plot_venn_de_directions(de_dmr_de_logfc_tss, set_colours_dict, ax=ax)
    fig.savefig(os.path.join(outdir, "de_from_group_spec_de_dmrs_tss.png"), dpi=200)

    # assess concordance between DM and DE direction
    # start with scatterplot

    de_data_all = {}
    dm_data_all = {}
    de_data_tss = {}
    dm_data_tss = {}

    ct_all = {}
    ct_tss = {}
    fisher_p = {}
    fisher_p_tss = {}

    for pid in pids:
        grp = groups_inv[pid]
        this_ix = de_dmrs_all[grp]
        this_ix_tss = de_dmrs_tss[grp]

        this_y = de_dmr_dmr_median_delta_all[grp].loc[[t[0] for t in this_ix]]
        this_x = de_dmr_de_logfc_all[grp].reindex([t[1] for t in this_ix])

        this_y_tss = de_dmr_dmr_median_delta_tss[grp].loc[[t[0] for t in this_ix_tss]]
        this_x_tss = de_dmr_de_logfc_tss[grp].reindex([t[1] for t in this_ix_tss])

        # why are there missing lookups? Are those genes we decided were duplicated??
        # if so: can we just remove these in the original function?
        missing_in_x = this_x.isnull().all(axis=1)
        this_x = this_x.loc[~missing_in_x.values]
        this_y = this_y.loc[~missing_in_x.values]

        missing_in_x_tss = this_x_tss.isnull().all(axis=1)
        this_x_tss = this_x_tss.loc[~missing_in_x_tss.values]
        this_y_tss = this_y_tss.loc[~missing_in_x_tss.values]

        ix_x = ~this_x[pid].isnull()
        ix_y = ~this_y[pid].isnull()

        ix = ix_x.values & ix_y.values

        de_data_all[pid] = this_x[pid].loc[ix]
        dm_data_all[pid] = this_y[pid].loc[ix]

        ct_all[pid] = basic.construct_contingency(this_x[pid].loc[ix].values, this_y[pid].loc[ix].values)
        fisher_p[pid] = stats.fisher_exact(ct_all[pid])

        # ensure we drop corresponding rows in x and y
        ix_x_tss = ~this_x_tss[pid].isnull()
        ix_y_tss = ~this_y_tss[pid].isnull()

        ix_tss = ix_x_tss.values & ix_y_tss.values

        de_data_tss[pid] = this_x_tss[pid].loc[ix_tss]
        dm_data_tss[pid] = this_y_tss[pid].loc[ix_tss]

        ct_tss[pid] = basic.construct_contingency(this_x_tss[pid].loc[ix_tss].values, this_y_tss[pid].loc[ix_tss].values)
        fisher_p_tss[pid] = stats.fisher_exact(ct_tss[pid])

    fig, axs = de_dmr_concordance_scatter(
        de_data_all,
        dm_data_all,
        groups,
        de_edge=de_data_tss,
        dm_edge=dm_data_tss,
    )

    fig.savefig(os.path.join(outdir, "gs_de_dmr_scatter.png"), dpi=200)

    pvals = pd.DataFrame(index=pids, columns=['All', 'TSS'])
    pvals.loc[:, 'All'] = [-np.log10(fisher_p[pid][1]) for pid in pids]
    pvals.loc[:, 'TSS'] = [-np.log10(fisher_p_tss[pid][1]) for pid in pids]

    concords = pd.DataFrame(index=pids, columns=['All', 'TSS'])
    concords.loc[:, 'All'] = pd.Series(dict([
        (pid, (ct_all[pid][1, 0] + ct_all[pid][0, 1]) / float(ct_all[pid].sum()) * 100.) for pid in pids
    ]))
    concords.loc[:, 'TSS'] = pd.Series(dict([
        (pid, (ct_tss[pid][1, 0] + ct_tss[pid][0, 1]) / float(ct_tss[pid].sum()) * 100.) for pid in pids
    ]))

    plt_dict = line_plot_pvalues_concordance(pvals, concords)
    plt_dict['fig'].savefig(os.path.join(outdir, "fisher_p_concordance_line_full_dmr_de_genes.png"), dpi=200)

    raise StopIteration

    # ncol = max([len(t) for t in groups.values()])
    # fig, axs = plt.subplots(
    #     nrows=len(groups),
    #     ncols=ncol,
    #     sharex=True,
    #     sharey=True,
    #     figsize=(1.4 * ncol, 1.8 * len(groups))
    # )
    # big_ax = common.add_big_ax_to_subplot_fig(fig)
    #
    # subplots_filled = set()
    # for i, grp in enumerate(groups):
    #     # combination of DE/DMRs that are group-specific
    #     this_ix = de_dmrs_all[grp]
    #     this_ix_tss = de_dmrs_tss[grp]
    #
    #     this_y = de_dmr_dmr_median_delta_all[grp].loc[[t[0] for t in this_ix]]
    #     this_x = de_dmr_de_logfc_all[grp].reindex([t[1] for t in this_ix])
    #
    #     this_y_tss = de_dmr_dmr_median_delta_tss[grp].loc[[t[0] for t in this_ix_tss]]
    #     this_x_tss = de_dmr_de_logfc_tss[grp].reindex([t[1] for t in this_ix_tss])
    #
    #     # why are there missing lookups? Are those genes we decided were duplicated??
    #     # if so: can we just remove these in the original function?
    #     missing_in_x = this_x.isnull().all(axis=1)
    #     this_x = this_x.loc[~missing_in_x.values]
    #     this_y = this_y.loc[~missing_in_x.values]
    #
    #     missing_in_x_tss = this_x_tss.isnull().all(axis=1)
    #     this_x_tss = this_x_tss.loc[~missing_in_x_tss.values]
    #     this_y_tss = this_y_tss.loc[~missing_in_x_tss.values]
    #
    #     for j, pid in enumerate(groups[grp]):
    #         # ensure we drop corresponding rows in x and y
    #         ix_x = ~this_x[pid].isnull()
    #         ix_y = ~this_y[pid].isnull()
    #
    #         ix = ix_x.values & ix_y.values
    #
    #         de_data_all[pid] = this_x[pid].loc[ix]
    #         dm_data_all[pid] = this_y[pid].loc[ix]
    #
    #         subplots_filled.add((i, j))
    #         ax = axs[i, j]
    #
    #         ax.scatter(
    #             this_x[pid].loc[ix],
    #             this_y[pid].loc[ix],
    #             c=consts.METHYLATION_DIRECTION_COLOURS[grp.lower()],
    #             alpha=0.5
    #         )
    #         # fisher test
    #         ct = basic.construct_contingency(this_x[pid].loc[ix].values, this_y[pid].loc[ix].values)
    #         fisher_p[pid] = stats.fisher_exact(ct)
    #
    #         # ensure we drop corresponding rows in x and y
    #         ix_x_tss = ~this_x_tss[pid].isnull()
    #         ix_y_tss = ~this_y_tss[pid].isnull()
    #
    #         ix_tss = ix_x_tss.values & ix_y_tss.values
    #
    #         de_data_tss[pid] = this_x_tss[pid].loc[ix_tss]
    #         dm_data_tss[pid] = this_y_tss[pid].loc[ix_tss]
    #
    #         # highlight TSS-linked (only)
    #         ax.scatter(
    #             this_x_tss[pid].loc[ix_tss],
    #             this_y_tss[pid].loc[ix_tss],
    #             facecolors='none',
    #             edgecolor='k',
    #             linewidth=1.,
    #             alpha=0.5
    #         )
    #         # fisher test
    #         ct_tss = basic.construct_contingency(this_x_tss[pid].loc[ix_tss].values, this_y_tss[pid].loc[ix_tss].values)
    #         fisher_p_tss[pid] = stats.fisher_exact(ct_tss)
    #
    #         ax.axhline(0., c='k', linestyle='--')
    #         ax.axvline(0., c='k', linestyle='--')
    #         ax.set_title(pid)
    # big_ax.set_xlabel(r'DE logFC')
    # big_ax.set_ylabel(r'DMR median $\Delta M$')
    # fig.tight_layout()
    #
    # for i in range(axs.shape[0]):
    #     for j in range(axs.shape[1]):
    #         if (i, j) not in subplots_filled:
    #             axs[i, j].set_visible(False)
    # axs[0,0].set_xlim([-18, 18])
    #
    # fig.subplots_adjust(left=0.07, bottom=0.15, top=0.92, right=0.98, wspace=0.04, hspace=0.18)
    #
    # fig.savefig(os.path.join(outdir, "gs_de_dmr_scatter.png"), dpi=200)

    # the issue with comparing between patients here is that the number of data points changes
    # therefore 'downsample' all to give a common number

    n_iter = 1000
    lowest_n_all = min([len(v) for v in de_data_all.values()])
    lowest_n_tss = min([len(v) for v in de_data_tss.values()])

    fisher_p_ds_all = {}
    fisher_p_ds_tss = {}

    for pid in pids:
        this_x = de_data_all[pid]
        this_x_tss = de_data_tss[pid]
        this_y = dm_data_all[pid]
        this_y_tss = dm_data_tss[pid]

        try:
            fisher_p_ds_all[pid] = fishers_downsample(this_x, this_y, lowest_n_all, n_iter=n_iter)
            fisher_p_ds_tss[pid] = fishers_downsample(this_x_tss, this_y_tss, lowest_n_tss, n_iter=n_iter)
        except ValueError:
            # this will happen if we manually override the lowest_n
            pass

    # the `if` clause here is unnecessary unless we override lowest_n_{all,tss}, but does no harm
    df_all = pd.DataFrame([np.array(fisher_p_ds_all[pid])[:, 1] if pid in fisher_p_ds_all else [1.0] for pid in pids], index=pids).transpose()
    df_tss = pd.DataFrame([np.array(fisher_p_ds_tss[pid])[:, 1] if pid in fisher_p_ds_tss else [1.0] for pid in pids], index=pids).transpose()

    cols = [consts.METHYLATION_DIRECTION_COLOURS[groups_inv[pid].lower()] for pid in pids]
    plt_kws = dict(
        color=cols,
        edgecolor='k',
        linewidth=1.
    )
    plt_alpha = 0.01
    fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(5., 5.5))
    ((df_tss < plt_alpha).sum() / float(n_iter) * 100.).loc[pids].plot.bar(ax=axs[0], **plt_kws)
    ((df_all < plt_alpha).sum() / float(n_iter) * 100.).loc[pids].plot.bar(ax=axs[1], **plt_kws)
    axs[0].set_ylabel("All relations % significant")
    axs[1].set_ylabel("TSS relations % significant")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "permutation_fisher_pct_significant.png"), dpi=200)

    # quantify the % concordance
    pct_conc = pd.DataFrame(index=pids, columns=['All', 'TSS'])

    for pid in pids:
        this_x = de_data_all[pid].values
        this_x_tss = de_data_tss[pid].values
        this_y = dm_data_all[pid].values
        this_y_tss = dm_data_tss[pid].values

        this_joint = joint_de_dmr_s1[pid].loc[joint_de_dmr_s1[pid].cluster_id.isin(groups[grp])]
        pct_conc.loc[pid, 'All'] = (np.sign(this_x) != np.sign(this_y)).sum() / float(this_x.size) * 100.

        this_joint = this_joint.loc[(this_joint.dmr_TSS1500 | this_joint.dmr_TSS200).astype(bool)]
        pct_conc.loc[pid, 'TSS'] = (np.sign(this_x_tss) != np.sign(this_y_tss)).sum() / float(this_x_tss.size) * 100.

    # permutation test: do we see the same level of concordance with a 'random' background selection?
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

        # now look at the group-specific DMRs to determine how many we need to draw
        n_hypo = (dm_data_all[pid] < 0).sum()
        n_hyper = (dm_data_all[pid] > 0).sum()

        logger.info(
            "Patient %s. Drawing %d hypo and %d hyper DMRs from background (%d hypo / %d hyper).",
            pid,
            n_hypo,
            n_hyper,
            bg_hypo.shape[0],
            bg_hyper.shape[0]
        )

    #     this_perm_res = []
    #     for i in range(n_iter_bg):
    #
    #
    #         # multiple random draws
    #         pass
