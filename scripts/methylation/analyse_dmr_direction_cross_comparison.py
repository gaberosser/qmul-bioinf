from plotting import bar, common, pie
from methylation import loader, dmr, process
import pandas as pd
from stats import nht
from utils import output, setops, genomics, log
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
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts
from scripts.methylation.analyse_dmr_direction_validation_cohort import run_dmr_analyses, bar_plot
from scripts.methylation import analyse_dmr_direction_and_distribution as addd

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
logger = log.get_console_logger()

"""
Here we analyse the direction (and genomic locations?) of DMRs in a cross-comparison.
We compare our GIC lines with (non-matching) iNSC lines.
"""
def get_dmr_locations(
    dmr_res,
    clusters,
    coord_summary_method='first'
):
    """

    :param dmr_res:
    :param clusters:
    :param chrom_lengths:
    :param window_size:
    :param split_by_direction:
    :param coord_summary_method: Method used to reduce the list of CpG coordinates to a single one.
    Default is 'first', meaning take the 5'-most coordinate (first in the list). Other options: 'last', 'median', 'mean'
    :return:
    """
    dmr_loci_hypo = {}
    dmr_loci_hyper = {}

    for pid in dmr_res:
        this_loci_hyper = collections.defaultdict(dict)
        this_loci_hypo = collections.defaultdict(dict)

        for cluster_id, cl in dmr_res[pid].items():
            # get the chrom and locus
            pc = clusters[cluster_id]
            # use the requested method to get a representative coordinate from the list
            if coord_summary_method == 'first':
                the_coord = pc.coord_list[0]
            elif coord_summary_method == 'last':
                the_coord = pc.coord_list[-1]
            elif coord_summary_method == 'median':
                the_coord = np.median(pc.coord_list)
            elif coord_summary_method == 'mean':
                the_coord = np.mean(pc.coord_list)
            else:
                raise ValueError("Unsupported coordinate summary method '%s'." % coord_summary_method)



            if cl['median_change'] > 0:
                this_loci_hyper[pc.chr][the_coord] = cl['median_change']
            else:
                this_loci_hypo[pc.chr][the_coord] = cl['median_change']

        dmr_loci_hyper[pid] = this_loci_hyper
        dmr_loci_hypo[pid] = this_loci_hypo

    return {
        'hyper': dmr_loci_hyper,
        'hypo': dmr_loci_hypo,
    }


def pairwise_ks_test(dat):
    """
    Run the 2 sample KS test on every possible pair of data
    :param dat: DataFrame containing the inputs. Columns represent samples. Rows represent the features.
    The values are integers, representing the number of times each feature is observed.
    :return:
    """
    res = pd.DataFrame(
        index=pd.Index(dat.columns, name='cmp1'),
        columns=pd.Index(dat.columns, name='cmp2'),
        dtype=float
    )
    for a in res.index:
        xa = dat[a]
        # tile the data for KS test
        xa = xa.index.repeat(xa.values)
        for b in res.columns:
            xb = dat[b]
            # tile the data for KS test
            xb = xb.index.repeat(xb.values)
            res.loc[a, b] = stats.ks_2samp(xa, xb)[1]
    return res


def mht_correction_to_pairwise_pvals(df):
    """
    Apply a MHT correction to the symmetric pairwise distance matrix df. NB we assume symmetry and only use the upper
    triangular values
    :param df:
    :return:
    """
    ud_ix = np.triu_indices_from(df, k=1)
    ud = df.values[ud_ix]
    ud_mht = nht.mht_p_adjust(ud, method='BH')
    res = pd.DataFrame(
        index=df.index,
        columns=df.columns,
        dtype=float
    )
    res.values[ud_ix] = ud_mht
    return res


def log_pval_pcolor(df, figsize=None, pk_alpha=-np.log10(0.05)):
    to_plot = -np.log10(df)
    to_plot_masked = np.ma.masked_where(np.isnan(to_plot), to_plot)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    h1 = ax.pcolor(
        to_plot_masked,
        vmin=pk_alpha,
        cmap='Reds',
        edgecolor='w',
        linewidth=1.,
        zorder=5,
    )

    to_plot_masked = np.ma.masked_where(np.isnan(to_plot) | (to_plot < pk_alpha), to_plot)
    h2 = ax.pcolor(
        to_plot_masked,
        facecolor='none',
        edgecolor='k',
        linewidth=1.,
        zorder=10,
    )
    # hack - not sure why this is necessary
    h2.set_edgecolor('k')
    ax.invert_yaxis()
    ax.yaxis.tick_right()
    ax.xaxis.tick_top()
    ax.xaxis.set_ticks(np.arange(0, len(pids) + 1) + 0.5)
    ax.xaxis.set_ticklabels(df.columns, rotation=90)
    ax.yaxis.set_ticks(np.arange(0, len(pids) + 1) + 0.5)
    ax.yaxis.set_ticklabels(df.index, rotation=0)
    cbar = fig.colorbar(h1, pad=0.1, shrink=0.8)
    cbar.set_label(r'$-\log_{10}(\mathrm{FDR})$')

    fig.tight_layout()

    return {
        'fig': fig,
        'ax': ax,
        'cbar': cbar,
        'h1': h1,
        'h2': h2
    }


if __name__ == "__main__":
    pids = consts.PIDS
    norm_method = 'swan'
    alpha = 0.05
    pk_alpha = -np.log10(alpha)


    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    # load our data
    our_obj = loader.load_by_patient(pids, norm_method=norm_method, samples=consts.S1_METHYL_SAMPLES)
    anno = loader.load_illumina_methylationepic_annotation()
    our_obj.meta.insert(0, 'patient_id', our_obj.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    dat = process.m_from_beta(our_obj.data)
    meta = our_obj.meta
    common_probes = anno.index.intersection(dat.index)
    dat = dat.reindex(common_probes)
    anno = anno.reindex(common_probes)

    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method

    the_hash = tsgd.dmr_results_hash(meta.sort_index().index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_cross_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        dmr_res = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        logger.info("Unable to locate pre-existing results. Computing from scratch (this can take a while).")
        comparisons = {}
        for pid in pids:  # GIC
            for pid2 in pids:  # iNSC
                gics = meta.index[(meta.type == 'GBM') & (meta.patient_id == pid)]
                inscs = meta.index[(meta.type == 'iNSC') & (meta.patient_id == pid2)]
                comparisons['-'.join([pid, pid2])] = [gics, inscs]
        dmr_res = run_dmr_analyses(dat, comparisons, anno, dmr_params)
        # Save DMR results to disk
        dmr_res.to_pickle(fn, include_annotation=False)
        logger.info("Saved DMR results to %s", fn)

    dmr_res_all = dmr_res.results_significant

    # number of DMRs, split by direction

    n_dmr = pd.DataFrame(index=pids, columns=pids)
    n_dmr.index.name = 'GIC'
    n_dmr.columns.name = 'iNSC'
    n_hyper = n_dmr.copy()
    n_hypo = n_dmr.copy()

    for k, v in dmr_res_all.items():
        p1, p2 = k.split('-')
        n_dmr.loc[p1, p2] = len(v)
        n_hyper.loc[p1, p2] = len([t for t in v.values() if t['median_change'] > 0])
        n_hypo.loc[p1, p2] = len([t for t in v.values() if t['median_change'] < 0])

    n_hyper_pct= n_hyper / n_dmr * 100.

    # box whisker plot

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = sns.boxplot(n_hyper_pct.transpose(), whis=5, ax=ax)
    ax.scatter(range(len(pids)), np.diagonal(n_hyper_pct), facecolor='k')
    ax.set_ylim([0, 100])

    # run down the rows or columns and generate an 'overlap spectrum' for each one
    # rows: check the effect of varying the iNSC line (CONSISTENCY)
    # cols: check the effect of varying the GIC line (non-syngeneic DIFFERENCE)
    # also repeat for the columns, which is just the S1 approach (SYNGENEIC)

    row_collapse = pd.DataFrame(
        dict([
            (
                pid,
                setops.quantify_feature_membership(
                    setops.venn_from_arrays(
                        *[dmr_res_all['%s-%s' % (pid, p)].keys() for p in pids]
                    )[1]
                )
            )
            for pid in pids
        ])
    )[pids]

    col_collapse = pd.DataFrame(
        dict([
            (
                pid,
                setops.quantify_feature_membership(
                    setops.venn_from_arrays(
                        *[dmr_res_all['%s-%s' % (p, pid)].keys() for p in pids]
                    )[1]
                )
            )
            for pid in pids
        ])
    )[pids]

    syn_dist = setops.quantify_feature_membership(
        setops.venn_from_arrays(
            *[dmr_res_all['%s-%s' % (p, p)].keys() for p in pids]
        )[1]
    )

    # bar charts
    fig, axs = plt.subplots(len(pids), sharex=True, sharey=True, figsize=(3.5, 7.5))
    big_ax = fig.add_subplot(111, frameon=False)
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.grid(False)
    for i, pid in enumerate(pids):
        axs[i].bar(range(1, len(pids) + 1), row_collapse[pid])
        axs[i].set_ylabel(pid)
    big_ax.set_ylabel('Number in set')
    fig.subplots_adjust(left=0.25, right=0.98, bottom=0.07, top=0.98)
    big_ax.set_position([0.15, 0, 1, 1])
    axs[-1].xaxis.set_ticks(range(1, len(pids) + 1))
    axs[-1].set_xlabel('Number of intersecting comparisons')
    fig.savefig(os.path.join(outdir, "row_consistency_intersection_freq.png"), dpi=200)

    fig, axs = plt.subplots(len(pids) + 1, sharex=True, sharey=True, figsize=(3.5, 7.5))
    big_ax = fig.add_subplot(111, frameon=False)
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.grid(False)
    axs[0].bar(range(1, len(pids) + 1), syn_dist, color='k')
    axs[0].set_ylabel('Syngeneic')
    for i, pid in enumerate(pids):
        axs[i + 1].bar(range(1, len(pids) + 1), col_collapse[pid], edgecolor='k')
        axs[i + 1].set_ylabel(pid)
    big_ax.set_ylabel('Number in set')
    fig.subplots_adjust(left=0.25, right=0.98, bottom=0.07, top=0.98)
    big_ax.set_position([0.15, 0, 1, 1])
    axs[-1].xaxis.set_ticks(range(1, len(pids) + 1))
    axs[-1].set_xlabel('Number of intersecting comparisons')
    fig.savefig(os.path.join(outdir, "col_nonsyn_intersection_freq.png"), dpi=200)

    # pairwise KS test for the column collapsed data
    for_testing = col_collapse.copy()
    for_testing.insert(0, 'syngeneic', syn_dist)
    ks_2samp = pairwise_ks_test(for_testing)
    ks_2samp_mht = mht_correction_to_pairwise_pvals(ks_2samp)

    # plot heatmap of logged pvals
    plt_dict = log_pval_pcolor(ks_2samp_mht)
    ## FIXME: this only works when run interactively?! Something to do with the choice of backend?
    plt_dict['h2'].set_edgecolor('k')
    plt_dict['fig'].savefig(os.path.join(outdir, "padj_syn_nonsyn.png"), dpi=200)

    # we will later show that 026 is less consistent than the others
    # let's repeat but without it?
    # plot heatmap of logged pvals
    plt_dict = log_pval_pcolor(ks_2samp_mht.drop('026', axis=0).drop('026', axis=1))
    ## FIXME: this only works when run interactively?! Something to do with the choice of backend?
    plt_dict['h2'].set_edgecolor('k')
    plt_dict['fig'].savefig(os.path.join(outdir, "padj_syn_nonsyn_no026.png"), dpi=200)

    # repeat but with row collapsed samples ('consistency')
    ks_2samp_row = pairwise_ks_test(row_collapse)
    ks_2samp_row_mht = mht_correction_to_pairwise_pvals(ks_2samp_row)

    # plot heatmap of logged pvals
    plt_dict = log_pval_pcolor(ks_2samp_row_mht)
    ## FIXME: this only works when run interactively?! Something to do with the choice of backend?
    plt_dict['h2'].set_edgecolor('k')
    plt_dict['fig'].savefig(os.path.join(outdir, "padj_row_consistency.png"), dpi=200)

    # repeat without 026
    plt_dict = log_pval_pcolor(ks_2samp_row_mht.drop('026', axis=0).drop('026', axis=1))
    ## FIXME: this only works when run interactively?! Something to do with the choice of backend?
    plt_dict['h2'].set_edgecolor('k')
    plt_dict['fig'].savefig(os.path.join(outdir, "padj_row_consistency_no026.png"), dpi=200)

    # clustermap of columns
    cc = col_collapse.copy()
    cc.insert(0, 'syngeneic', syn_dist)
    cg = sns.clustermap(
        cc,
        row_cluster=False,
        cmap='Reds',
        metric='correlation',
        z_score=None,
        method='average',
        figsize=(6.7, 7.7)
    )
    fig = cg.ax_heatmap.figure
    cg.gs.update(left=0.05, top=0.98, bottom=0.05, right=0.95, wspace=-.1)
    fig.savefig(os.path.join(outdir, "syn_nonsyn_clustermap.png"), dpi=200)

    # clustermap of rows
    # here we see that 026 is an outlier of sorts (at least for this analysis?)
    cg = sns.clustermap(
        row_collapse,
        row_cluster=False,
        cmap='Reds',
        metric='correlation',
        z_score=None,
        method='average',
        figsize=(6.7, 7.7)
    )
    fig = cg.ax_heatmap.figure
    cg.gs.update(left=0.05, top=0.98, bottom=0.05, right=0.95, wspace=-.1)
    fig.savefig(os.path.join(outdir, "row_consistency_clustermap.png"), dpi=200)

    # clustermap of rows
    # remove the outlier 026
    cg = sns.clustermap(
        row_collapse.drop('026', axis=1),
        row_cluster=False,
        cmap='Reds',
        metric='correlation',
        z_score=None,
        method='average',
        figsize=(6.7, 7.7)
    )
    fig = cg.ax_heatmap.figure
    cg.gs.update(left=0.05, top=0.98, bottom=0.05, right=0.95, wspace=-.1)
    fig.savefig(os.path.join(outdir, "row_consistency_clustermap_no026.png"), dpi=200)

    ## TODO: we really need some kind of estimate of the 'normal' amount of variation in these categories
    # this will provide a reference frame for the observed differences.
    # The issue with this idea is that any simulation of the sets results in dramatic departures from the observed
    # numbers of DMRs in each region, so we can't compare it directly?

    # focus on one patient: are we seeing the same loci in syngeneic and non??
    # start with a single GIC line...
    pid = '018'
    # ...and iterate over the possible comparators

    # first: we need chromosome lengths
    chroms = [str(t) for t in range(1, 23)]
    fa_fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human/ensembl/GRCh37/fa/Homo_sapiens.GRCh37.dna.primary_assembly.fa'
    )
    window_size = int(2e5)
    coord_summary_method = 'first'
    chrom_length = genomics.feature_lengths_from_fasta(fa_fn, features=chroms)

    # get probe density
    bg_density = {}
    for ch in chroms:
        this_coord = anno.loc[anno.CHR == ch].MAPINFO
        tt = np.histogram(anno.loc[anno.CHR == ch].MAPINFO, np.arange(0, chrom_length[ch], window_size))
        bg_density[ch] = pd.Series(tt[0], index=tt[1][:-1])

    # extract the DMRs for this GIC line
    this_res = dict([
        (p, dmr_res_all['-'.join([pid, p])]) for p in pids
    ])
    tmp = get_dmr_locations(
        this_res,
        dmr_res.clusters,
    )

    dmr_loci_hyper = tmp['hyper']
    dmr_loci_hypo = tmp['hypo']

    for p in pids:
        hyper = dmr_loci_hyper[p]
        hypo = dmr_loci_hypo[p]
        hyper_for_plot = dict([(ch, sorted(hyper[ch].keys())) for ch in chroms])
        hypo_for_plot = dict([(ch, sorted(hypo[ch].keys())) for ch in chroms])
        hyper_hist = dict([
            (
                ch,
                np.histogram(
                    hyper_for_plot[ch],
                    range(0, chrom_length[ch], window_size) + [chrom_length[ch]]
                )[0]
            ) for ch in chroms
        ])
        hypo_hist = dict([
            (
                ch,
                np.histogram(
                    hypo_for_plot[ch],
                    range(0, chrom_length[ch], window_size) + [chrom_length[ch]]
                )[0]
            ) for ch in chroms
        ])

    plt_dict = addd.polar_distribution_plot_with_kdes(
        hyper_for_plot,
        hypo_for_plot,
        chrom_length,
        bg_density=bg_density,
        window_size=window_size
    )


