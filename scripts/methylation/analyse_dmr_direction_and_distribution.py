from plotting import bar, common
from methylation import loader, dmr, process
import pandas as pd
from statsmodels.sandbox.stats import multicomp
from utils import output, setops, genomics, log
import multiprocessing as mp
import os
import collections
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt, patches
from matplotlib.colors import Normalize
from sklearn.neighbors import KernelDensity

import seaborn as sns
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts
from plotting import common

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR
logger = log.get_console_logger()

"""
Here we analyse the direction and genomic locations of DMRs in the individual patients.
We note a bias in the direction, which is amplified considerably when we consider patient-specific DMRs.
"""

def get_binned_dmr_locations(
    dmr_res,
    clusters,
    chrom_lengths,
    window_size=20000,
    split_by_direction=False,
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
    if split_by_direction:
        dmr_loci_hypo = {}
        dmr_loci_binned_hypo = {}
        dmr_loci_hyper = {}
        dmr_loci_binned_hyper = {}
    else:
        dmr_loci = {}
        dmr_loci_binned = {}

    for pid in dmr_res:
        if split_by_direction:
            dmr_loci_binned_hypo[pid] = {}
            this_loci_hypo = collections.defaultdict(list)
            dmr_loci_binned_hyper[pid] = {}
            this_loci_hyper = collections.defaultdict(list)
        else:
            dmr_loci_binned[pid] = {}
            this_loci = collections.defaultdict(list)
        # this_loci = dict([(chrom, []) for chrom in chroms])
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

            if split_by_direction:
                if cl['median_change'] > 0:
                    this_loci_hyper[pc.chr].append(the_coord)
                else:
                    this_loci_hypo[pc.chr].append(the_coord)
            else:
                this_loci[pc.chr].append(the_coord)
        if split_by_direction:
            dmr_loci_hyper[pid] = this_loci_hyper
            dmr_loci_hypo[pid] = this_loci_hypo
        else:
            dmr_loci[pid] = this_loci
        # run the histogram process on each chrom
        if split_by_direction:
            for chrom, arr in this_loci_hyper.items():
                edges = range(1, chrom_lengths[chrom] + 1, window_size)
                if edges[-1] != chrom_lengths[chrom]:
                    edges.append(chrom_lengths[chrom])
                    this_counts, _ = np.histogram(arr, edges)
                    dmr_loci_binned_hyper[pid][chrom] = pd.Series(this_counts, index=edges[:-1])
            for chrom, arr in this_loci_hypo.items():
                edges = range(1, chrom_lengths[chrom] + 1, window_size)
                if edges[-1] != chrom_lengths[chrom]:
                    edges.append(chrom_lengths[chrom])
                    this_counts, _ = np.histogram(arr, edges)
                    dmr_loci_binned_hypo[pid][chrom] = pd.Series(this_counts, index=edges[:-1])
        else:
            for chrom, arr in this_loci.items():
                edges = range(1, chrom_lengths[chrom] + 1, window_size)
                if edges[-1] != chrom_lengths[chrom]:
                    edges.append(chrom_lengths[chrom])
                    this_counts, _ = np.histogram(arr, edges)
                    dmr_loci_binned[pid][chrom] = pd.Series(this_counts, index=edges[:-1])
    if split_by_direction:
        return {
            'dmr_loci_hyper': dmr_loci_hyper,
            'dmr_loci_hypo': dmr_loci_hypo,
            'dmr_binned_hyper': dmr_loci_binned_hyper,
            'dmr_binned_hypo': dmr_loci_binned_hypo,
        }
    else:
        return {
            'dmr_loci': dmr_loci,
            'dmr_binned': dmr_loci_binned
        }


def ks_test_dmr_locations(
        dmr_locs,
        null_locs,
        mht_method='fdr_bh',
):
    """
    Run the 2 sample KS test on DMR locations, using the supplied null locations as the reference sample.
    :param dmr_locs: Nested dictionary keyed by PID then feature (chromosome name)
    :param null_locs: Dictionary keyed by feature (chromosome name)
    :param mht_method: The method to use for MHT correction. Passed to multicomp.multipletests.
    If None, skip MHT and use the unadjusted pvalues.
    :return:
    """
    ks_pval_flat = []
    ks_stat_flat = []
    ks_flat_map = []
    for pid, val in dmr_locs.items():
        for feat in val:
            tmp = stats.ks_2samp(dmr_locs[pid][feat], null_locs[feat])
            ks_stat_flat.append(tmp[0])
            ks_pval_flat.append(tmp[1])
            ks_flat_map.append((pid, feat))

    # correct for MHT
    if mht_method is None:
        padj = ks_pval_flat
    else:
        _, padj, _, _ = multicomp.multipletests(ks_pval_flat, alpha=0.05, method=mht_method)

    res_fdr = dict([(pid, {}) for pid in dmr_locs.keys()])
    res_stat = dict([(pid, {}) for pid in dmr_locs.keys()])

    i = 0
    for pid, feat in ks_flat_map:
        res_fdr[pid][feat] = padj[i]
        res_stat[pid][feat] = ks_stat_flat[i]
        i += 1

    return res_fdr, res_stat


def fit_kde_dmr_location(locs, xi, bandwidth, normed=True):
    # the KDE library expects a 2D array, even if data are not multidimensional
    locs = np.array(locs)
    if len(locs.shape) == 1:
        locs = locs[:, None]
    elif len(locs.shape) == 2:
        if locs.shape[0] == 1 and locs.shape[1] > 1:
            locs = locs.transpose()
        elif locs.shape[1] == 1 and locs.shape[0] > 1:
            pass
        else:
            raise ValueError("locations array must be 1D or (quasi-)2D")
    else:
        raise ValueError("locations array must be 1D or (quasi-)2D")

    xi = np.array(xi)
    if len(xi.shape) == 1:
        xi = xi[:, None]
    elif len(xi.shape) == 2:
        if xi.shape[0] == 1 and xi.shape[1] > 1:
            xi = xi.transpose()
        elif xi.shape[1] == 1 and xi.shape[0] > 1:
            pass
        else:
            raise ValueError("xi array must be 1D or (quasi-)2D")
    else:
        raise ValueError("xi array must be 1D or (quasi-)2D")

    k = KernelDensity(bandwidth=bandwidth, kernel='gaussian')
    k.fit(locs)  # this increases the dim of the 1D array
    logd_score = k.score_samples(xi)
    # blank out baseline (for plotting purposes)
    logd_score[logd_score < -100] = -np.inf
    score = np.ma.masked_equal(np.exp(logd_score), 0.)
    if not normed:
        # convert to unnormalised density
        score *= locs.size
    return score


def polar_distribution_plot(
    dmr_loci_hyper,
    dmr_loci_hypo,
    chrom_length,
    bg_density=None,
    unmapped_density=None,
    window_size=2e5,
    inner_r=3.,
    central_width = 0.25,
    gap_radians = 2 * np.pi / 200.,
    kde_width = 0.5,
    bg_fmin=0.,
    bg_fmax=0.99,
    bg_vmin=None,
    bg_vmax=None,
    hyper_colour=consts.METHYLATION_DIRECTION_COLOURS['hyper'],
    hypo_colour = consts.METHYLATION_DIRECTION_COLOURS['hypo'],
    unmap_threshold_pct=10.,
    plot_border = True,
    bg_cmap=plt.cm.get_cmap('RdYlBu_r'),
):
    chroms = chrom_length.keys()

    if bg_vmin is None and bg_fmin is None:
        bg_vmin = 0.
    if bg_vmax is None and bg_fmax is None:
        raise ValueError("Must supply one of bg_vmax or bg_fmax.")

    if bg_density is not None:
        all_bg = []
        for chrom in chroms:
            this_bg = bg_density[chrom]
            all_bg.extend(this_bg.values)
        all_bg = sorted(all_bg)
        bg_mean = np.mean(all_bg)
        if bg_vmin is None:
            bg_vmin = all_bg[int(bg_fmin * len(all_bg))]
        if bg_vmax is None:
            bg_vmax = all_bg[int(bg_fmax * len(all_bg))]
        if bg_vmin > bg_mean:
            raise ValueError("The vmin value calculated for the background density is greater than the mean.")
        if bg_vmax < bg_mean:
            raise ValueError("The vmax value calculated for the background density is less than the mean.")

    sum_chrom_length = float(sum(chrom_length.values()))
    radians_per_bp = (2 * np.pi - gap_radians * len(chroms)) / sum_chrom_length

    outer_r = inner_r + central_width

    # generate all KDEs first, so we can rescale them correctly at plot time
    hypo_kdes = {}
    hyper_kdes = {}
    kde_th = {}

    curr_theta = 0.
    for chrom in chroms:
        xx = np.arange(1, chrom_length[chrom], window_size)
        if xx[-1] != chrom_length[chrom]:
            xx = np.concatenate((xx, [chrom_length[chrom]]))
        # this_bg = bg_density[chrom]
        # xx = np.array(this_bg.index.tolist() + [chrom_length[chrom]])
        dummy_res = np.ma.masked_all(xx.shape)
        this_hypo = dmr_loci_hypo[chrom]
        this_hyper = dmr_loci_hyper[chrom]

        if len(this_hypo):
            hypo_kdes[chrom] = fit_kde_dmr_location(this_hypo, xx, window_size, normed=False)
        else:
            hypo_kdes[chrom] = dummy_res

        if len(this_hyper):
            hyper_kdes[chrom] = fit_kde_dmr_location(this_hyper, xx, window_size, normed=False)
        else:
            hyper_kdes[chrom] = dummy_res

        th = curr_theta + xx * radians_per_bp
        kde_th[chrom] = th
        curr_theta = th[-1] + gap_radians


    # rescale KDEs
    # we want the density to reflect the number of DMRs

    hypo_kdes_rescaled = {}
    hyper_kdes_rescaled = {}

    hypo_max = max([t.max() for t in hypo_kdes.values()])
    hyper_max = max([t.max() for t in hyper_kdes.values()])

    for chrom in chroms:
        d_hypo = hypo_kdes[chrom]
        hypo_kdes_rescaled[chrom] = d_hypo / hypo_max * kde_width

        d_hyper = hyper_kdes[chrom]
        hyper_kdes_rescaled[chrom] = d_hyper / hyper_max * kde_width

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='polar')
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    # hack required to convert patches correctly to polar coords
    # https://github.com/matplotlib/matplotlib/issues/8521
    ax.bar(0, 1).remove()

    text_handles = {}
    border_patches = {}
    bg_handles = {}
    unmapped_handles = {}

    curr_theta = 0
    for chrom in chroms:
        if bg_density is not None:
            this_bg = bg_density[chrom]

            d_hypo = hypo_kdes_rescaled[chrom]
            d_hyper = hyper_kdes_rescaled[chrom]

            th = curr_theta + np.array(this_bg.index.tolist() + [chrom_length[chrom]]) * radians_per_bp
            tt = np.array([th, th])
            rr = np.zeros_like(tt) + inner_r
            rr[1] = outer_r

            cc = np.array([this_bg.values])

            norm = common.MidpointNormalize(midpoint=bg_mean, vmin=bg_vmin, vmax=bg_vmax)
            h_bg = ax.pcolor(tt, rr, cc, cmap=bg_cmap, norm=norm, vmin=bg_vmin, vmax=bg_vmax)
            bg_handles[chrom] = h_bg

        if unmapped_density is not None:
            this_unmapped = unmapped_density[chrom]
            this_unmapped_pct = this_unmapped / float(window_size) * 100.
            th = curr_theta + np.array(this_unmapped.index.tolist() + [chrom_length[chrom]]) * radians_per_bp

            uu = np.ma.masked_less(np.array([this_unmapped_pct.values]), unmap_threshold_pct)
            # since we don't care about the extent of unmapping, replace all values with a single one
            uu[~uu.mask] = 0.5
            h_unmapped = ax.pcolor(tt, rr, uu, cmap='Greys', vmax=1., vmin=0.)
            unmapped_handles[chrom] = h_unmapped

        if plot_border:
            this_patch = plt.Rectangle(
                [curr_theta, inner_r],
                width=chrom_length[chrom] * radians_per_bp,
                height=central_width,
                edgecolor='k',
                facecolor='none',
                linewidth=1.,
                zorder=999,
            )
            ax.add_artist(this_patch)
            border_patches[chrom] = this_patch

        # add chromosome name
        th_midpoint = (curr_theta + chrom_length[chrom] * radians_per_bp / 2.)
        h_text = ax.text(
            th_midpoint,
            outer_r + kde_width / 2.,
            chrom,
            horizontalalignment='center',
            verticalalignment='center'
        )
        text_handles[chrom] = h_text

        # plot the KDEs
        th = kde_th[chrom]
        ax.fill_between(th, y1=d_hyper + outer_r, y2=outer_r, alpha=0.9, color=hyper_colour)
        ax.fill_between(th, y1=inner_r, y2=inner_r - d_hypo, alpha=0.9, color=hypo_colour)

        ax.set_facecolor('w')

        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        curr_theta = th[-1] + gap_radians

    fig.tight_layout()

    return {
        'fig': fig,
        'ax': ax,
        'border_patches': border_patches,
        'text_handles': text_handles,
        'unmapped_handles': unmapped_handles,
        'bg_handles': bg_handles,
        'hypo_kde_rescaled': hypo_kdes_rescaled,
        'hyper_kde_rescaled': hyper_kdes_rescaled,
        'hypo_kde': hypo_kdes,
        'hyper_kde': hyper_kdes,
    }


# plot just one chromosome from one patient at a time for more detailed inspection
def plot_one_chrom_location(
        hyper_dmr_loci,
        hypo_dmr_loci,
        bg_density,
        chrom_len,
        unmapped_density=None,
        unmapped_threshold=10,
        window_size=5e5,
        bg_vmin=1,
        bg_vmax=10,
        cmap=plt.cm.get_cmap('RdYlBu_r'),
        plot_border=True,
):
    gs = plt.GridSpec(ncols=1, nrows=3, height_ratios=[6, 1, 6])
    fig = plt.figure(figsize=(6, 3.5))
    ax_hyper = fig.add_subplot(gs[0])
    ax_bg = fig.add_subplot(gs[1], sharex=ax_hyper)
    ax_hypo = fig.add_subplot(gs[2], sharex=ax_hyper)

    # plot bg density
    xx = np.array([bg_density.index.tolist() + [chrom_len]] * 2)
    yy = np.zeros_like(xx);
    yy[1] = 1.
    cc = np.ma.masked_less([bg_density.values], bg_vmin)

    norm = common.MidpointNormalize(midpoint=bg_density.mean(), vmin=bg_vmin, vmax=bg_vmax)
    ax_bg.pcolor(xx, yy, cc, cmap=cmap, norm=norm)
    ax_bg.set_xlim([-10, xx[0, -1] + 10])

    if unmapped_density is not None:
        uu = np.ma.masked_less(np.array([unmapped_density.values]), unmapped_threshold)
        # since we don't care about the extent of unmapping, replace all values with a single one
        uu[~uu.mask] = 0.8
        ax_bg.pcolor(xx, yy, uu, cmap='Greys', vmax=1., vmin=0.)

    if plot_border:
        # draw a border around the extent of the chromosome
        border = plt.Rectangle(
            [0, 0.01],
            chrom_len,
            .98,
            edgecolor='k',
            facecolor='none',
            linewidth=1.,
            zorder=100.
        )
        ax_bg.add_patch(border)

    # plot the histograms
    edges = np.arange(1, chrom_len, window_size)
    if edges[-1] != chrom_len:
        edges = np.concatenate((edges, [chrom_len]))
    hyper_to_plot, _ = np.histogram(hyper_dmr_loci, edges)
    ax_hyper.bar(
        edges[:-1],
        hyper_to_plot,
        width=edges[1:] - edges[:-1],
        color=consts.METHYLATION_DIRECTION_COLOURS['hyper'],
        edgecolor='none'
    )

    hypo_to_plot, _ = np.histogram(hypo_dmr_loci, edges)
    ax_hypo.bar(
        edges[:-1],
        hypo_to_plot,
        width=edges[1:] - edges[:-1],
        color=consts.METHYLATION_DIRECTION_COLOURS['hypo'],
        edgecolor='none'
    )

    ymax = max(hyper_to_plot.max(), hypo_to_plot.max())

    ax_bg.set_ylim([-.02, 1.02])
    ax_hypo.set_ylim([0, ymax * 1.02])
    ax_hyper.set_ylim([0, ymax * 1.02])
    plt.setp([ax_hypo, ax_hyper, ax_bg], facecolor='w')
    ax_hyper.set_ylabel('# Hyper DMRs')
    ax_hypo.set_ylabel('# Hypo DMRs')

    ax_bg.yaxis.set_visible(False)
    ax_bg.xaxis.set_visible(False)
    ax_hypo.xaxis.set_visible(False)
    ax_hyper.xaxis.set_visible(False)

    # hypo: needs to be plotted upside down
    ax_hypo.invert_yaxis()

    gs.update(hspace=0.02, right=0.99)

    return {
        'fig': fig,
        'gs': gs,
        'ax_bg': ax_bg,
        'ax_hyper': ax_hyper,
        'ax_hypo': ax_hypo,
    }


if __name__ == "__main__":
    pids = consts.PIDS
    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    # set this to True if output bed files are required (this is quite slow due to the large number of combinations)
    write_bed_files = False

    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    me_obj, anno = tsgd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
    me_data = me_obj.data
    me_meta = me_obj.meta
    me_meta.insert(0, 'patient_id', me_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

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

    # patient-specific DMRs
    pu_sets = list(setops.binary_combinations_sum_eq(len(pids), 1))[::-1]  # reverse order to get same order as pids
    venn_set, venn_ct = setops.venn_from_arrays(*[dmr_res_all[pid] for pid in pids])
    vs = dict([
                  (p, venn_set[q]) for p, q in zip(pids, pu_sets)
                  ])
    dmr_res_specific = dict([
                                (
                                    pid,
                                    dict([(t, dmr_res_all[pid][t]) for t in vs[pid]])
                                ) for pid in pids
                                ])

    chroms = [str(t) for t in range(1, 23)]
    fa_fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human/ensembl/GRCh38.release87/fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    )
    window_size = int(2e5)
    coord_summary_method = 'first'
    cg_locs, chrom_length = genomics.motif_locations(fa_fn, features=chroms)

    # reorder
    chrom_length = collections.OrderedDict([(k, chrom_length[k]) for k in chroms])

    dmr_s1_clusters = dmr_res_s1[pids[0]].clusters

    unmapped_density, _ = genomics.cg_content_windowed(fa_fn, features=chroms, window_size=window_size, motif='N')

    # get all DMRs (whether significant or not) to use as a null hypothesis
    tmp = get_binned_dmr_locations(
        dmr_res_s1.results,
        dmr_s1_clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=False,
        coord_summary_method=coord_summary_method
    )

    # these are the same for all patients, so just take the first
    dmr_loci_all = tmp['dmr_loci'][pids[0]]
    dmr_binned_all = tmp['dmr_binned'][pids[0]]

    #1 ) All DMRs
    tmp = get_binned_dmr_locations(
        dmr_res_all,
        dmr_s1_clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=True,
        coord_summary_method=coord_summary_method
    )
    dmr_loci_hyper = tmp['dmr_loci_hyper']
    dmr_loci_hypo = tmp['dmr_loci_hypo']


    hyper_all_fdr, hyper_all_stat = ks_test_dmr_locations(
        dmr_loci_hyper,
        dmr_loci_all
    )
    hyper_all_fdr = pd.DataFrame(hyper_all_fdr).loc[chroms, pids]

    hypo_all_fdr, hypo_all_stat = ks_test_dmr_locations(
        dmr_loci_hypo,
        dmr_loci_all
    )
    hypo_all_fdr = pd.DataFrame(hypo_all_fdr).loc[chroms, pids]

    vmin = .5
    vmax = 3.
    alpha = 0.1
    gs = plt.GridSpec(1, 3, width_ratios=[9, 9, 1])
    fig = plt.figure(figsize=(9.5, 5.5))
    ax_hyper = fig.add_subplot(gs[0])
    ax_hypo = fig.add_subplot(gs[1], sharex=ax_hyper, sharey=ax_hyper)
    ax_cbar = fig.add_subplot(gs[2])

    hyper_to_plot = -np.log10(hyper_all_fdr)
    hyper_to_plot = hyper_to_plot[hyper_to_plot > -np.log10(alpha)]
    hypo_to_plot = -np.log10(hypo_all_fdr)
    hypo_to_plot = hypo_to_plot[hypo_to_plot > -np.log10(alpha)]

    sns.heatmap(
        hyper_to_plot,
        cmap='Reds',
        vmin=vmin,
        vmax=vmax,
        ax=ax_hyper,
        cbar=False,
        linewidths=1.,
        linecolor='w'
    )
    sns.heatmap(
        hypo_to_plot,
        cmap='Reds',
        vmin=vmin,
        vmax=vmax,
        ax=ax_hypo,
        linewidths=1.,
        linecolor='w',
        cbar_ax=ax_cbar
    )
    plt.setp(ax_hypo.yaxis.get_ticklabels(), visible=False)
    ax_hyper.set_ylabel('Chromosome')
    ax_hyper.set_xlabel('Patient')
    ax_hyper.set_title('Hypermethylated DMRs')
    ax_hypo.set_xlabel('Patient')
    ax_hypo.set_title('Hypomethylated DMRs')
    ax_cbar.set_ylabel('$-\log_{10}(p)$')
    plt.setp(ax_hyper.yaxis.get_ticklabels(), rotation=0)

    gs.update(left=0.06, right=0.94, top=0.95, wspace=0.1)
    fig.savefig(os.path.join(outdir, "ks_fdr_all_dmr_location.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "ks_fdr_all_dmr_location.tiff"), dpi=200)

    # repeat for patient-specific DMRs

    tmp = get_binned_dmr_locations(
        dmr_res_specific,
        dmr_s1_clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=True,
        coord_summary_method=coord_summary_method
    )

    dmr_loci_hyper_specific = tmp['dmr_loci_hyper']
    dmr_loci_hypo_specific = tmp['dmr_loci_hypo']

    hyper_specific_fdr, hyper_specific_stat = ks_test_dmr_locations(
        dmr_loci_hyper_specific,
        dmr_loci_all
    )
    hyper_specific_fdr = pd.DataFrame(hyper_specific_fdr).loc[chroms, pids]

    hypo_specific_fdr, hypo_specific_stat = ks_test_dmr_locations(
        dmr_loci_hypo_specific,
        dmr_loci_all
    )
    hypo_specific_fdr = pd.DataFrame(hypo_specific_fdr).loc[chroms, pids]

    vmin = .5
    vmax = 3.
    alpha = 0.1
    gs = plt.GridSpec(1, 3, width_ratios=[9, 9, 1])
    fig = plt.figure(figsize=(9.5, 5.5))
    ax_hyper = fig.add_subplot(gs[0])
    ax_hypo = fig.add_subplot(gs[1], sharex=ax_hyper, sharey=ax_hyper)
    ax_cbar = fig.add_subplot(gs[2])

    hyper_to_plot_specific = -np.log10(hyper_specific_fdr)
    hyper_to_plot_specific = hyper_to_plot_specific[hyper_to_plot_specific > -np.log10(alpha)]
    hypo_to_plot_specific = -np.log10(hypo_specific_fdr)
    hypo_to_plot_specific = hypo_to_plot_specific[hypo_to_plot_specific > -np.log10(alpha)]

    sns.heatmap(
        hyper_to_plot_specific,
        cmap='Reds',
        vmin=vmin,
        vmax=vmax,
        ax=ax_hyper,
        cbar=False,
        linewidths=1.,
        linecolor='w'
    )
    sns.heatmap(
        hypo_to_plot_specific,
        cmap='Reds',
        vmin=vmin,
        vmax=vmax,
        ax=ax_hypo,
        linewidths=1.,
        linecolor='w',
        cbar_ax=ax_cbar
    )
    plt.setp(ax_hypo.yaxis.get_ticklabels(), visible=False)
    ax_hyper.set_ylabel('Chromosome')
    ax_hyper.set_xlabel('Patient')
    ax_hyper.set_title('Hypermethylated DMRs')
    ax_hypo.set_xlabel('Patient')
    ax_hypo.set_title('Hypomethylated DMRs')
    ax_cbar.set_ylabel('$-\log_{10}(p)$')
    plt.setp(ax_hyper.yaxis.get_ticklabels(), rotation=0)

    gs.update(left=0.06, right=0.94, top=0.95, wspace=0.1)
    fig.savefig(os.path.join(outdir, "ks_fdr_specific_dmr_location.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "ks_fdr_specific_dmr_location.tiff"), dpi=200)

    # generate a polar plot illustrating DMR locations

    # truncate the intended colourmap
    # this softens the appearance by removing the deepest colours at either end of the map
    cmap = plt.cm.get_cmap('RdYlBu_r')
    new_cmap_list = [cmap(t) for t in np.linspace(0.15, 0.85, 256)]
    new_cmap = cmap.from_list('RdYlBu_r_truncated', new_cmap_list)

    # bands for annotating statistical significance
    fdr_annot = [
        (0.001, "***"),
        (0.01, "**"),
        (0.05, "*"),
    ]
    annot_offset = 0.25
    inner_r = 3.
    kde_width = 0.5
    annot_fontdict = {
        'fontsize': 16,
        'weight': 'bold',
        'horizontalalignment': 'center',
        'verticalalignment': 'center'
    }
    hyper_kdes = {}
    hypo_kdes = {}

    for pid in pids:
        plt_dict = polar_distribution_plot(
            dmr_loci_hyper[pid],
            dmr_loci_hypo[pid],
            chrom_length,
            bg_density=dmr_binned_all,
            unmapped_density=unmapped_density,
            window_size=5e5,
            bg_cmap=new_cmap,
            central_width=.3,
            inner_r=inner_r,
            kde_width=kde_width
        )
        hyper_kdes[pid] = plt_dict['hyper_kde']
        hypo_kdes[pid] = plt_dict['hypo_kde']

        ax = plt_dict['ax']
        this_fdr_hyper = hyper_all_fdr[pid]
        this_fdr_hypo = hypo_all_fdr[pid]
        for chrom in chroms:
            for thresh, annot in fdr_annot:
                if this_fdr_hyper[chrom] < thresh:
                    # add annotation
                    loc = plt_dict['text_handles'][chrom].get_position()
                    ax.text(
                        loc[0],
                        loc[1] + annot_offset,
                        annot,
                        color=consts.METHYLATION_DIRECTION_COLOURS['hyper'],
                        **annot_fontdict
                    )
                    break
            for thresh, annot in fdr_annot:
                if this_fdr_hypo[chrom] < thresh:
                    # add annotation
                    loc = plt_dict['text_handles'][chrom].get_position()
                    ax.text(
                        loc[0],
                        inner_r - kde_width / 2.,
                        annot,
                        color=consts.METHYLATION_DIRECTION_COLOURS['hypo'],
                        **annot_fontdict
                    )
                    break

        fig = plt_dict['fig']

        fig.savefig(os.path.join(outdir, "all_dmrs_polar_distribution_plot_%s.png" % pid), dpi=200)
        fig.savefig(os.path.join(outdir, "all_dmrs_polar_distribution_plot_%s.tiff" % pid), dpi=200)
        fig.savefig(os.path.join(outdir, "all_dmrs_polar_distribution_plot_%s.pdf" % pid), dpi=200)

    # for each chromosome / patient of interest, try to determine whether there are any 'hotspots' causing significant
    # non-randomness
    subwindow_size = int(5e6)
    # rolling window shift
    roll_dw = 5e5
    ks_subwindows_hyper = {}
    ks_subwindows_hypo = {}

    pool = mp.Pool()
    jobs = {}

    for pid in pids:
        ks_subwindows_hyper[pid] = {}
        ks_subwindows_hypo[pid] = {}
        for chrom in chroms:
            ks_subwindows_hyper[pid][chrom] = {}
            ks_subwindows_hypo[pid][chrom] = {}
            this_all = np.array(dmr_loci_all[chrom])
            this_hyper = np.array(dmr_loci_hyper[pid][chrom])
            this_hypo = np.array(dmr_loci_hypo[pid][chrom])

            this_res_hyper = {}
            this_res_hypo = {}

            # rolling window approach
            s = 0
            e = s + subwindow_size
            while e < chrom_length[chrom]:
                sub_all = this_all[(this_all >= s) & (this_all < e)]
                sub_hyper = this_hyper[(this_hyper >= s) & (this_hyper < e)]
                sub_hypo = this_hypo[(this_hypo >= s) & (this_hypo < e)]
                if len(sub_hypo) > 0:
                    jobs[(pid, chrom, (s, e), 'hypo')] = pool.apply_async(stats.ks_2samp, args=(sub_hypo, sub_all))
                    # this_res_hypo[(s, e)] = stats.ks_2samp(sub_hypo, sub_all)[1]
                if len(sub_hyper) > 0:
                    jobs[(pid, chrom, (s, e), 'hyper')] = pool.apply_async(stats.ks_2samp, args=(sub_hyper, sub_all))
                    # this_res_hyper[(s, e)] = stats.ks_2samp(sub_hyper, sub_all)[1]
                s += roll_dw
                e = s + subwindow_size

            # ks_subwindows_hyper[pid][chrom] = pd.Series(this_res_hyper)
            # ks_subwindows_hypo[pid][chrom] = pd.Series(this_res_hypo)
    pool.close()
    pool.join()

    for k, v in jobs.iteritems():
        if k[-1] == 'hyper':
            ks_subwindows_hyper[k[0]][k[1]][k[2]] = v.get()[1]
        elif k[-1] == 'hypo':
            ks_subwindows_hypo[k[0]][k[1]][k[2]] = v.get()[1]
        else:
            raise KeyError("Unrecognised final key %s." % str(k[-1]))

    for pid in pids:
        for chrom in chroms:
            ks_subwindows_hyper[pid][chrom] = pd.Series(ks_subwindows_hyper[pid][chrom])
            ks_subwindows_hypo[pid][chrom] = pd.Series(ks_subwindows_hypo[pid][chrom])

    # overlay DMRs with features?
    # use utils.genomics.GtfAnnotation (.region method)
    # or just lookup in the annotation file?


    # plot single chromosomes in individual patients to illustrate examples of interest.

    pid = '018'
    chrom = '14'

    pid_chroms = [
        ('018', '6'),
        ('018', '14'),
    ]

    for pid, chrom in pid_chroms:

        plt_dict = plot_one_chrom_location(
            dmr_loci_hyper[pid][chrom],
            dmr_loci_hypo[pid][chrom],
            dmr_binned_all[chrom],
            chrom_length[chrom],
            unmapped_density=unmapped_density[chrom] / window_size * 100.,
            cmap=new_cmap
        )
        fig = plt_dict['fig']
        fig.savefig(os.path.join(outdir, "patient_%s_chrom_%s_locations.png" % (pid, chrom)), dpi=200)

        # second plot showing the KS 'substatistic'
        fig, axs = plt.subplots(nrows=2, figsize=(8, 4), sharex=True)
        x_hyper = [np.mean(t) for t in ks_subwindows_hyper[pid][chrom].index]
        y_hyper = -np.log10(ks_subwindows_hyper[pid][chrom].values)
        x_hypo = [np.mean(t) for t in ks_subwindows_hypo[pid][chrom].index]
        y_hypo = -np.log10(ks_subwindows_hypo[pid][chrom].values)

        axs[0].plot(x_hyper, y_hyper, color='0.7', zorder=11)
        axs[0].scatter(x_hyper, y_hyper, marker='o', c=y_hyper, edgecolor='none', cmap='copper_r', zorder=12)
        axs[1].plot(x_hypo, y_hypo, color='0.7', zorder=11)
        axs[1].scatter(x_hypo, y_hypo, marker='o', c=y_hypo, edgecolor='none', cmap='copper_r', zorder=12)
        # set same ylims
        ymax = max(y_hyper.max(), y_hypo.max())
        plt.setp(axs, ylim=[-.5, 1.1 * ymax])
        axs[1].invert_yaxis()

        axs[0].set_xlim([0, 1.02 * max(max(x_hyper), max(x_hypo))])

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "patient_%s_chrom_%s_sub_ks.png" % (pid, chrom)), dpi=200)
