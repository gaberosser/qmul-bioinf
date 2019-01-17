from plotting import bar, common, pie
from methylation import loader, dmr, process
import pandas as pd
from statsmodels.sandbox.stats import multicomp
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
from plotting import common

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
logger = log.get_console_logger()

# FIXME: this is a hack to try to avoid non-thread safe TKinter issue, which results in a segfault when we try
# to use multiprocessing. This requires ipython NOT to be run with --matplotlib
matplotlib.use('Agg')


"""
Here we analyse the direction and genomic locations of DMRs in the individual patients.
We note a bias in the direction, which is amplified considerably when we consider patient-specific DMRs.
"""


def dmr_median_delta_by_chromosome(
        dmr_res,
        clusters,
        pids=consts.PIDS,
):
    median_deltas = {}
    for pid in pids:
        median_deltas[pid] = collections.defaultdict(list)
        this_dmrs = dmr_res[pid]
        for cid, attr in dmr_res[pid].items():
            ftr = clusters[cid].chr
            median_deltas[pid][ftr].append(attr['median_change'])

    return median_deltas


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

        dmr_regions_hypo = {}
        dmr_regions_hyper = {}
    else:
        dmr_loci = {}
        dmr_loci_binned = {}
        dmr_regions = {}

    for pid in dmr_res:
        if split_by_direction:
            dmr_loci_binned_hypo[pid] = {}
            this_loci_hypo = collections.defaultdict(list)
            dmr_loci_binned_hyper[pid] = {}
            this_loci_hyper = collections.defaultdict(list)

            this_regions_hypo = collections.defaultdict(list)
            this_regions_hyper = collections.defaultdict(list)
        else:
            dmr_loci_binned[pid] = {}
            this_loci = collections.defaultdict(list)
            this_regions = collections.defaultdict(list)

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
                    this_regions_hyper[pc.chr].append((min(pc.coord_list), max(pc.coord_list)))
                else:
                    this_loci_hypo[pc.chr].append(the_coord)
                    this_regions_hypo[pc.chr].append((min(pc.coord_list), max(pc.coord_list)))
            else:
                this_loci[pc.chr].append(the_coord)
                this_regions[pc.chr].append((min(pc.coord_list), max(pc.coord_list)))
        if split_by_direction:
            dmr_loci_hyper[pid] = this_loci_hyper
            dmr_loci_hypo[pid] = this_loci_hypo
            dmr_regions_hyper[pid] = this_regions_hyper
            dmr_regions_hypo[pid] = this_regions_hypo
        else:
            dmr_loci[pid] = this_loci
            dmr_regions[pid] = this_regions

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
            'dmr_regions_hyper': dmr_regions_hyper,
            'dmr_regions_hypo': dmr_regions_hypo,
        }
    else:
        return {
            'dmr_loci': dmr_loci,
            'dmr_binned': dmr_loci_binned,
            'dmr_regions': dmr_regions,
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
        pval_hyper=None,
        pval_hypo=None,
        unmapped_threshold=10,
        window_size=5e5,
        bg_vmin=1,
        bg_vmax=10,
        cmap=plt.cm.get_cmap('RdYlBu_r'),
        plot_border=True,
):
    """

    :param hyper_dmr_loci:
    :param hypo_dmr_loci:
    :param bg_density:
    :param chrom_len:
    :param unmapped_density:
    :param pval_hyper, pval_hypo: If supplied, these are plotted on the same axes as the DMR density histogram (with a
    2nd y axis).
    :param unmapped_threshold:
    :param window_size:
    :param bg_vmin:
    :param bg_vmax:
    :param cmap:
    :param plot_border:
    :return:
    """
    gs = plt.GridSpec(ncols=1, nrows=3, height_ratios=[6, 1, 6])
    fig = plt.figure(figsize=(6, 3.5))
    ax_hyper = fig.add_subplot(gs[0])
    ax_bg = fig.add_subplot(gs[1], sharex=ax_hyper)
    ax_hypo = fig.add_subplot(gs[2], sharex=ax_hyper)

    twin_axes = False
    out = {
        'ax_bg': ax_bg,
        'ax_hyper': ax_hyper,
        'ax_hypo': ax_hypo,
        'fig': fig,
        'gs': gs,
    }

    # plot bg density
    xx = np.array([bg_density.index.tolist() + [chrom_len]] * 2)
    yy = np.zeros_like(xx);
    yy[1] = 1.
    cc = np.ma.masked_less([bg_density.values], bg_vmin)

    norm = common.MidpointNormalize(midpoint=bg_density.mean(), vmin=bg_vmin, vmax=bg_vmax)
    ax_bg.pcolor(xx, yy, cc, cmap=cmap, norm=norm)
    ax_bg.set_xlim([-10, chrom_len + 10])

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

    pval_ymax = None
    if pval_hyper is not None and pval_hypo is not None:
        # unify the colour map and scale
        y_hyper = -np.log10(pval_hyper.values)
        y_hypo = -np.log10(pval_hypo.values)
        pval_ymax = max(np.nanmax(y_hyper), np.nanmax(y_hypo))

    if pval_hyper is not None:
        twin_axes = True
        ax_p_hyper = ax_hyper.twinx()
        out['ax_p_hyper'] = ax_p_hyper
        # turn off gridlines
        ax_p_hyper.grid(False)
        x_hyper = [np.mean(t) for t in pval_hyper.index]
        y_hyper = -np.log10(pval_hyper.values)
        ax_p_hyper.plot(x_hyper, y_hyper, color='k', alpha=0.35, linewidth=1.)
        vmax = None
        if pval_ymax is not None:
            vmax = pval_ymax
        ax_p_hyper.scatter(
            x_hyper,
            y_hyper,
            c=y_hyper,
            alpha=0.5,
            edgecolor='none',
            s=10,
            cmap='copper_r',
            vmin=0,
            vmax=vmax
        )
        ax_p_hyper.set_ylabel("$-\log_{10}(p_{\mathrm{local}})$")

    hypo_to_plot, _ = np.histogram(hypo_dmr_loci, edges)
    ax_hypo.bar(
        edges[:-1],
        hypo_to_plot,
        width=edges[1:] - edges[:-1],
        color=consts.METHYLATION_DIRECTION_COLOURS['hypo'],
        edgecolor='none'
    )

    if pval_hypo is not None:
        twin_axes = True
        ax_p_hypo = ax_hypo.twinx()
        out['ax_p_hypo'] = ax_p_hypo
        # turn off gridlines
        ax_p_hypo.grid(False)
        x_hypo = [np.mean(t) for t in pval_hypo.index]
        y_hypo = -np.log10(pval_hypo.values)
        ax_p_hypo.plot(x_hypo, y_hypo, color='k', alpha=0.35, linewidth=1.)
        vmax = None
        if pval_ymax is not None:
            vmax = pval_ymax
        ax_p_hypo.scatter(
            x_hypo,
            y_hypo,
            c=y_hypo,
            alpha=0.5,
            edgecolor='none',
            s=10,
            cmap='copper_r',
            vmin=0,
            vmax=vmax
        )
        ax_p_hypo.set_ylabel("$-\log_{10}(p_{\mathrm{local}})$")

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

    if pval_ymax is not None:
        # maintain symmetric y axis limits
        ax_p_hyper.set_ylim([0, pval_ymax * 1.02])
        ax_p_hypo.set_ylim([0, pval_ymax * 1.02])

    # hypo: needs to be plotted upside down
    ax_hypo.invert_yaxis()
    if pval_hypo is not None:
        ax_p_hypo.invert_yaxis()

    # set x limits
    ax_bg.set_xlim([-10, chrom_len + 10])

    if twin_axes:
        gs.update(hspace=0.02, right=0.9)
    else:
        gs.update(hspace=0.02, right=0.99)

    return out


def classify_dmr_regions(
        gtf_db,
        clusters,
        gene_biotypes=('protein_coding',),
        features=tuple([str(t) for t in range(1, 23)]),
        njobs=None
):
    gene_biotypes = set(gene_biotypes)
    features = set(features)

    region_to_cluster_id = {}
    regions = []

    for i, pc in clusters.items():
        if pc.chr in features:
            foo = (pc.chr,) + pc.coord_range
            regions.append(foo)
            region_to_cluster_id[foo] = i

    cluster_id_to_region = dict([x[::-1] for x in region_to_cluster_id.items()])

    lookup_result = genomics.multiple_region_lookup(regions, gtf_db, njobs=njobs)

    # reindex
    lookup_result = dict([
        (region_to_cluster_id[reg], val) for reg, val in lookup_result.items()
    ])

    classification = {}
    for cid, ftr_arr in lookup_result.items():
        ftr_arr = [x for x in ftr_arr if x['gene_biotype'][0] in gene_biotypes]
        if len(ftr_arr) == 0:
            classification[cid] = 'intergenic'
        else:
            typs = np.array([x.featuretype for x in ftr_arr])
            if (typs == 'exon').any():
                for i in np.where(typs == 'exon')[0]:
                    if ftr_arr[i]['exon_number'][0] == '1':
                        classification[cid] = 'first_exon'
                        continue
                    classification[cid] = 'exon'
            elif (typs == 'gene').any():
                classification[cid] = 'intron'
    return lookup_result, classification, cluster_id_to_region


def pie_chart_location(
        class_bg,
        class_patients,
        lookup_cats,
        colours_by_lookup,
        cat_labels,
        pids=consts.PIDS,
        chroms=tuple([str(t) for t in range(1, 23)]),
        scale_by_number='area',
        inner_radius=.5,
        segment_width=0.5,
        min_segment_width=0.3
):
    min_inner_radius = inner_radius * min_segment_width / float(segment_width)
    # construct the combination over all chromosomes for each patient
    class_all = collections.defaultdict(dict)
    for pid, d in class_patients.items():
        for d2 in d.values():
            class_all[pid].update(d2)

    # construct the combination over all chromosomes for the BG
    class_bg_all = {}
    for d in class_bg.values():
        class_bg_all.update(d)

    if scale_by_number is not None:
        n_all_chrom = [pd.Series(class_all[pid]).value_counts().reindex(lookup_cats).fillna(0.).sum() for pid in pids]
        nmax_all_chrom = float(max(n_all_chrom))

        n_indiv = np.array([
            [pd.Series(class_patients[pid][chrom]).value_counts().reindex(lookup_cats).fillna(0.).sum() for pid in pids]
            for chrom in chroms
        ])
        nmax = float(n_indiv.max().max())

        n_bg = [pd.Series(class_bg[chrom]).value_counts().reindex(lookup_cats).fillna(0.).sum() for chrom in chroms]
        nmax_bg = float(max(n_bg))

    fig = plt.figure(figsize=(5.5, 6.8))
    gs = plt.GridSpec(nrows=len(chroms) + 1, ncols=len(pids) + 1)
    axs = [[]]

    # top left axis
    ax = fig.add_subplot(gs[0, 0], facecolor='none')

    the_dat = pd.Series(class_bg_all).value_counts().reindex(lookup_cats).fillna(0.)
    if the_dat.sum() > 0:
        # no size scaling for this pie chart (too aggregated)
        ax.pie(
            the_dat,
            colors=[colours_by_lookup[t] for t in lookup_cats],
            radius=inner_radius + segment_width,
            wedgeprops={'edgecolor': 'k', 'width': segment_width},
            startangle=90,
        )
        ax.set_aspect('equal')
    else:
        ax.yaxis.set_ticks([])
        ax.tick_params(top='off', bottom='off', left='off', right='off', labelcolor='none')
        ax.grid(False)
    ax.set_title('BG')
    ax.set_ylabel('All', rotation=0, verticalalignment='center', horizontalalignment='right')
    axs[0].append(ax)

    # first row: summary over all chromosomes
    for i, pid in enumerate(pids):
        ax = fig.add_subplot(gs[0, i + 1])
        axs[0].append(ax)
        the_dat = pd.Series(class_all[pid]).value_counts().reindex(lookup_cats).fillna(0.)
        if the_dat.sum() > 0:
            # compute normalisation const (relative to 'all')
            k = 1.
            if scale_by_number == 'area':
                k = (the_dat.sum() / nmax_all_chrom) ** .5
            elif scale_by_number == 'radius':
                k = the_dat.sum() / nmax_all_chrom

            w_eff = max(segment_width * k, min_segment_width)
            inner_r_eff = max(inner_radius * k, min_inner_radius)
            ax.pie(
                the_dat,
                colors=[colours_by_lookup[t] for t in lookup_cats],
                radius=inner_r_eff + w_eff,
                wedgeprops={'edgecolor': 'k', 'width': w_eff},
                startangle=90,
            )
            ax.set_aspect('equal')
            ax.set_title(pid)
        else:
            ax.set_visible(False)

    for j, chrom in enumerate(chroms):
        axs.append([])
        ax = fig.add_subplot(gs[j + 1, 0])
        axs[j + 1].append(ax)
        ax.set_ylabel(chrom, rotation=0, verticalalignment='center', horizontalalignment='right')

        the_dat = pd.Series(class_bg[chrom]).value_counts().reindex(lookup_cats).fillna(0.)
        if the_dat.sum() > 0:
            # compute normalisation const (relative to 'bg')
            k = 1.
            if scale_by_number == 'area':
                k = (the_dat.sum() / nmax_bg) ** .5
            elif scale_by_number == 'radius':
                k = the_dat.sum() / nmax_bg

            w_eff = max(segment_width * k, min_segment_width)
            inner_r_eff = max(inner_radius * k, min_inner_radius)

            ax.pie(
                the_dat,
                colors=[colours_by_lookup[t] for t in lookup_cats],
                radius=inner_r_eff + w_eff,
                wedgeprops={'edgecolor': 'k', 'width': w_eff},
                startangle=90,
            )
            ax.set_aspect('equal')
        else:
            ax.set_visible(False)

        for i, pid in enumerate(pids):
            ax = fig.add_subplot(gs[j + 1, i + 1])
            axs[j + 1].append(ax)

            the_dat = pd.Series(class_patients[pid][chrom]).value_counts().reindex(lookup_cats).fillna(0.)
            if the_dat.sum() > 0:
                # compute normalisation const (relative to 'patients')
                k = 1.
                if scale_by_number == 'area':
                    k = (the_dat.sum() / nmax) ** .5
                elif scale_by_number == 'radius':
                    k = the_dat.sum() / nmax

                w_eff = max(segment_width * k, min_segment_width)
                inner_r_eff = max(inner_radius * k, min_inner_radius)

                ax.pie(
                    the_dat,
                    colors=[colours_by_lookup[t] for t in lookup_cats],
                    radius=inner_r_eff + w_eff,
                    wedgeprops={'edgecolor': 'k', 'width': w_eff},
                    startangle=90,
                )
                ax.set_aspect('equal')
            else:
                ax.set_visible(False)

    gs.update(hspace=0, wspace=0, bottom=0.01, left=0.06, top=0.95, right=0.78)
    # custom legend
    legend_dict = {
        '': collections.OrderedDict(
            [(
                 cat_labels[cat],
                 {
                     'class': 'patch',
                     'facecolor': colours_by_lookup[cat],
                     'edgecolor': 'k',
                     'linewidth': 0.5
                 }
             ) for cat in lookup_cats
             ])
    }
    common.add_custom_legend(
        axs[0][-1],
        legend_dict,
        loc_outside=True,
    )

    return fig, axs


def pairwise_distance_by_class_count(
        count_dat,
        inner_dist_fun=cluster.hierarchy.distance.cosine,
        outer_metric=np.linalg.norm
):
    rows = count_dat.keys()
    cols = count_dat.values()[0].columns

    pdist_row = pd.DataFrame(
        index=rows,
        columns=rows,
    )

    pdist_col = pd.DataFrame(
        index=cols,
        columns=cols
    )

    # loop: keys
    for k1 in rows:
        for k2 in rows:
            df1 = count_dat[k1]
            df2 = count_dat[k2]
            # reduce two matrices to a vector using inner distance func
            inner_vec = pd.Series(
                dict([
                         (col, inner_dist_fun(df1[col], df2[col])) for col in cols
                         ])
            )
            # reduce this to a single distance value using outer distance fun
            pdist_row.loc[k1, k2] = outer_metric(inner_vec.dropna())

    # loop: columns
    for k1 in cols:
        for k2 in cols:
            df1 = pd.DataFrame(
                dict([
                         (r, count_dat[r][k1]) for r in rows
                         ])
            )
            df2 = pd.DataFrame(
                dict([
                         (r, count_dat[r][k2]) for r in rows
                         ])
            )
            # reduce two matrices to a vector using inner distance func
            inner_vec = pd.Series(
                dict([
                         (col, inner_dist_fun(df1[row], df2[row])) for row in rows
                         ])
            )
            # reduce this to a single distance value using outer distance fun
            pdist_col.loc[k1, k2] = outer_metric(inner_vec.dropna())

    return pdist_row, pdist_col


def dmr_direction_by_chromosome_pie_array(
        dat,
        number,
        outer_sm_hypo,
        outer_sm_hyper,
        pids=consts.PIDS,
        chroms=tuple([str(t) for t in range(1, 23)]),
        edges=None,
        inner_colours=None,
        inner_radius=0.5,
        segment_width=0.5,
        scale_by_number='area',
        min_segment_width=0.2,
        wedgeprops=None
):
    if wedgeprops is None:
        wedgeprops = {
            'edgecolor': 'k',
            'linewidth': .5
        }
    if edges is None:
        edges = np.array([-np.inf, -5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5., np.inf])
    if inner_colours is None:
        inner_colours = {
            'Hypermethylated': consts.METHYLATION_DIRECTION_COLOURS['hyper'],
            'Hypomethylated': consts.METHYLATION_DIRECTION_COLOURS['hypo'],
        }

    # result of a call to pie
    # initially None and set only once
    res = None

    # combine results over chromosomes
    dat_all = {}
    for pid in pids:
        dat_all[pid] = []
        for v in dat[pid].values():
            dat_all[pid].extend(v)

    number_all = number.sum(axis=0)

    # construct a dictionary of all colours
    colours = dict(inner_colours)
    for e0, e1 in zip(edges[:-1], edges[1:]):
        if e0 < 0:
            sm = outer_sm_hypo
        else:
            sm = outer_sm_hyper
        if e0 == edges[0]:
            lbl = r"$\Delta M < %d$" % e1
        elif e1 == edges[-1]:
            lbl = r"$\Delta M > %d$" % e0
        else:
            lbl = r"$%d \leq \Delta M < %d$" % (e0, e1)
        colours[lbl] = sm.to_rgba(np.sign(e0) * np.abs([e0, e1]).min())

    def bin_one(vals, edges):
        binned, _ = np.histogram(vals, edges)
        to_plot = collections.OrderedDict()
        to_plot['Hypermethylated'] = collections.OrderedDict()
        to_plot['Hypomethylated'] = collections.OrderedDict()
        for e0, e1, b in zip(edges[:-1], edges[1:], binned):
            if e0 == edges[0]:
                lbl = r"$\Delta M < %d$" % e1
            elif e1 == edges[-1]:
                lbl = r"$\Delta M > %d$" % e0
            else:
                lbl = r"$%d \leq \Delta M < %d$" % (e0, e1)

            if e0 > 0:
                dest = to_plot['Hypermethylated']
            elif e1 < 0:
                dest = to_plot['Hypomethylated']
            else:
                continue
            dest[lbl] = b
        return to_plot

    if scale_by_number is not None:
        nmax = float(number.max().max())
        nmax_all = float(number.sum(axis=0).max())

    gs = plt.GridSpec(
        nrows=len(chroms) + 1,
        ncols=len(pids),
    )
    fig = plt.figure(figsize=(7., 10.))
    axs = [[] for i in range(len(chroms) + 1)]
    sharex = None
    # first row: all
    for i, pid in enumerate(pids):
        ax = fig.add_subplot(gs[0, i], sharex=sharex, sharey=sharex)
        ax.set_aspect('equal')
        ax.set_title(pid, fontsize=12, horizontalalignment='center')
        if i == 0:
            ax.set_ylabel(
                'All',
                fontsize=12,
                verticalalignment='center',
                rotation=0,
                horizontalalignment='right'
            )
        axs[0].append(ax)
        if sharex is None:
            sharex = ax

        if len(dat_all[pid]) > 0:
            to_plot = bin_one(dat_all[pid], edges)
            # compute normalisation const
            k = 1.
            if scale_by_number == 'area':
                k = (number_all.loc[pid] / nmax_all) ** .5
            elif scale_by_number == 'radius':
                k = number_all.loc[pid] / nmax_all

            w_eff = max(segment_width * k, min_segment_width)
            inner_r_eff = inner_radius * k

            tmp = pie.nested_pie_chart(
                to_plot,
                colours,
                inner_radius=inner_r_eff,
                width_per_level=w_eff,
                ax=ax,
                legend_entries=None,
                **wedgeprops
            )
            if res is None:
                res = tmp
        else:
            ax.set_visible(False)

        for j, chrom in enumerate(chroms):
            ax = fig.add_subplot(gs[j + 1, i], sharex=sharex, sharey=sharex)
            ax.set_aspect('equal')
            if i == 0:
                ax.set_ylabel(
                    chrom,
                    fontsize=12,
                    verticalalignment='center',
                    rotation=0,
                    horizontalalignment='right'
                )
            axs[j + 1].append(ax)

            if len(dat[pid][chrom]) > 0:
                to_plot = bin_one(dat[pid][chrom], edges)
                # compute normalisation const
                k = 1.
                if scale_by_number == 'area':
                    k = (number.loc[chrom, pid] / nmax) ** .5
                elif scale_by_number == 'radius':
                    k = number.sum(axis=0).loc[pid] / nmax

                w_eff = max(segment_width * k, min_segment_width)
                inner_r_eff = inner_radius * k

                tmp = pie.nested_pie_chart(
                    to_plot,
                    colours,
                    inner_radius=inner_r_eff,
                    width_per_level=w_eff,
                    ax=ax,
                    legend_entries=None,
                    **wedgeprops
                )
                if res is None:
                    res = tmp
            else:
                ax.set_visible(False)

    # fix axis limits across all plots
    axs[0][0].set_xlim(np.array([-1, 1]) * (inner_radius + 2.5 * segment_width))

    # add legend outside axes
    patches = res['patches']
    legend_dict = collections.OrderedDict()
    legend_dict['Hyper'] = collections.OrderedDict()
    legend_dict['Hypo'] = collections.OrderedDict()
    for k in to_plot['Hypermethylated']:
        legend_dict['Hyper'][k] = {
            'class': 'patch',
            'facecolor': patches[1][k].get_facecolor(),
            'edgecolor': patches[1][k].get_edgecolor(),
            'linewidth': patches[1][k].get_linewidth(),
        }
    for k in to_plot['Hypomethylated'].keys()[::-1]:
        legend_dict['Hypo'][k] = {
            'class': 'patch',
            'facecolor': patches[1][k].get_facecolor(),
            'edgecolor': patches[1][k].get_edgecolor(),
            'linewidth': patches[1][k].get_linewidth(),
        }

    common.add_custom_legend(
        axs[int((len(chroms) + 1) / 2.)][-1],
        legend_dict,
        loc_outside=True
    )

    gs.update(left=0.05, bottom=0.02, top=0.97, right=0.72, hspace=0.06, wspace=0.06)

    return {
        'fig': fig,
        'axs': axs,
        'gs': gs,
        'legend_dict': legend_dict,
        'pie_res': res,
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

    clusters = dmr_res_s1[pids[0]].clusters

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

    # look at the distribution of DMR directions by chromosome
    median_delta_all_by_chrom = dmr_median_delta_by_chromosome(dmr_res_all, clusters)
    median_delta_specific_by_chrom = dmr_median_delta_by_chromosome(dmr_res_specific, clusters)

    num_dmr_all_by_chrom = pd.DataFrame(median_delta_all_by_chrom).applymap(len).loc[chroms, pids]
    num_dmr_specific_by_chrom = pd.DataFrame(median_delta_specific_by_chrom).applymap(
        lambda x: len(x) if hasattr(x, '__len__') else 0
    ).loc[chroms, pids]

    # integrate over chromosomes for overall picture
    median_delta_all = dict([
        (pid, reduce(lambda x, y: x + y, median_delta_all_by_chrom[pid].values())) for pid in pids
    ])
    median_delta_specific = dict([
        (pid, reduce(lambda x, y: x + y, median_delta_specific_by_chrom[pid].values())) for pid in pids
    ])

    hypo_cmap = plt.get_cmap('Greens_r')
    hyper_cmap = plt.get_cmap('Reds')
    hypo_norm = Normalize(vmin=-6, vmax=0)
    hyper_norm = Normalize(vmin=0, vmax=6)
    hypo_sm = cm.ScalarMappable(norm=hypo_norm, cmap=hypo_cmap)
    hyper_sm = cm.ScalarMappable(norm=hyper_norm, cmap=hyper_cmap)

    plot_dict = dmr_direction_by_chromosome_pie_array(
        median_delta_all_by_chrom,
        num_dmr_all_by_chrom,
        outer_sm_hyper=hyper_sm,
        outer_sm_hypo=hypo_sm,
    )
    fig = plot_dict['fig']
    fig.savefig(os.path.join(outdir, "dmr_direction_by_chrom_pie_chart_array.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "dmr_direction_by_chrom_pie_chart_array.pdf"), dpi=200)

    plot_dict = dmr_direction_by_chromosome_pie_array(
        median_delta_specific_by_chrom,
        num_dmr_specific_by_chrom,
        outer_sm_hyper=hyper_sm,
        outer_sm_hypo=hypo_sm,
    )
    fig = plot_dict['fig']
    fig.savefig(os.path.join(outdir, "specific_dmr_direction_by_chrom_pie_chart_array.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "specific_dmr_direction_by_chrom_pie_chart_array.pdf"), dpi=200)

    # small bar chart showing the number of DMRs in each chrom / patient for inset (?)
    # by patient
    n_by_patient_by_direction = pd.DataFrame(index=['hypo', 'hyper'], columns=pids)
    n_by_patient_by_direction.loc['hypo'] = [sum(np.array(median_delta_all[pid]) < 0) for pid in pids]
    n_by_patient_by_direction.loc['hyper'] = [sum(np.array(median_delta_all[pid]) > 0) for pid in pids]
    fig = plt.figure(figsize=(4.5, 2.8))
    ax = fig.add_subplot(111)
    bar.stacked_bar_chart(n_by_patient_by_direction, ax=ax, colours=consts.METHYLATION_DIRECTION_COLOURS, ec='k', lw=1., legend=False)
    ax.set_ylabel('Number DMRs')
    ax.set_xlabel('Patient')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_direction_number_by_patient.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "dmr_direction_number_by_patient.pdf"), dpi=200)

    # by chromosome
    n_by_chrom_by_direction = pd.DataFrame(index=['hypo', 'hyper'], columns=chroms, dtype=float)
    n_by_chrom_by_direction.loc['hypo'] = [
        np.mean(
            [sum(np.array(median_delta_all_by_chrom[pid][chrom]) < 0) for pid in pids]
        ) for chrom in chroms
    ]
    n_by_chrom_by_direction.loc['hyper'] = [
        np.mean(
            [sum(np.array(median_delta_all_by_chrom[pid][chrom]) > 0) for pid in pids]
        ) for chrom in chroms
    ]
    fig = plt.figure(figsize=(6., 2.8))
    ax = fig.add_subplot(111)
    bar.stacked_bar_chart(n_by_chrom_by_direction, ax=ax, colours=consts.METHYLATION_DIRECTION_COLOURS, ec='k', lw=1., legend=False)
    ax.set_ylabel('Number DMRs')
    ax.set_xlabel('Chromosome')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_direction_number_by_chrom.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "dmr_direction_number_by_chrom.pdf"), dpi=200)

    # by chromosome background level of DMRs
    n_by_chrom_bg = pd.Series([t.chr for t in clusters.values()]).value_counts()[chroms]
    fig = plt.figure(figsize=(6., 2.8))
    ax = fig.add_subplot(111)
    ax.bar(range(len(chroms)), n_by_chrom_bg, ec='k', lw=1., color='lightgrey')
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms)
    ax.set_ylabel('Number clusters')
    ax.set_xlabel('Chromosome')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_bg_number_by_chrom.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "dmr_bg_number_by_chrom.pdf"), dpi=200)

    # by chromosome as a pct of the BG
    n_by_chrom_by_direction_pct = n_by_chrom_by_direction.divide(n_by_chrom_bg, axis=1) * 100.
    fig = plt.figure(figsize=(6., 2.8))
    ax = fig.add_subplot(111)
    bar.stacked_bar_chart(n_by_chrom_by_direction_pct, ax=ax, colours=consts.METHYLATION_DIRECTION_COLOURS, ec='k', lw=1., legend=False)
    ax.set_ylabel('% of clusters that are DM')
    ax.set_xlabel('Chromosome')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_direction_number_by_chrom_pct.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "dmr_direction_number_by_chrom_pct.pdf"), dpi=200)

    # same but for specific DMRs
    n_by_patient_by_direction_specific = pd.DataFrame(index=['hypo', 'hyper'], columns=pids)
    n_by_patient_by_direction_specific.loc['hypo'] = [sum(np.array(median_delta_specific[pid]) < 0) for pid in pids]
    n_by_patient_by_direction_specific.loc['hyper'] = [sum(np.array(median_delta_specific[pid]) > 0) for pid in pids]
    fig = plt.figure(figsize=(4.5, 2.8))
    ax = fig.add_subplot(111)
    bar.stacked_bar_chart(n_by_patient_by_direction_specific, ax=ax, colours=consts.METHYLATION_DIRECTION_COLOURS, ec='k', lw=1., legend=False)
    ax.set_ylabel('Number DMRs')
    ax.set_xlabel('Patient')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "specific_dmr_direction_number_by_patient.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "specific_dmr_direction_number_by_patient.pdf"), dpi=200)

    # by chromosome
    n_by_chrom_by_direction_specific = pd.DataFrame(index=['hypo', 'hyper'], columns=chroms, dtype=float)
    n_by_chrom_by_direction_specific.loc['hypo'] = [
        np.mean(
            [sum(np.array(median_delta_specific_by_chrom[pid][chrom]) < 0) for pid in pids]
        ) for chrom in chroms
    ]
    n_by_chrom_by_direction_specific.loc['hyper'] = [
        np.mean(
            [sum(np.array(median_delta_specific_by_chrom[pid][chrom]) > 0) for pid in pids]
        ) for chrom in chroms
    ]
    fig = plt.figure(figsize=(6., 2.8))
    ax = fig.add_subplot(111)
    bar.stacked_bar_chart(n_by_chrom_by_direction_specific, ax=ax, colours=consts.METHYLATION_DIRECTION_COLOURS, ec='k', lw=1., legend=False)
    ax.set_ylabel('Number DMRs')
    ax.set_xlabel('Chromosome')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "specific_dmr_direction_number_by_chrom.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "specific_dmr_direction_number_by_chrom.pdf"), dpi=200)

    # by chromosome as a pct of the BG
    n_by_chrom_by_direction_specific_pct = (n_by_chrom_by_direction_specific.divide(n_by_chrom_bg, axis=1) * 100.).loc[:, chroms]
    fig = plt.figure(figsize=(6., 2.8))
    ax = fig.add_subplot(111)
    bar.stacked_bar_chart(n_by_chrom_by_direction_specific_pct, ax=ax, colours=consts.METHYLATION_DIRECTION_COLOURS, ec='k', lw=1., legend=False)
    ax.set_ylabel('% of clusters that are DM')
    ax.set_xlabel('Chromosome')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "specific_dmr_direction_number_by_chrom_pct.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "specific_dmr_direction_number_by_chrom_pct.pdf"), dpi=200)

    # compare against the (normalised) number of oncogenes in each chrom
    oncogene_fn = os.path.join(GIT_LFS_DATA_DIR, 'oncogene_database', 'ongene_human.txt')
    df = pd.read_csv(oncogene_fn, sep='\t', index_col=0, header=0)
    on_chroms = df['Cytoband'].str.replace(r'^(?P<n>[0-9]*)[pq].*', '\g<n>')
    on_chroms = on_chroms.loc[on_chroms.isin(chroms)]
    on_by_chrom = on_chroms.value_counts().loc[chroms]
    on_by_chrom_n = on_by_chrom / pd.Series(chrom_length).loc[chroms] * 1e7

    fig = plt.figure(figsize=(6., 2.8))
    ax = fig.add_subplot(111)
    ax.bar(range(len(chroms)), on_by_chrom_n, facecolor='lightgrey', edgecolor='k', linewidth=1.)
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Oncogene density (a.u.)')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "oncogene_normalised_density_by_chrom.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "oncogene_normalised_density_by_chrom.pdf"), dpi=200)

    unmapped_density, _ = genomics.cg_content_windowed(fa_fn, features=chroms, window_size=window_size, motif='N')

    # get all DMRs (whether significant or not) to use as a null hypothesis
    tmp = get_binned_dmr_locations(
        dmr_res_s1.results,
        clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=False,
        coord_summary_method=coord_summary_method
    )

    # these are the same for all patients, so just take the first
    dmr_loci_all = tmp['dmr_loci'][pids[0]]
    dmr_binned_all = tmp['dmr_binned'][pids[0]]
    dmr_regions_all = tmp['dmr_regions'][pids[0]]

    #1 ) All DMRs
    tmp = get_binned_dmr_locations(
        dmr_res_all,
        clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=True,
        coord_summary_method=coord_summary_method
    )
    dmr_loci_hyper = tmp['dmr_loci_hyper']
    dmr_loci_hypo = tmp['dmr_loci_hypo']
    dmr_regions_hyper = tmp['dmr_regions_hyper']
    dmr_regions_hypo = tmp['dmr_regions_hypo']


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
        clusters,
        chrom_length,
        window_size=window_size,
        split_by_direction=True,
        coord_summary_method=coord_summary_method
    )

    dmr_loci_hyper_specific = tmp['dmr_loci_hyper']
    dmr_loci_hypo_specific = tmp['dmr_loci_hypo']
    dmr_regions_hyper_specific = tmp['dmr_regions_hyper']
    dmr_regions_hypo_specific = tmp['dmr_regions_hypo']

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

    plt.close('all')
    plt.pause(1)

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
                ks_subwindows_hyper[pid][chrom][(s, e)] = None
                ks_subwindows_hypo[pid][chrom][(s, e)] = None

                sub_all = this_all[(this_all >= s) & (this_all < e)]
                sub_hyper = this_hyper[(this_hyper >= s) & (this_hyper < e)]
                sub_hypo = this_hypo[(this_hypo >= s) & (this_hypo < e)]
                if len(sub_hypo) > 0:
                    jobs[(pid, chrom, (s, e), 'hypo')] = pool.apply_async(stats.ks_2samp, args=(sub_hypo, sub_all))
                if len(sub_hyper) > 0:
                    jobs[(pid, chrom, (s, e), 'hyper')] = pool.apply_async(stats.ks_2samp, args=(sub_hyper, sub_all))
                s += roll_dw
                e = s + subwindow_size

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
            ks_subwindows_hyper[pid][chrom] = pd.Series(ks_subwindows_hyper[pid][chrom], dtype=float)
            ks_subwindows_hypo[pid][chrom] = pd.Series(ks_subwindows_hypo[pid][chrom], dtype=float)

    # overlay DMRs with features?
    # use utils.genomics.GtfAnnotation (.region method)
    # or just lookup in the annotation file?


    # plot single chromosomes in individual patients to illustrate examples of interest.
    pid_chroms = [
        ('018', '2'),
        ('018', '5'),
        ('018', '6'),
        ('018', '14'),
        ('019', '2'),
        ('019', '5'),
        ('019', '6'),
        ('019', '14'),
        ('018', '3'),
        ('017', '3'),
    ]

    for pid, chrom in pid_chroms:

        plt_dict = plot_one_chrom_location(
            dmr_loci_hyper[pid][chrom],
            dmr_loci_hypo[pid][chrom],
            dmr_binned_all[chrom],
            chrom_length[chrom],
            pval_hyper=ks_subwindows_hyper[pid][chrom],
            pval_hypo=ks_subwindows_hypo[pid][chrom],
            unmapped_density=unmapped_density[chrom] / window_size * 100.,
            cmap=new_cmap
        )
        fig = plt_dict['fig']
        fig.savefig(os.path.join(outdir, "patient_%s_chrom_%s_locations.png" % (pid, chrom)), dpi=200)

    plt.close('all')
    plt.pause(1)

    # lookup intersecting features in selected DMRs (and the background - all clusters)
    gtf_fn = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh37', 'gtf', 'Homo_sapiens.GRCh37.87.gtf.gz')

    gene_biotypes = {'protein_coding'}

    hash_dict = {}
    hash_dict.update(dmr_hash_dict)
    hash_dict['gene_biotypes'] = tuple(sorted(gene_biotypes))
    hash_dict['gtf_fn'] = gtf_fn.replace(LOCAL_DATA_DIR, '')
    the_hash = tsgd.dmr_results_hash(
        me_obj.meta.index.tolist(),
        dmr_hash_dict,
        gene_biotypes=tuple(sorted(gene_biotypes)),
        gtf_fn=gtf_fn.replace(LOCAL_DATA_DIR, '')
    )

    lookup_fn = os.path.join(DMR_LOAD_DIR, 'dmr_gtf_lookup.%s.pkl' % the_hash)
    if os.path.isfile(lookup_fn):
        logger.info("Loading pre-computed DMR / GTF lookup data from file %s.", lookup_fn)
        with open(lookup_fn, 'rb') as f:
            tmp = pickle.load(f)
            lookup_all = tmp['region']
            class_all = tmp['classification']
            cluster_id_to_region = tmp['cluster_id_to_region']
            gtf_fn = os.path.join(LOCAL_DATA_DIR, tmp['gtf_fn'])
            db = genomics.GtfAnnotation(gtf_fn)
    else:
        logger.info("No pre-computed DMR / GTF lookup found. Generating results now (this can take a long time).")
        db = genomics.GtfAnnotation(gtf_fn, logger=None)
        lookup_all, class_all, cluster_id_to_region = classify_dmr_regions(db, clusters, gene_biotypes=gene_biotypes)
        # dump to pickle
        logger.info("Dumping pre-computed DMR / GTF lookup data to file %s.", lookup_fn)
        with open(lookup_fn, 'wb') as f:
            pickle.dump(
                {
                    'region': lookup_all,
                    'classification': class_all,
                    'gtf_fn': gtf_fn.replace(LOCAL_DATA_DIR, ''),
                    'cluster_id_to_region': cluster_id_to_region
                },
                f
            )

    region_to_cluster_id = dict([x[::-1] for x in cluster_id_to_region.items()])

    # split into chromosomes for easier access
    tmp_class_all = collections.defaultdict(dict)
    tmp_lookup_all = collections.defaultdict(dict)
    for cid in class_all.keys():
        reg = cluster_id_to_region[cid]
        tmp_class_all[reg[0]][cid] = class_all[cid]
        tmp_lookup_all[reg[0]][cid] = lookup_all[cid]

    class_all = tmp_class_all
    lookup_all = tmp_lookup_all

    # now use this long list to populate the short lists
    class_hyper = {}
    class_hypo = {}
    class_patients = {}

    for pid in pids:
        class_hyper[pid] = {}
        class_hypo[pid] = {}
        class_patients[pid] = {}
        for chrom in chroms:
            class_hyper[pid][chrom] = {}
            class_hypo[pid][chrom] = {}
            class_patients[pid][chrom] = {}
            for reg in dmr_regions_hyper[pid][chrom]:
                k = (chrom,) + reg
                cid = region_to_cluster_id[k]
                class_hyper[pid][chrom][cid] = class_all[chrom][cid]
                class_patients[pid][chrom][cid] = class_all[chrom][cid]
            for reg in dmr_regions_hypo[pid][chrom]:
                k = (chrom,) + reg
                cid = region_to_cluster_id[k]
                class_hypo[pid][chrom][cid] = class_all[chrom][cid]
                class_patients[pid][chrom][cid] = class_all[chrom][cid]

    # repeat for specific DMRs
    class_hyper_specific = {}
    class_hypo_specific = {}
    class_patients_specific = {}

    for pid in pids:
        class_hyper_specific[pid] = {}
        class_hypo_specific[pid] = {}
        class_patients_specific[pid] = {}
        for chrom in chroms:
            class_hyper_specific[pid][chrom] = {}
            class_hypo_specific[pid][chrom] = {}
            class_patients_specific[pid][chrom] = {}
            for reg in dmr_regions_hyper_specific[pid][chrom]:
                k = (chrom,) + reg
                cid = region_to_cluster_id[k]
                class_hyper_specific[pid][chrom][cid] = class_all[chrom][cid]
                class_patients_specific[pid][chrom][cid] = class_all[chrom][cid]
            for reg in dmr_regions_hypo_specific[pid][chrom]:
                k = (chrom,) + reg
                cid = region_to_cluster_id[k]
                class_hypo_specific[pid][chrom][cid] = class_all[chrom][cid]
                class_patients_specific[pid][chrom][cid] = class_all[chrom][cid]

    # analyse the distributions
    lookup_cats = ['intergenic', 'intron', 'exon', 'first_exon']
    cat_labels = dict([
        (cat, cat.capitalize().replace('_', ' ')) for cat in lookup_cats
    ])

    colours_by_lookup = {
        'intergenic': '#abdda4',
        'intron': '#2b83ba',
        'exon': '#fdae61',
        'first_exon': '#d7191c',
    }

    def dmr_class_to_counts(dat):
        res = {}
        for k1, v1 in dat.items():
            res[k1] = {}
            for k2, v2 in v1.items():
                res[k1][k2] = pd.Series(v2).value_counts().reindex(lookup_cats).fillna(0)
            res[k1] = pd.DataFrame(res[k1])
        return res

    # reduce to counts
    class_counts_hyper = dmr_class_to_counts(class_hyper)
    class_counts_hypo = dmr_class_to_counts(class_hypo)
    class_counts_hyper_specific = dmr_class_to_counts(class_hyper_specific)
    class_counts_hypo_specific = dmr_class_to_counts(class_hypo_specific)

    # pairwise distances
    pdist_hyper = pairwise_distance_by_class_count(class_counts_hyper)
    pdist_hypo = pairwise_distance_by_class_count(class_counts_hypo)
    pdist_hyper_specific = pairwise_distance_by_class_count(class_counts_hyper_specific)
    pdist_hypo_specific = pairwise_distance_by_class_count(class_counts_hypo_specific)

    # convert to vector form and compute linkage
    # linkage_method = 'average'
    # pdist_hyper_lkg = [
    #     cluster.hierarchy.linkage(
    #         cluster.hierarchy.distance.squareform(x), method=linkage_method
    #     )  for x in pdist_hyper
    # ]
    # pdist_hypo_lkg = [
    #     cluster.hierarchy.linkage(
    #         cluster.hierarchy.distance.squareform(x), method=linkage_method
    #     )  for x in pdist_hypo
    # ]
    # pdist_hyper_specific_lkg = [
    #     cluster.hierarchy.linkage(
    #         cluster.hierarchy.distance.squareform(x), method=linkage_method
    #     )  for x in pdist_hyper_specific
    # ]
    # pdist_hypo_specific_lkg = [
    #     cluster.hierarchy.linkage(
    #         cluster.hierarchy.distance.squareform(x), method=linkage_method
    #     )  for x in pdist_hypo_specific
    # ]

    fig, axs = pie_chart_location(
        class_all,
        class_patients,
        lookup_cats,
        colours_by_lookup,
        cat_labels,
    )

    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_both.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_both.pdf'), dpi=200)

    fig, axs = pie_chart_location(
        class_all,
        class_hyper,
        lookup_cats,
        colours_by_lookup,
        cat_labels,
    )

    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_hyper.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_hyper.pdf'), dpi=200)

    fig, axs = pie_chart_location(
        class_all,
        class_hypo,
        lookup_cats,
        colours_by_lookup,
        cat_labels,
    )

    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_hypo.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_hypo.pdf'), dpi=200)

    fig, axs = pie_chart_location(
        class_all,
        class_patients_specific,
        lookup_cats,
        colours_by_lookup,
        cat_labels,
    )

    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_both_specific.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_both_specific.pdf'), dpi=200)

    fig, axs = pie_chart_location(
        class_all,
        class_hyper_specific,
        lookup_cats,
        colours_by_lookup,
        cat_labels,
    )

    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_hyper_specific.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_hyper_specific.pdf'), dpi=200)

    fig, axs = pie_chart_location(
        class_all,
        class_hypo_specific,
        lookup_cats,
        colours_by_lookup,
        cat_labels,
    )

    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_hypo_specific.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'dmr_classification_pie_array_hypo_specific.pdf'), dpi=200)
