import pandas as pd
import os
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import collections

from methylation import loader
from utils import setops, genomics, output, excel
from plotting import common
from settings import LOCAL_DATA_DIR
from scripts.hgic_final import consts


def polar_distribution_plot_with_kdes(
    logfc_vals,
    chrom_length,
    bg_density=None,
    unmapped_density=None,
    window_size=2e5,
    inner_r=3.,
    central_width = 0.25,
    gap_radians = 2 * np.pi / 200.,
    bar_scaling=0.03,
    bar_width_bp=1000,
    bg_fmin=0.,
    bg_fmax=0.99,
    bg_vmin=None,
    bg_vmax=None,
    colour=None,
    hyper_colour=consts.METHYLATION_DIRECTION_COLOURS['hyper'],
    hypo_colour = consts.METHYLATION_DIRECTION_COLOURS['hypo'],
    unmap_threshold_pct=10.,
    plot_border = True,
    bg_cmap=plt.cm.get_cmap('RdYlBu_r'),
    ax=None,
):
    chroms = chrom_length.keys()

    if colour is not None:
        hyper_colour = colour
        hypo_colour = colour

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

    new_ax = False
    if ax is None:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='polar')
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location('N')
        # hack required to convert patches correctly to polar coords
        # https://github.com/matplotlib/matplotlib/issues/8521
        ax.bar(0, 1).remove()
        new_ax = True

    text_handles = {}
    border_patches = {}
    bg_handles = {}
    unmapped_handles = {}

    curr_theta = 0
    for chrom in chroms:
        th_set = False
        this_vals = logfc_vals[chrom].astype(float)
        this_up = this_vals[this_vals > 0]
        this_down = this_vals[this_vals < 0]

        if bg_density is not None:
            this_bg = bg_density[chrom]

            th = curr_theta + np.array(this_bg.index.tolist() + [chrom_length[chrom]]) * radians_per_bp
            tt = np.array([th, th])
            rr = np.zeros_like(tt) + inner_r
            rr[1] = outer_r

            cc = np.array([this_bg.values])

            norm = common.MidpointNormalize(midpoint=bg_mean, vmin=bg_vmin, vmax=bg_vmax)
            h_bg = ax.pcolor(tt, rr, cc, cmap=bg_cmap, norm=norm, vmin=bg_vmin, vmax=bg_vmax)
            bg_handles[chrom] = h_bg
            th_set = True

        if unmapped_density is not None:
            this_unmapped = unmapped_density[chrom]
            this_unmapped_pct = this_unmapped / float(window_size) * 100.
            th = curr_theta + np.array(this_unmapped.index.tolist() + [chrom_length[chrom]]) * radians_per_bp

            uu = np.ma.masked_less(np.array([this_unmapped_pct.values]), unmap_threshold_pct)
            # since we don't care about the extent of unmapping, replace all values with a single one
            uu[~uu.mask] = 0.5
            h_unmapped = ax.pcolor(tt, rr, uu, cmap='Greys', vmax=1., vmin=0.)
            unmapped_handles[chrom] = h_unmapped
            th_set = True

        if not th_set:
            th = [curr_theta + chrom_length[chrom] * radians_per_bp]

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
            # outer_r + kde_width / 2.,
            outer_r + bar_scaling * 5.,
            chrom,
            horizontalalignment='center',
            verticalalignment='center'
        )
        text_handles[chrom] = h_text

        # plot the values
        th_up = curr_theta + np.array(this_up.index) * radians_per_bp
        th_down = curr_theta + np.array(this_down.index) * radians_per_bp

        ax.bar(
            th_up,
            this_up.values * bar_scaling,
            width=radians_per_bp * bar_width_bp,
            bottom=outer_r,
            color=hyper_colour,
            )
        ax.bar(
            th_down,
            this_down.values * bar_scaling,
            width=radians_per_bp * bar_width_bp,
            bottom=inner_r,
            color=hypo_colour,
            )
        # ax.scatter(th_up, outer_r + this_up.values, color=hyper_colour)
        # ax.scatter(th_down, inner_r + this_down.values, color=hypo_colour)

        ax.set_facecolor('w')

        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        curr_theta = th[-1] + gap_radians

    if new_ax:
        fig.tight_layout()

    return {
        'fig': ax.figure,
        'ax': ax,
        'border_patches': border_patches,
        'text_handles': text_handles,
        'unmapped_handles': unmapped_handles,
        'bg_handles': bg_handles,
    }


if __name__ == "__main__":
    outdir = output.unique_output_dir()
    anno = loader.load_illumina_methylationepic_annotation()
    indir = '/home/gabriel/Dropbox/research/qmul/MB_project/current/core_pipeline/methylation/'

    # first, add annotation to results and save out again
    file_list = [
        'dmps_1299_noob.xlsx',
        'dmps_3021_noob.xlsx',
        'dmps_1299_swan.xlsx',
        'dmps_3021_swan.xlsx',
    ]
    for fn in file_list:
        ff = os.path.join(indir, fn)
        dat = pd.read_excel(ff, header=0, index_col=0, sheet_name=None)
        new_dat = {}
        for k, df in dat.items():
            to_add = anno.reindex(df.index)[['CHR', 'MAPINFO', 'UCSC_RefGene_Name', 'UCSC_RefGene_Group']]
            # to_add.columns = ['chrom', 'coord', 'genes', 'gene_relation']
            df.insert(df.shape[1], 'chrom', to_add.CHR)
            df.insert(df.shape[1], 'coord', to_add.MAPINFO.fillna(-1).astype(int))
            df.insert(df.shape[1], 'gene', [','.join(t) if hasattr(t, '__iter__') else '' for t in to_add.UCSC_RefGene_Name])
            df.insert(df.shape[1], 'gene_relation', [','.join(t) if hasattr(t, '__iter__') else '' for t in to_add.UCSC_RefGene_Group])
            new_dat[k] = df
        excel.pandas_to_excel(new_dat, os.path.join(outdir, fn.replace('.xlsx', '.annotated.xlsx')))

    dmp_fn = os.path.join(indir, 'dmps_3021_swan.xlsx')
    dmps = pd.read_excel(dmp_fn, header=0, index_col=0, sheet_name=None)

    # combine all DMPs into a single wideform
    cols = reduce(lambda x, y: x + y, [['%s' % t, '%s_logFC' % t, '%s_FDR' % t] for t in dmps])
    all_probes = setops.reduce_union(*[v.loc[v['adj.P.Val'] < 0.05].index for v in dmps.values()])
    all_probes = all_probes.intersection(anno.index)

    dmps_all = pd.DataFrame(index=all_probes, columns=['CHR', 'coord', 'genes'] + cols)
    dmps_all.loc[:, 'CHR'] = anno.loc[dmps_all.index, 'CHR']
    dmps_all.loc[:, 'coord'] = anno.loc[dmps_all.index, 'MAPINFO']
    dmps_all.loc[:, 'genes'] = anno.loc[dmps_all.index, 'UCSC_RefGene_Name']
    dmps_all.loc[:, dmps.keys()] = False

    for k, v in dmps.items():
        this = v.loc[v['adj.P.Val'] < 0.05]
        this = this.loc[this.index.intersection(all_probes)]
        dmps_all.loc[this.index, k] = True
        dmps_all.loc[this.index, "%s_logFC" % k] = this['logFC']
        dmps_all.loc[this.index, "%s_FDR" % k] = this['adj.P.Val']

    # reorder
    dmps_all = dmps_all.loc[anno.loc[all_probes].sort_values(by=['CHR', 'MAPINFO']).index]

    # polar plot
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

    # get probe density
    bg_density = {}
    for ch in chroms:
        this_coord = anno.loc[anno.CHR == ch].MAPINFO
        tt = np.histogram(anno.loc[anno.CHR == ch].MAPINFO, np.arange(0, chrom_length[ch], window_size))
        bg_density[ch] = pd.Series(tt[0], index=tt[1][:-1])

    # full polar plot
    # plot_dict = polar_distribution_plot_with_kdes(
    #         dmps_all,
    #         chrom_length,
    #         bg_density=None,

    # plot
    for_plot = {}
    for k in dmps.keys():
        ix = dmps_all[k].values
        this_logfc = dmps_all.loc[ix, '%s_logFC' % k]
        this_chr = dmps_all.loc[ix, 'CHR']
        this_coord = dmps_all.loc[ix, 'coord']
        logfc_vals = {}
        for ch in chroms:
            tt = this_logfc.loc[this_chr == ch]
            tt.index = this_coord.loc[this_chr == ch]
            logfc_vals[ch] = tt
        for_plot[k] = logfc_vals

    for cond in for_plot:
        plot_dict = polar_distribution_plot_with_kdes(
            for_plot[cond],
            chrom_length,
            bg_density=bg_density,
            bar_width_bp=3e6,
            bar_scaling=0.1
        )
        plot_dict['fig'].savefig(os.path.join(outdir, "%s_dmps_genomewide.png" % cond.replace(' ', '')), dpi=200)

    gs = plt.GridSpec(nrows=3, ncols=1, height_ratios=[3, 1, 3], hspace=0.01, left=0.05, right=0.97, bottom=0.15, top=0.99)

    cond_colours = {
        'shBMI1 - Scr': 'r',
        'shCHD7 - Scr': 'b',
        'shBMI1shCHD7 - shCHD7': 'k',
    }
    width_frac = 1 / 500.

    for chrom in chroms:
        width = chrom_length[chrom] * width_frac
        this_bg = bg_density[chrom]
        xx = this_bg.index.tolist() + [chrom_length[chrom]]
        xx = np.array([xx, xx])
        yy = np.zeros_like(xx)
        yy[1] = 1.
        cc = np.array([this_bg.values])

        fig = plt.figure(figsize=(12, 3.5))
        ax_ctr = fig.add_subplot(gs[1])
        axs = [
            fig.add_subplot(gs[0], sharex=ax_ctr),
            ax_ctr,
            fig.add_subplot(gs[2], sharex=ax_ctr),
        ]
        ax_ctr.pcolor(xx, yy, cc, cmap='RdYlBu_r')

        for i, (cnd, col) in enumerate(cond_colours.items()):
            a = for_plot[cnd][chrom]
            hyper = a[a > 0]
            hypo = a[a < 0]
            axs[0].bar(hyper.index, hyper.values, bottom=0, width=width, color=col, alpha=0.6 if i > 0 else 1, zorder=10+i)
            axs[2].bar(hypo.index, hypo.values, bottom=0, width=width, color=col, alpha=0.6 if i > 0 else 1, zorder=10+i)

        # axs[1].pcolor(xx, yy, cc, cmap='RdYlBu_r')
        # axs[0].bar(hyper1.index, hyper1.values, bottom=0, width=width, color='r', alpha=0.6, zorder=10)
        # axs[0].bar(hyper2.index, hyper2.values, bottom=0, width=width, color='k', zorder=9)
        # axs[2].bar(hypo1.index, hypo1.values, bottom=0, width=width, color='r', alpha=0.6, zorder=10)
        # axs[2].bar(hypo2.index, hypo2.values, bottom=0, width=width, color='k', zorder=9)

        axs[0].xaxis.set_visible(False)
        axs[1].xaxis.set_visible(False)
        axs[1].yaxis.set_visible(False)
        axs[2].xaxis.set_ticks([chrom_length[chrom]])
        axs[0].set_ylim([0, 1.1 * max(hyper.max(), hypo.abs().max())])

        axs[1].set_xlim([0, chrom_length[chrom]])
        fig.savefig(os.path.join(outdir, "chrom%s.png" % chrom), dpi=200)



