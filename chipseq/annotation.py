import pandas as pd
import numpy as np
import multiprocessing as mp
import tabix
from chipseq import feature_enrichment
from chipseq import loader
from settings import LOCAL_DATA_DIR
import collections
import os
import subprocess
from matplotlib import pyplot as plt
import seaborn as sns
from utils import output

SOURCES = {'ensembl', 'havana', 'ensembl_havana'}


def prepare_tabix_indexed_gtf(gtf_fn):
    """
    It is necessary to run this procedure to index a GTF file for use with the tabix application.
    This gtf must be sorted first
    The GTF must first be SORTED, then BGZIPPED, then INDEXED
    """
    filestem, ext = os.path.splitext(gtf_fn)
    if ext.lower() == '.gz':
        cmd = "zcat {fn}"
        outfile = filestem + ".bgz"
    else:
        cmd = "cat {fn}"
        outfile = gtf_fn + ".bgz"
    cmd += " | bedtools sort | bgzip > {outfn}"
    cmd = cmd.format(fn=gtf_fn, outfn=outfile)
    subprocess.call(cmd, shell=True)
    cmd = "tabix {fn}".format(fn=outfile)
    subprocess.call(cmd, shell=True)


def assign_peaks_to_basic_features(peak_dat, gtf_fn, tss_pad=500):
    """
    Given the peak data (from MACS), get basic peak assignment: TSS, exon, intron, intergenic
    (in that order of priority)
    :param peak_dat: Pandas DataFrame containing ChIP peaks of interest. Must have the columns start, end and chrom.
    :param gtf_fn: Path to a sorted, BGzipped, tabix-indexed GTF file.
    :return:
    """
    tb = tabix.open(gtf_fn)
    peak_assignment = pd.Series(index=peak_dat.index)
    assignment_count = collections.defaultdict(float)
    for i, row in peak_dat.iterrows():
        qry = tb.query(str(row.chrom), row.start, row.end)
        hit_gene = False
        hit_exon = False
        hit_tss = False
        for t in qry:
            _, typ, ftr, i0, i1, _, strand, _, attr = t
            if typ in SOURCES:
                if ftr in {'gene', 'transcript'}:
                    hit_gene = True
                if ftr == 'exon':
                    if 'exon_number "1"' in attr:
                        tss_loc = int(i0) if strand == '+' else int(i1)
                        if (row.start - tss_pad) <= tss_loc <= (row.end + tss_pad):
                            hit_tss = True
                        else:
                            hit_exon = True
                    else:
                        hit_exon = True
        if hit_tss:
            assignment_count['tss'] += 1.
            peak_assignment[i] = 'tss'
        elif hit_exon:
            assignment_count['exon'] += 1.
            peak_assignment[i] = 'exon'
        elif hit_gene:
            assignment_count['intron'] += 1.
            peak_assignment[i] = 'intron'
        else:
            assignment_count['intergenic'] += 1.
            peak_assignment[i] = 'intergenic'
    return peak_assignment


def get_closest_tss(peak_loc, tss_loc, tss_names):
    """
    For every peak in peak_loc, find the nearest TSS in tss_loc. Provide the name (from tss_names) and distance.
    Distance is defined so that negative values are upstream of the TSS.
    NB not considering chromosomes here - so split by those in an outer loop.
    :param peak_loc: Numpy array containing the coordinates of the peaks (all in the same chromosome).
    :param tss_loc: Numpy array containing the coordinates of each TSS (all in the matching chromosome).
    :param tss_names: Numpy array containing the name(s) of each TSS (all in the matching chromosome).
    :return: Pandas DataFrame with two columns: 'gene' and 'distance_to_tss'.
    """
    N_pk = len(peak_loc)
    N = len(tss_loc)
    ix = np.searchsorted(tss_loc, peak_loc, side='left')

    closest_tss = pd.DataFrame(index=range(N_pk), columns=['gene', 'distance_to_tss'])

    # ix == 0 and ix == N are special cases (closest TSS is unambiguous)
    closest_tss.loc[ix == 0, 'gene'] = tss_names[0]
    closest_tss.loc[ix == N, 'gene'] = tss_names[N - 1]
    closest_tss.loc[ix == 0, 'distance_to_tss'] = peak_loc[ix == 0] - tss_loc[0]
    closest_tss.loc[ix == N, 'distance_to_tss'] = peak_loc[ix == N] - tss_loc[N - 1]

    # everything else: we need to check "one up and one down" for the closest match
    for i in range(1, N):
        n_in_segment = (ix == i).sum()
        dist_down = peak_loc[ix == i] - tss_loc[i - 1]
        dist_up = tss_loc[i] - peak_loc[ix == i]
        winning_ix_rel = np.array([0, 1])[(dist_up < dist_down).astype(int)]
        # winning_ix_abs = np.array([i - 1, i])[(dist_up < dist_down).astype(int)]
        winning_ix_abs = winning_ix_rel + i - 1
        closest_tss.loc[ix == i, 'gene'] = tss_names[winning_ix_abs]
        closest_tss.loc[ix == i, 'distance_to_tss'] = np.vstack((dist_down, -dist_up)).transpose()[
            range(n_in_segment),
            winning_ix_rel
        ]

    return closest_tss



if __name__ == '__main__':
    # tss_pad controls how far from a TSS is still counted as TSS
    tss_pad = 500
    run_types = ['default', 'broad']
    chip_target_colours = {
        'BMI1': '#1b9e77',
        'CHD7': '#d95f02',
        'H3K27me3': '#7570b3',
        'H3K27ac': '#e7298a',
        'H3K36me3': '#66a61e',
        'H3K4me3': '#e6ab02',
    }
    outdir = output.unique_output_dir("chipseq_peak_annotation", reuse_empty=True)
    # in order to use tabix (fast genomic region indexer / looker-upper), we need a sorted, bgzipped GTF file
    # if this doesn't exist, use prepare_tabix_indexed_gtf on a regular GTF
    fn_bzip = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human',
        'ensembl',
        'GRCh38.release87',
        'gtf',
        'Homo_sapiens.GRCh38.87.sorted.gtf.bgz'
    )
    fn_gzip = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human',
        'ensembl',
        'GRCh38.release87',
        'gtf',
        'Homo_sapiens.GRCh38.87.gtf.gz'
    )

    reg, names = feature_enrichment.get_gene_tss_from_gtf(fn_gzip, distance=0)

    peak_assignment = {}
    peak_sizes = {}
    peak_data = {}
    closest_tss = {}

    for run_type in run_types:

        # recreate the essential operation of Homer suite's annotatePeaks.pl
        peaks = loader.load_macs2_by_patient('all', run_type=run_type)

        # 1) overlap of peaks with simple genomic features (intron, exon, tss, intergenic)
        # define TSS as being exon 1 start +/- a padding distance
        # order of priority:
        # 1) TSS (or close to a TSS)
        # 2) Exon
        # 3) Intron
        # 4) Intergenic (none of the above)

        pac = {}
        ps = {}

        pool = mp.Pool()
        jobs = {}

        for k in peaks.data:
            print "Sample %s" % k

            this_assignment = collections.defaultdict(float)
            # add a new column to store the assignment
            peaks.data[k].insert(0, 'assignment', None)
            dat = peaks.data[k]

            ps[k] = (dat.end - dat.start).values

            for c in dat.chrom.unique():
                this_dat = dat.loc[dat.chrom == c]
                jobs[(k, c)] = pool.apply_async(assign_peaks_to_basic_features, args=(this_dat, fn_bzip), kwds=dict(tss_pad=tss_pad))

        pool.close()
        pool.join()

        for k, c in jobs:
            this_res = jobs[(k, c)].get()
            peaks.data[k].loc[this_res.index, 'assignment'] = this_res.values

            # get the number of peaks assigned to each class
            this_counts = this_res.groupby(this_res).size()

            pac.setdefault(k, collections.defaultdict(float))
            for d in this_counts.index:
                pac[k][d] += this_counts[d]

        pac = pd.DataFrame.from_dict(pac)
        pac.to_csv(os.path.join(outdir, "%s_peak_assignment_count.csv" % run_type))

        # peaks.data[k].to_csv(os.path.join(outdir, "%s_%s_annotated_peaks.csv" % (run_type, k)))

        peak_assignment[run_type] = pac
        peak_sizes[run_type] = ps
        peak_data[run_type] = peaks

        # PLOT: distribution of peak sizes by ChIP target
        chip_groups = [
            ('H3K27ac',),
            ('H3K36me3',),
            ('H3K4me3',),
            ('BMI1', 'CHD7',),
            ('H3K27me3',),
        ]

        for grp in chip_groups:
            snames_all = []
            x_all = []
            n = sum([(peaks.meta.chip_target == ct).sum() for ct in grp])
            fig = plt.figure(figsize=(5.2, 7.2 / 8. * n))
            ax = fig.add_subplot(111)
            i = 0
            for ct in grp:
                ix = peaks.meta.chip_target == ct
                snames = peaks.meta.index[ix]

                this_dat = [ps[sn] for sn in snames]
                this_log_dat = [np.log10(t) for t in this_dat]
                this_x = range(i, i + len(this_dat))

                i = len(this_dat)
                x_all.extend(this_x)
                snames_all.extend(snames)

                vp = ax.violinplot(this_log_dat, this_x, widths=0.8, vert=False)
                ax.plot(np.log10(np.array([np.median(t) for t in this_dat])), this_x, 'ko')
                for pc in vp['bodies']:
                    pc.set_facecolor(chip_target_colours[ct])
                vp['cbars'].set_edgecolor(chip_target_colours[ct])
                vp['cmins'].set_edgecolor(chip_target_colours[ct])
                vp['cmaxes'].set_edgecolor(chip_target_colours[ct])
            ax.set_yticks(x_all)
            ax.set_yticklabels(snames_all)
            ax.set_xlim([2, 6])
            ax.set_xlabel('log10(peak size)')
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, "%s_peak_size_distribution _%s.png" % (run_type, '_'.join(grp))), dpi=200)

        # PLOT: dual bar chart - number of peaks and assignment %

        # reorder by chip target
        idx = peaks.meta.sort_values(['chip_target', 'input_sample']).index

        ss = pac.sum().loc[idx]
        pct = pac.divide(ss, axis=1)[idx.tolist()] * 100.
        ss_colours = pd.Series(index=ss.index)
        for ct in peaks.meta.chip_target.unique():
            ss_colours[peaks.meta.chip_target == ct] = chip_target_colours[ct]

        # plot the numbers of peaks and the assignment distribution
        fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(7, 8))

        ss.plot(kind='barh', ax=axs[0], width=0.9, color=ss_colours.values)
        pct.transpose().plot(kind='barh', stacked=True, ax=axs[1], width=0.9)

        box = axs[0].get_position()
        axs[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])

        box = axs[1].get_position()
        axs[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # axs[0].set_ylabel('Number of ChIP peaks')
        axs[0].set_xlim([0, 150000])
        axs[0].set_xlabel('Number of ChIP peaks')
        # axs[1].set_ylabel('Percentage assigned to feature')
        axs[1].set_xlabel('% assigned to feature')
        axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # fig.subplots_adjust(right=0.85, hspace=0.05, top=0.95, bottom=0.28)
        fig.subplots_adjust(left=0.28, right=0.8, wspace=0.2, top=0.95, bottom=0.1)

        fig.savefig(os.path.join(outdir, "%s_peaks_feature_assignment.png" % run_type), dpi=200)

        pool = mp.Pool()
        jobs = {}

        for k in peaks.data:
            print "Sample %s" % k

            for c in dat.chrom.unique():
                c = str(c)
                this_dat = dat.loc[dat.chrom.astype(str) == c]
                if run_type == 'default':
                    this_peak_loc = this_dat.start + this_dat.rel_peak_pos
                else:
                    this_peak_loc = 0.5 * (this_dat.start + this_dat.end)
                this_peak_loc = this_peak_loc.values

                this_tss_loc = np.array([t[1] for t in reg if t[0] == c])
                this_tss_genes = np.array([';'.join(names[(c, t)]) for t in this_tss_loc])

                jobs[(k, c)] = pool.apply_async(get_closest_tss, args=(
                    this_peak_loc,
                    this_tss_loc,
                    this_tss_genes
                ))

        pool.close()
        pool.join()

        ct_temp = {}
        closest_tss_chroms = {}
        for k, c in jobs:
            res = jobs[(k, c)].get()
            ct_temp.setdefault(k, []).append(res)
            closest_tss_chroms.setdefault(k, []).extend([c] * res.shape[0])

        closest_tss[run_type] = {}
        for k in peaks.data:
            closest_tss[run_type] = pd.concat(
                ct_temp[k], axis=0
            )
            closest_tss[run_type].reset_index(inplace=True)
            closest_tss[run_type].insert(0, 'chrom', closest_tss_chroms[k])
