from chipseq import loader, annotation, feature_enrichment
from settings import LOCAL_DATA_DIR
from utils import output
import os
import collections
import pandas as pd
import numpy as np
import multiprocessing as mp

SOURCES = {'ensembl', 'havana', 'ensembl_havana'}


if __name__ == '__main__':
    # matplotlib import: check whether we have an X session
    if 'DISPLAY' in os.environ:
        from matplotlib import pyplot as plt
    else:
        print "No $DISPLAY env var, so we'll use Agg as the matplotlib backend."
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
    import seaborn as sns

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

        out_subdir = os.path.join(outdir, run_type)
        os.makedirs(out_subdir)

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
                jobs[(k, c)] = pool.apply_async(
                    annotation.assign_peaks_to_basic_features,
                    args=(this_dat, fn_bzip),
                    kwds=dict(tss_pad=tss_pad)
                )

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
        pac.to_csv(os.path.join(out_subdir, "peak_assignment_count.csv"))

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
            fig.savefig(os.path.join(out_subdir, "peak_size_distribution_%s.png" % '_'.join(grp)), dpi=200)

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

        fig.savefig(os.path.join(out_subdir, "peaks_feature_assignment.png"), dpi=200)

        pool = mp.Pool()
        jobs = {}

        # get the closest TSS

        for k in peaks.data:
            print "Sample %s" % k
            dat = peaks.data[k]

            # add new columns to store the information
            dat.insert(0, 'distance_to_tss', None)
            dat.insert(0, 'closest_tss', None)

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

                jobs[(k, c)] = pool.apply_async(
                    annotation.get_closest_tss,
                    args=(this_peak_loc, this_tss_loc, this_tss_genes)
                )

        pool.close()
        pool.join()

        ct_temp = {}
        closest_tss_chroms = {}
        # insert the results into the original df in chunks
        for k, c in jobs:
            dat = peaks.data[k]
            res = jobs[(k, c)].get()
            # ensure that this is the SAME indexing that we used earlier to generate results
            dat.loc[dat.chrom.astype(str) == c, ['closest_tss', 'distance_to_tss']] = res.values

        # save data
        for k in peaks.data:
            peaks.data[k].to_csv(os.path.join(out_subdir, "%s_annotated_peaks.csv" % k))