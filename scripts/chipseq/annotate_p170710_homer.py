from chipseq import loader, annotation, feature_enrichment
from settings import LOCAL_DATA_DIR
from utils import output
import os
import pandas as pd
import numpy as np
from rnaseq import general


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

    sources = {'ensembl', 'havana', 'ensembl_havana'}

    run_types = ['homer_default', 'homer_broad']
    chip_target_colours = {
        'BMI1': '#1b9e77',
        'CHD7': '#d95f02',
        'H3K27me3': '#7570b3',
        'H3K27ac': '#e7298a',
        'H3K36me3': '#66a61e',
        'H3K4me3': '#e6ab02',
    }
    chip_groups = [
        ('H3K27ac',),
        ('H3K36me3',),
        ('H3K4me3',),
        ('BMI1', 'CHD7',),
        ('H3K27me3',),
    ]
    annot_types = [
        ('tss', 'promoter'),
        ('tts', 'TTS'),
        ('exon', 'exon'),
        ('intron', 'intron'),
        ('intergenic', 'Intergenic'),
    ]
    outdir = output.unique_output_dir("chipseq_peak_annotation_homer", reuse_empty=True)

    macs2_peaks = {}
    peaks = {}
    annot = {}
    annot_counts = {}

    for rt in run_types:
        rt_m = rt.replace('homer_', '')
        macs2_obj = loader.load_macs2_by_patient('all', run_type=rt_m)
        macs2_dat = macs2_obj.data
        for k in macs2_dat:
            macs2_dat[k] = macs2_dat[k].set_index('peak_name')
        macs2_peaks[rt_m] = macs2_dat
        peaks[rt] = loader.load_macs2_by_patient('all', run_type=rt)
        annot[rt] = {}
        this_annot_counts = {}
        for k, the_dat in peaks[rt].data.items():
            this_annot = pd.Series(index=the_dat.index)

            for u, v in annot_types:
                this_annot[the_dat.Annotation.str.contains(v)] = u
            annot[rt][k] = this_annot
            this_annot_counts[k] = this_annot.groupby(this_annot).size().to_dict()
        annot_counts[rt] = pd.DataFrame.from_dict(this_annot_counts)

    for grp in chip_groups:
        n = sum([(peaks.values()[0].meta.chip_target == ct).sum() for ct in grp])
        fig, axs = plt.subplots(
            ncols=2,
            figsize=(10., 7.2 / 8. * n),
            sharey=True
        )
        for i, rt in enumerate(run_types):
            ax = axs[i]
            # ax.set_title(rt)
            snames_all = []
            x_all = []
            j = 0  # counter for the y axis position
            for ct in grp:
                ix = peaks[rt].meta.chip_target == ct
                snames = peaks[rt].meta.index[ix].sort_values()

                this_dat = [(peaks[rt].data[sn].End - peaks[rt].data[sn].Start).values for sn in snames]
                this_log_dat = [np.log10(t) for t in this_dat]
                this_x = range(j, j + len(this_dat))

                j = len(this_dat)
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
        fig.savefig(os.path.join(outdir, "peak_size_distribution_%s.png" % '_'.join(grp)), dpi=200)

    # PLOT: dual bar chart - number of peaks and assignment %
    for rt in run_types:
        # reorder by chip target
        idx = peaks[rt].meta.sort_values(['chip_target', 'input_sample']).index
        this_annot = annot_counts[rt]

        ss = this_annot.sum().loc[idx]
        pct = this_annot.divide(ss, axis=1)[idx.tolist()] * 100.
        ss_colours = pd.Series(index=ss.index)
        for ct in peaks[rt].meta.chip_target.unique():
            ss_colours[peaks[rt].meta.chip_target == ct] = chip_target_colours[ct]

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

        fig.subplots_adjust(left=0.28, right=0.8, wspace=0.2, top=0.95, bottom=0.1)

        fig.savefig(os.path.join(outdir, "peaks_feature_assignment_%s.png" % rt), dpi=200)

    # add gene symbol and save annotated data
    peaks_data = {}
    gene_lookup = general.transcript_to_gene_lookup(tax_id=9606)

    for rt in run_types:
        rt_m = rt.replace('homer_', '')
        peaks_data[rt] = {}
        out_subdir = os.path.join(outdir, rt)
        os.makedirs(out_subdir)

        for k in peaks[rt].data:
            the_dat = peaks[rt].data[k].copy()
            the_dat.drop('Strand', axis=1, inplace=True)
            the_dat.dropna(how='all', axis=1, inplace=True)
            orig_dat = macs2_peaks[rt_m][k].loc[the_dat.index]
            the_dat.insert(5, 'basic_annotation', annot[rt][k])
            the_dat.insert(4, 'fold_change', orig_dat['fc'])
            the_dat.insert(4, '-log10(fdr)', orig_dat['-log10q'])

            the_genes = gene_lookup.loc[the_dat['Nearest PromoterID'].values]
            the_dat.insert(the_dat.shape[1], 'ensembl_gene_id', the_genes['Gene stable ID'].values)
            the_dat.insert(the_dat.shape[1], 'gene_symbol', the_genes['Gene name'].values)

            the_dat.to_csv(os.path.join(out_subdir, "%s_annotated_peaks.csv" % k))
            peaks_data[rt][k] = the_dat

