import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

import settings
from utils import output

"""
Generate some pertinent tables and plots relating to the QC of the ChIPSeq, by batches
"""

if __name__ == "__main__":
    outdir = output.unique_output_dir()

    p170710_dirs = [
        os.path.join(settings.CHIPSEQ_DIR, 'wtchg_p170710_pilot'),
        os.path.join(settings.CHIPSEQ_DIR, 'wtchg_p170710b'),
        os.path.join(settings.CHIPSEQ_DIR, 'wtchg_p170710c'),
        os.path.join(settings.CHIPSEQ_DIR, 'wtchg_p170710'),
    ]

    p180649_dirs = [
        os.path.join(settings.CHIPSEQ_DIR, 'wtchg_p180649')
    ]

    all_dirs = {
        'P170710': p170710_dirs,
        'P180649': p180649_dirs
    }

    # load metadata
    # this holds most of the relevant information
    meta = {}
    for k, arr in all_dirs.items():
        these_metas = []
        for the_dir in arr:
            this_meta = pd.read_csv(
                os.path.join(the_dir, 'sources.csv'),
                header=0,
                index_col=0
            )
            the_batch = os.path.basename(the_dir)
            this_meta.insert(0, 'batch', the_batch)
            this_meta.index = ["%d_%s" % (t, the_batch) for t in this_meta.index]
            these_metas.append(this_meta)
        meta[k] = pd.concat(these_metas, axis=0, sort=True)

    meta_combined = pd.concat(meta.values(), axis=0, sort=True)

    # raw reads
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.swarmplot(
        x='batch',
        y='read_count',
        hue='chip_target',
        data=meta_combined,
        ax=ax
    )
    ax.set_ylabel('Raw reads (x 10 mi)')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "raw_read_counts.png"), dpi=200)

    # reads mapped (Q>10)
    meta_combined.insert(
        meta_combined.shape[1],
        'pct_mapped_q10',
        meta_combined.flagstat_reads_mapped_q10 / meta_combined.read_count * 100
    )
    ix = meta_combined.flagstat_reads_mapped_q10.isnull()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.swarmplot(
        x='batch',
        y='pct_mapped_q10',
        hue='chip_target',
        data=meta_combined.loc[~ix],
        ax=ax
    )
    ax.set_ylabel(r'% mapped $(Q \geq 10)$')
    ax.set_ylim([0, 100])
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "pct_mapped_q10.png"), dpi=200)

    # % GC in all the samples
    pct_gc_full_unmapped = {}
    pct_gc_full_q10 = {}

    pct_gc_unmapped = {}
    pct_gc_q10 = {}

    for k, arr in all_dirs.items():
        this_meta = meta[k]

        this_res_full_un = {}
        this_res_full_q10 = {}

        this_res_un = {}
        this_res_q10 = {}

        for the_dir in arr:
            the_batch = os.path.basename(the_dir)
            # check whether the gc dir and files are there
            indir_un = os.path.join(the_dir, 'human', 'bt2_alignment', 'gc_frac_unmapped')
            indir_q10 = os.path.join(the_dir, 'human', 'bt2_alignment', 'gc_frac_q10')
            if os.path.isdir(indir_un):
                for fn in os.listdir(indir_un):
                    ff = os.path.join(indir_un, fn)
                    bb = fn.replace('.sorted.bam.gc', '')
                    tmp = pd.read_csv(ff, header=None, index_col=0).index.astype(float)
                    this_res_full_un[bb] = np.array(tmp)
                    this_res_un[bb] = (np.mean(tmp), np.std(tmp))
            if os.path.isdir(indir_q10):
                for fn in os.listdir(indir_q10):
                    ff = os.path.join(indir_q10, fn)
                    bb = fn.replace('.sorted.bam.gc', '')
                    tmp = pd.read_csv(ff, header=None, index_col=0).index.astype(float)
                    this_res_full_q10[bb] = np.array(tmp)
                    this_res_q10[bb] = (np.mean(tmp), np.std(tmp))

        aa = pd.DataFrame(this_res_un, index=['mean', 'stdev']).transpose()
        aa.index = ["%s_%s" % (t, the_batch) for t in aa.index]
        pct_gc_unmapped[k] = aa
        pct_gc_unmapped[k].insert(0, 'batch', meta[k].loc[aa.index, 'batch'])
        pct_gc_unmapped[k].insert(0, 'chip_target', meta[k].loc[aa.index, 'chip_target'])
        pct_gc_full_unmapped[k] = pd.DataFrame(this_res_full_un)

        aa = pd.DataFrame(this_res_q10, index=['mean', 'stdev']).transpose()
        aa.index = ["%s_%s" % (t, the_batch) for t in aa.index]
        pct_gc_q10[k] = aa
        pct_gc_q10[k].insert(0, 'batch', meta[k].loc[aa.index, 'batch'])
        pct_gc_q10[k].insert(0, 'chip_target', meta[k].loc[aa.index, 'chip_target'])
        pct_gc_full_q10[k] = pd.DataFrame(this_res_full_q10)

    pct_gc_q10_combined = pd.concat(pct_gc_q10.values(), axis=0)
    pct_gc_unmapped_combined = pd.concat(pct_gc_unmapped.values(), axis=0)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.swarmplot(
        x='batch',
        y='mean',
        hue='chip_target',
        data=pct_gc_unmapped_combined,
        ax=ax
    )
    ax.set_ylabel('% GC before mapping')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "pct_gc_before.png"), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.swarmplot(
        x='batch',
        y='mean',
        hue='chip_target',
        data=pct_gc_q10_combined,
        ax=ax
    )
    ax.set_ylabel('% GC AFTER mapping')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "pct_gc_after.png"), dpi=200)