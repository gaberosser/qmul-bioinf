import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

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
            this_meta.insert(0, 'batch', os.path.basename(the_dir))
            these_metas.append(this_meta)
        meta[k] = pd.concat(these_metas, axis=0, sort=True)

    meta_combined = pd.concat(meta.values(), axis=0, sort=True)

    # raw reads
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.swarmplot(
        x='batch',
        y='read_count',
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


