from rnaseq import loader
from utils import output

import pandas as pd
import os
from matplotlib import pyplot as plt
import seaborn as sns


if __name__ == "__main__":
    """
    Aim:

    Help decide whether to opt for an additional lane of sequencing in the latest batch of SmartSeq2 RNA-Seq data.
    To support this decision, we want to consider three variables across some Poly(A) and two SS2 batches:
    a) Raw read count
    b) Mapped read count
    c) Assigned read count (to genes)
    All of these are based on the STAR pipeline and are already stored in the metadata
    """
    outdir = output.unique_output_dir()

    loaders = {
        'P180059': loader.wtchg_p180059,
        'P180347': loader.wtchg_p180347,
        'P170503': loader.wtchg_p170503,
        'P180349': loader.wtchg_p180349,
    }

    # load all metadata
    meta = dict([
        (k, pd.read_csv(v.meta_file, sep=',', header=0, index_col=0)) for k, v in loaders.items()
    ])

    # combine relevant data into a single big output table
    cols_to_keep = [
        'sample',
        'read_count',
        'uniquely_mapped_GRCh38',
        'pct_assigned_GRCh38',
        'pct_remain_after_dedupe_1',
        'pct_remain_after_dedupe_2'
    ]

    meta_all = []
    for k, v in meta.items():
        this = v[cols_to_keep]
        this.insert(1, 'batch', k)
        meta_all.append(this)
    meta_all = pd.concat(meta_all, axis=0)
    meta_all.insert(5, 'assigned_GRCh38', meta_all.read_count * meta_all.pct_assigned_GRCh38 / 100.)
    meta_all.to_excel(os.path.join(outdir, "metadata_all.xlsx"))

    # compare raw read counts
    ax = sns.swarmplot(data=meta_all, x='batch', y='read_count')
    ax.figure.savefig(os.path.join(outdir, 'raw_read_counts.png'), dpi=200)
    ax.cla()

    # mapped read counts
    ax = sns.swarmplot(data=meta_all, x='batch', y='uniquely_mapped_GRCh38')
    ax.figure.savefig(os.path.join(outdir, 'uniquely_mapped_read_counts.png'), dpi=200)
    ax.cla()

    # assigned read counts
    ax = sns.swarmplot(data=meta_all, x='batch', y='assigned_GRCh38')
    ax.figure.savefig(os.path.join(outdir, 'assigned_read_counts.png'), dpi=200)
    ax.cla()

    # compare assigned pct
    ax = sns.swarmplot(data=meta_all, x='batch', y='pct_assigned_GRCh38')
    ax.set_ylim([0, 100])
    ax.figure.savefig(os.path.join(outdir, 'pct_assigned.png'), dpi=200)
    ax.cla()
