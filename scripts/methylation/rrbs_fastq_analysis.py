import os
import collections
import gzip
import numpy as np
from scipy import stats
import re
from settings import DATA_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
import pysam
from matplotlib import pyplot as plt
import pandas as pd
import multiprocessing as mp
from utils import log, genomics, output
logger = log.get_console_logger(__name__)


def get_first_n(fn, n=3):
    res = collections.Counter()
    fq = pysam.FastqFile(fn, 'rb')
    for rd in fq:
        res[rd.sequence[:n]] += 1
    return res



if __name__ == '__main__':
    outdir = output.unique_output_dir("rrbs_fastq_analysis")
    # basedir = os.path.join(DATA_DIR_NON_GIT, 'rrbseq', 'GC-CV-7163')
    # fq_pattern = "GC-CV-7163-{sid}_S{sid}_L{0:03d}_{rid}.fastq.gz"
    # sample_names = [
    #     'eNSC3',
    #     'eNSC5',
    #     'eNSC6',
    #     'mDura3Human',
    #     'mDura5Human',
    #     'mDura6Human',
    # ]

    basedir = os.path.join(DATA_DIR, 'rrbseq', 'GC-CV-8176-b')
    fq_pattern = "{sid}-GC-CV-8176/{sid}-GC-CV-8176_S{sid}_L{0:03d}_R{rid}_001.fastq.gz"
    sample_names = [
        'm3_choi_1',
        'm3_choi_2',
        'm6_gibco_1',
        'm6_gibco_2',
        'm5_gibco_1',
        'm5_gibco_2',
        'm3_gibco_1',
        'm3_gibco_2',
        'm6_choi_1',
        'm6_choi_2',
        'm5_choi_1',
        'm5_choi_2',
    ]

    lanes = range(1, 5)
    samples = range(1, len(sample_names) + 1)

    first3 = {}
    jobs = {}

    for s in samples:
        first3[s] = {}

        for rid in [1, 2]:
            res = collections.Counter()

            for l in lanes:
                fq_fn = os.path.join(basedir, fq_pattern.format(l, sid=s, rid=rid))
                this_res = get_first_n(fq_fn)
                for k, v in this_res.items():
                    res[k] += v

            first3[s][rid] = res

    # first 3 BP distributions
    n_triplet = 15

    # first 3 bp distribution

    for rid in [1, 2]:
        ncol = min(3, len(sample_names))
        nrow = int(np.ceil(len(sample_names) / float(ncol)))
        fig, axs = plt.subplots(nrows=nrow, ncols=ncol, sharey=True, figsize=(3 * ncol, 2.5 * nrow))
        for s in samples:
            i = s - 1
            ax = axs.flat[i]
            ttl = sample_names[i]
            pairs = sorted(first3[s][rid].items(), key=lambda x: -x[1])[:n_triplet]
            tot = sum([t[1] for t in pairs])
            pct = [t[1] / float(tot) * 100. for t in pairs]
            ax.bar(range(n_triplet), pct)
            ax.set_xticks(range(n_triplet))
            ax.set_xticklabels([t[0] for t in pairs], rotation=90)
            ax.set_title(ttl)
        axs[0, 0].set_ylabel('% reads')
        axs[1, 0].set_ylabel('% reads')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "read_%d_first3_dist.png" % rid), dpi=200)
