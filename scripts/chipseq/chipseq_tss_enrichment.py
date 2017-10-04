import pandas as pd
import numpy as np
import gzip
import os
import csv


def coverage_reader(cov_fn, exclude=None):
    if isinstance(exclude, str):
        exclude = {exclude}
    elif exclude is not None:
        exclude = set(list(exclude))

    with gzip.GzipFile(cov_fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        current_ch = None
        current_res = []
        current_count = 0
        for line in c:
            if exclude is not None and line[0] in exclude:
                continue
            if line[0] != current_ch:
                if current_count > 0:
                    yield (current_ch, np.array(current_res).astype(np.uint32))
                # new chromosome
                current_ch = line[0]
                current_count = 0
                print "Start reading depth for chromosome %s" % current_ch
                current_res = []
            current_res.append(line[1:])  # ints are more space efficient and useful
            current_count += 1
            if (current_count % 1000000) == 0:
                print "Line %d" % current_count
        # yield final result
        yield (current_ch, np.array(current_res).astype(np.uint32))


if __name__ == "__main__":
    distance = 2000
    cov_fn = '/home/gabriel/Documents/qmul_data/samtools_depth_tss_pm_2000.cov.gz'
    gene_bed_fn = '/home/gabriel/Documents/qmul_data/tss_pm_2000.bed'
    samples = [
        ("GBM026", "H3K27me3"),
        ("GBM026", "H3K36me3"),
        ("Icb1299", "H3K4me3"),
        ("Icb1299", "H3K27me3"),
        ("Icb1299", "BMI1"),
        ("Icb1299", "CHD7"),
        ("Dura026_NSC", "H3K4me3"),
        ("Dura026_NSC", "H3K27me3"),
        ("Dura026_NSC", "H3K36me3"),
        ("GBM026", "H3K4me3"),
    ]

    # load gene bed file
    bed = []
    with open(gene_bed_fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for line in c:
            bed.append([line[0], int(line[1]), int(line[2]), line[3]])

    # exclusion list of non-standard chroms
    excl_chroms = [t for t in sorted(set([t[0] for t in bed])) if t[0] in ('G', 'K')]
    trace = {}

    for chrom, depth in coverage_reader(cov_fn=cov_fn, exclude=excl_chroms):

        # this_trace = np.zeros((2 * distance + 1, len(samples)))
        this_trace = np.zeros((2 * distance, len(samples)))
        this_n = 0.
        # get the coverage for each sample in each tss region
        for t in bed:
            if t[0] == chrom:
                # ix = (t[1] <= depth[:, 0]) & (t[2] >= depth[:, 0])
                ix = (t[1] <= depth[:, 0]) & (t[2] >= depth[:, 0])
                this_trace += depth[ix, 1:]
                this_n += 0.
        trace[chrom] = this_trace / float(this_n)

