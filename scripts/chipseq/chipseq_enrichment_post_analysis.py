import os
import pickle
from settings import OUTPUT_DIR
import numpy as np
import csv
import pandas as pd


if __name__ == "__main__":

    # get region overlap
    CHROMS = ['%d' % i for i in range(1, 24)]
    gene_bed_fn = '/home/gabriel/bioinf/tss_pm_2000.bed'
    bed = []
    with open(gene_bed_fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for line in c:
            bed.append([line[0], int(line[1]), int(line[2]), line[3]])
    bed = [t for t in bed if t[0] in CHROMS]

    overlaps = []
    for chrom in CHROMS:
        this_bed = [t for t in bed if t[0] == chrom]
        this_ol = []
        for i, (_, s0, e0, ftr0) in enumerate(this_bed):
            ol = dict([(j, 1.) for j in range(s0, e0 + 1)])
            # ol = 1.
            other_idx = range(i) + range(i + 1, len(this_bed))
            for j in other_idx:
                _, s1, e1, ftr1 = this_bed[j]
                if (s1 <= s0 <= e1) or (s0 <= s1 <= e0) or (s1 <= e0 <= e1) or (s0 <= e1 <= e0):
                    # some overlap
                    for k in range(s1, e1 + 1):
                        if k in ol:
                            ol[k] += 1

            this_ol.append(ol)
        overlaps.append(this_ol)

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
    distance = 2000
    min_reads = 10  # minimum number of assigned reads for any single sample
    # only work with non-GBM samples
    keep_idx = ['GBM' not in t[0] for t in samples]

    coords = range(-distance, distance + 1)

    all_traces = None
    all_traces_norm = None

    indir = os.path.join(OUTPUT_DIR, 'chipseq_enrichment_traces')
    chroms = [str(t) for t in range(1, 23)]
    for i in range(1, 23):
        chrom = str(i)
        fn = os.path.join(indir, 'trace_%d.pkl' % i)
        with open(fn, 'rb') as f:
            res = pickle.load(f)
        tra = res['traces']
        tra = tra[:, keep_idx, :]

        n = tra.shape[-1]
        print "Chrom %s. %d raw traces." % (chrom, n)

        # reassemble features
        features = reduce(lambda x, y: x + y, res['features'])

        # only keep traces that have some minmum number of reads
        keep_enough_reads = [i for i in range(n) if sum(tra[..., i].sum(axis=0) >= min_reads) > 0]

        print "Keeping %d traces with sufficient reads (%.2f %%)." % (
            len(keep_enough_reads),
            len(keep_enough_reads) / float(n) * 100.
        )

        if len(keep_enough_reads) == 0:
            continue

        tra = tra[..., keep_enough_reads]

        # normalise each trace??
        tra_norm = np.zeros_like(tra).astype(float)
        for i in range(tra.shape[-1]):
            tra_norm[..., i] = tra[..., i] / (tra[..., i].sum(axis=0).astype(float) + 1.)

        if all_traces is None:
            all_traces = tra
        else:
            all_traces = np.concatenate((all_traces, tra), axis=2)

        if all_traces_norm is None:
            all_traces_norm = tra_norm
        else:
            all_traces_norm = np.concatenate((all_traces_norm, tra_norm), axis=2)