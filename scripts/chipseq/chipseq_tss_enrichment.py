import pandas as pd
import numpy as np
import gzip
import os
import csv
import multiprocessing as mp
import collections
from utils.output import unique_output_dir
import pickle
import time


def coverage_reader(cov_fn, exclude=None, include=None):
    if os.path.splitext(cov_fn)[1].lower() == '.gz':
        opener = gzip.open
    else:
        opener = open

    if isinstance(exclude, str):
        exclude = {exclude}
    elif exclude is not None:
        exclude = set(list(exclude))

    if isinstance(include, str):
        include = {include}
    elif include is not None:
        include = set(list(include))

    with opener(cov_fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        current_ch = None
        current_res = []
        current_count = 0
        for line in c:
            if exclude is not None and line[0] in exclude:
                continue
            if include is not None and line[0] not in include:
                continue
            if line[0] != current_ch:
                if current_count > 0:
                    print "Complete"
                    yield (current_ch, np.array(current_res).astype(np.uint32))
                # new chromosome
                current_ch = line[0]
                current_count = 0
                print "Start reading depth for chromosome %s" % current_ch
                current_res = []
            current_res.append(line[1:])
            current_count += 1
            if (current_count % 1000000) == 0:
                print "Line %d" % current_count
        # yield final result
        print "Reached end of file"
        yield (current_ch, np.array(current_res).astype(np.uint32))


def depth_to_trace(depth, bed_regions):
    """
    :param depth: Numpy array with N + 1 columns, where N is the number of samples.
    The first col corresponds to genomic coordinate. The remaining cols give the coverage of each sample.
    :param bed_regions: List of regions, used to call samtools depth and generate the depth data. Each region is
    specified by 4 parameters: start coord, end coord, name, strand
    """
    bp = depth[:, 0]
    n_trace = len(bed_regions)
    n_sample = depth.shape[1] - 1
    n_depth = depth.shape[0]

    # use these indices to carry out a fast check whether to skip a region
    lower = depth[:, 0].min()
    upper = depth[:, 0].max()

    # infer the trace size from the first traces
    trace_length = max([
        reg[1] - reg[0] for reg in bed_regions
    ])

    this_traces = np.zeros((trace_length, n_sample, n_trace), dtype=np.uint32)

    # some slots may not get filled, so we need to keep a running tally
    n = 0
    skipped = 0
    features = []

    # get the coverage for each sample in each tss region
    for i, (start, end, name, strand) in enumerate(bed_regions):
        if i % 500 == 0:
            print i
        # account for 0-indexing of BED file
        start += 1
        # only need to correct the start coord - the GTF file deliberately includes one extra base at the end

        # if any part of the region is not within the coordinates of the depth array, skip
        if (end > upper) or (start < lower):
            # region is definitely not within the depth array
            skipped += 1
            continue

        # NB: searchsorted will return the array length if the element is not found
        ix0 = np.searchsorted(bp, start)

        ix1 = ix0 + trace_length
        if ix1 > n_depth:
            # this means that we are looking beyond the end of the depth array
            print "Bed region (%d, %d) extends beyond the depth array" % (start, end)
            skipped += 1
            continue
        else:
            the_trace = depth[ix0:ix1, 1:]
            if strand == '-':
                the_trace = the_trace[::-1]
            this_traces[:, :, n] = the_trace
            n += 1
            features.append(name)

    print "Computed %d traces. Skipped %d." % (n, skipped)
    return this_traces[..., :n], features


if __name__ == "__main__":
    outdir = unique_output_dir("chipseq_enrichment_traces", reuse_empty=True)

    CHROMS = ['%d' % i for i in range(1, 24)]
    # CHROMS = ['22']
    try:
        NCPU = mp.cpu_count()
    except Exception as exc:
        print repr(exc)
        NCPU = 1

    # cov_fn = '/home/gabriel/Documents/qmul_data/samtools_depth_tss_pm_2000.cov.gz'
    # cov_fn = '/home/gabriel/Documents/samtools_depth_tss_pm_2000.cov.gz'
    # gene_bed_fn = '/home/gabriel/Documents/qmul_data/tss_pm_2000.bed'
    # gene_bed_fn = '/home/gabriel/bioinf/tss_pm_2000.bed'
    # samples = [
    #     ("GBM026", "H3K27me3"),
    #     ("GBM026", "H3K36me3"),
    #     ("Icb1299", "H3K4me3"),
    #     ("Icb1299", "H3K27me3"),
    #     ("Icb1299", "BMI1"),
    #     ("Icb1299", "CHD7"),
    #     ("Dura026_NSC", "H3K4me3"),
    #     ("Dura026_NSC", "H3K27me3"),
    #     ("Dura026_NSC", "H3K36me3"),
    #     ("GBM026", "H3K4me3"),
    # ]




    # load gene bed file, indexed by chromosome
    bed = []
    regions = collections.defaultdict(list)
    with open(gene_bed_fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for row in c:
            # start, end, name, strand
            regions[row[0]].append([int(row[1]), int(row[2]), row[3], row[5]])

    trace = {}

    pool = None
    if NCPU > 1:
        pool = mp.Pool(processes=NCPU)
        jobs = {}

    for chrom, depth in coverage_reader(cov_fn=cov_fn, include=CHROMS):
        this_regions = np.array(regions[chrom])

        if pool is None:
            this_traces, features = depth_to_trace(depth, this_regions)
            print "Got traces for chromosome %s" % chrom
            # dump this set of traces to disk
            with open(os.path.join(outdir, 'trace_%s.pkl') % chrom, 'wb') as f:
                pickle.dump({'traces': this_traces, 'bed': this_regions, 'features': features}, f)
            trace[chrom] = this_traces.mean(axis=2)
        else:
            jobs[chrom] = pool.apply_async(depth_to_trace, args=(depth, this_regions))


    if pool is not None:
        pool.close()
        tic = time.time()
        # keep looping over jobs until finished or a ludicrously long time has elapsed
        while (len(jobs) > 1) and (time.time() - tic) < 50000:
            for chrom in jobs.keys():
                j = jobs[chrom]
                if j.ready():
                    this_traces, features = j.get()
                    print "Got traces for chromosome %s" % chrom
                    jobs.pop(chrom)
                    with open(os.path.join(outdir, 'trace_%s.pkl') % chrom, 'wb') as f:
                        pickle.dump({'traces': this_traces, 'bed': np.array(regions[chrom]), 'features': features}, f)
                    trace[chrom] = this_traces.mean(axis=2)


        # for chrom, j in jobs.items():
        #     this_traces, features = j.get
        #
        #
        #     # split depth with overhangs
        #     n = depth.shape[0]
        #     fracs = np.linspace(0, 1, NCPU + 1)
        #     split_ind = [int(f * n) for f in fracs]
        #
        #     splits = []
        #     split_coords= []
        #     jobs = []
        #     pool = mp.Pool(processes=NCPU)
        #
        #     for i in range(len(split_ind) - 1):
        #         # generate splits with overhang of length ONE LESS THAN the length of the fragment
        #         ix0 = split_ind[i]
        #         coord0 = bp[ix0]
        #         ix1 = split_ind[i + 1]
        #         # if we aren't on the final interval, include an overhang if necessary
        #         if i != (len(split_ind) - 2):
        #             coord1 = bp[ix1] + 2 * DISTANCE - 1
        #             # if the overhang is included in this run, we'll include it in full.
        #             # Otherwise we will just include enough of an overhang to guarantee no regions are missed.
        #             ix1 = np.where(bp <= coord1)[0][-1] + 1
        #
        #         jobs.append(
        #             pool.apply_async(depth_to_trace, args=(depth[ix0:ix1], this_bed), kwds={'distance': DISTANCE})
        #         )
        #         splits.append((ix0, ix1))
        #         split_coords.append((bp[ix0], bp[ix1 - 1]))  # this is INCLUSIVE
        #
        #     pool.close()
        #     pool.join()
        #
        #     n_trace = len(this_bed)
        #     n_sample = depth.shape[1] - 1
        #
        #     this_traces = np.zeros((2 * DISTANCE + 1, n_sample, n_trace), dtype=np.uint32)
        #     features = []
        #     i = 0
        #     for j in jobs:
        #         tt, ftr = j.get()
        #         features.append(ftr)
        #         n_res = tt.shape[2]
        #         this_traces[..., i:(i + n_res)] = tt
        #         i += n_res
        #     if i != n_trace:
        #         print "WARNING: the number of traces returned (%d) is less than the number expected (%d)" % (i, n_trace)
        #     this_traces = this_traces[..., :i]
        #
        # # dump this set of traces to disk
        # with open(os.path.join(outdir, 'trace_%s.pkl') % chrom, 'wb') as f:
        #     pickle.dump({'traces': this_traces, 'bed': this_bed, 'features': features}, f)

        # trace[chrom] = this_traces.mean(axis=2)

