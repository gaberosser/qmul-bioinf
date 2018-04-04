#!/usr/bin/env python
import pandas as pd
import numpy as np
import logging
import gzip
import os
import csv
import multiprocessing as mp
import subprocess
import collections
import argparse
import sys
import pickle
import time
import re

LOG_FMT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')


def opener(fn, *args, **kwargs):
    if os.path.splitext(fn)[1].lower() == '.gz':
        return gzip.open(fn, *args, **kwargs)
    else:
        return open(fn, *args, **kwargs)


def coverage_reader(cov_fn, exclude=None, include=None):
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
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file", help="BAM file to analyse")
    parser.add_argument("-b", "--bed", help="BED file containing regions of interest", required=True)
    parser.add_argument("-p", "--threads", help="Number of CPU threads [1].", type=int, default=1)
    parser.add_argument("-o", "--outdir", help="Output directory [.]", default='.')

    args = parser.parse_args()
    outdir = os.path.abspath(args.outdir)

    if not os.path.exists(args.bam_file):
        raise ValueError("Invalid path to BAM file.")

    if not os.path.exists(args.bed):
        raise ValueError("Invalid path to BED file.")

    filestem = os.path.split(args.bam_file)[1]
    filestem = re.sub(r'(\.sorted)?\.bam', '', filestem)

    logger = logging.getLogger("chipseq_tss_enrichment")

    # reset handlers
    logger.handlers = []
    sh = logging.StreamHandler()
    fmt = logging.Formatter(LOG_FMT)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    logger.setLevel(logging.INFO)

    CHROMS = ['%d' % i for i in range(1, 24)]

    if args.threads > 1:
        pool = mp.Pool(processes=args.threads)
        jobs = {}
    else:
        pool = None

    # Step 1: generate coverage if required
    cov_fn = os.path.join(outdir, "%s.cov.bed.gz" % filestem)

    if not os.path.isfile(cov_fn):
        cmd = "samtools depth -b {bed_file} -aa {bam_file} | gzip > {out_cov_file}".format(
            bed_file=os.path.abspath(args.bed),
            out_cov_file=cov_fn,
            bam_file=os.path.abspath(args.bam_file)
        )
        logger.info("Calling samtools depth.")
        logger.info("%s", cmd)

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, cwd=outdir)
        stdout, stderr = p.communicate()
        logger.info("Stdout: %s", stdout)
        logger.info("Stderr: %s", stderr)
        if p.returncode != 0:
            logger.error("samtools depth gave a non-zero return code (%s)", str(p.returncode))
        logger.info("Done.")
    else:
        logger.info("Using existing coverage file %s", cov_fn)

    trace = {}

    # load gene bed file, indexed by chromosome
    bed = []
    regions = collections.defaultdict(list)
    with opener(args.bed, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for row in c:
            # start, end, name, strand
            regions[row[0]].append([int(row[1]), int(row[2]), row[3], row[5]])

    logger.info("Threads: %d. Pool: %s", args.threads, str(pool))
    traces = []

    for chrom, depth in coverage_reader(cov_fn=cov_fn, include=CHROMS):
        this_regions = regions[chrom]
        pkl_fn = os.path.join(outdir, '%s.trace.%s.pkl') % (filestem, chrom)

        if pool is None:
            logger.info("Compute traces for chromosome %s", chrom)
            this_traces, features = depth_to_trace(depth, this_regions)
            # dump this set of traces to disk
            with gzip.open(pkl_fn, 'wb') as f:
                pickle.dump({'traces': this_traces, 'bed': this_regions, 'features': features}, f)
            logger.info("Dumped traces to %s", pkl_fn)
            traces.append(this_traces)
        else:
            jobs[chrom] = pool.apply_async(depth_to_trace, args=(depth, this_regions))


    if pool is not None:
        pool.close()
        tic = time.time()
        # keep looping over jobs until finished
        while len(jobs) > 1:
            for chrom in jobs.keys():
                j = jobs[chrom]
                if j.ready():
                    try:
                        this_traces, features = j.get(600)
                        logger.info("Completed traces for chromosome %s", chrom)
                        pkl_fn = os.path.join(outdir, '%s.trace.%s.pkl') % (filestem, chrom)
                        with gzip.open(pkl_fn, 'wb') as f:
                            pickle.dump(
                                {'traces': this_traces, 'bed': regions[chrom], 'features': features},
                                f
                            )
                            logger.info("Dumped traces to %s", pkl_fn)
                            traces.append(this_traces)
                    except Exception:
                        logger.exception("Failed to compute trace for chromosome %s", chrom)
                    jobs.pop(chrom)

    # finally aggregate over all traces
    mean_trace = np.dstack(traces).mean(axis=2).squeeze()
    trace_fn = os.path.join(outdir, "%s.trace" % filestem)

    with open(trace_fn, 'wb') as f:
        f.write(", ".join(["%.6f" % t for t in mean_trace]))

    logger.info("Wrote mean trace to %s.", trace_fn)
