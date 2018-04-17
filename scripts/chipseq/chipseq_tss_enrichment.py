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
from chipseq import feature_enrichment


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file", help="BAM file to analyse")
    parser.add_argument("-b", "--bed", help="BED file containing regions of interest", required=True)
    parser.add_argument("-p", "--threads", help="Number of CPU threads [1].", type=int, default=1)
    parser.add_argument("-o", "--outdir", help="Output directory [.]", default='.')
    parser.add_argument("-Q", "--minmapqual", help="Minimum mapping quality for samtools depth [0]", default=0, type=int)

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

    CHROMS = ['%d' % i for i in range(1, 23)]

    if args.threads > 1:
        pool = mp.Pool(processes=args.threads)
        jobs = {}
    else:
        pool = None

    # Step 1: generate coverage if required
    cov_fn = os.path.join(outdir, "%s.cov.bed.gz" % filestem)

    if not os.path.isfile(cov_fn):
        cmd = "samtools depth -b {bed_file} -aa {bam_file} -Q {minmapqual}| gzip > {out_cov_file}".format(
            bed_file=os.path.abspath(args.bed),
            out_cov_file=cov_fn,
            bam_file=os.path.abspath(args.bam_file),
            minmapqual=args.minmapqual,
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
    with feature_enrichment.opener(args.bed, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for row in c:
            # start, end, name, strand
            regions[row[0]].append([int(row[1]), int(row[2]), row[3], row[5]])

    logger.info("Threads: %d. Pool: %s", args.threads, str(pool))
    traces = []

    for chrom, depth in feature_enrichment.coverage_reader(cov_fn=cov_fn, include=CHROMS):
        this_regions = regions[chrom]
        pkl_fn = os.path.join(outdir, '%s.trace.%s.pkl') % (filestem, chrom)

        if pool is None:
            logger.info("Compute traces for chromosome %s", chrom)
            this_traces, features = feature_enrichment.depth_to_trace(depth, this_regions)
            # dump this set of traces to disk
            with gzip.open(pkl_fn, 'wb') as f:
                pickle.dump({'traces': this_traces, 'bed': this_regions, 'features': features}, f)
            logger.info("Dumped traces to %s", pkl_fn)
            traces.append(this_traces)
        else:
            jobs[chrom] = pool.apply_async(feature_enrichment.depth_to_trace, args=(depth, this_regions))


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
