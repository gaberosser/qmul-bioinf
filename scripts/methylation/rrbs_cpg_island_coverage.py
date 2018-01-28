import sys
import os
import re
import pandas as pd
import multiprocessing as mp
import argparse
import json
import pysam
import datetime
from StringIO import StringIO

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from settings import DATA_DIR_NON_GIT, GIT_LFS_DATA_DIR
from utils import genomics, output, log

log_dir = os.path.join(os.environ['HOME'], 'log')
now_str = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
flogger = log.get_file_logger("rrbs_cpg_island_coverage", os.path.join(log_dir, "rrbs_cpg_island_coverage"))
# clogger = log.get_console_logger(__name__)


def get_one_coverage(bam_fn, region, region_pad=500, n_perm=100):
    chrom, start, end = region
    start = max(start - region_pad, 0)
    end += region_pad  # no problem if we run off the end

    perm_depths = [None] * n_perm

    # get the coverage over this region
    region_depth = pysam.depth(bam_fn, "-a", "-r", "%s:%d-%d" % (chrom, start, end))

    for j in range(n_perm):
        chrom, start, end = genomics.random_genomic_interval(chrom_lengths, n_bp=end - start + 1)
        perm_depths[j] = pysam.depth(bam_fn, "-a", "-r", "%s:%d-%d" % (chrom, start, end))

    return region_depth, perm_depths


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')

    required.add_argument("-i", "--input", help="Path to BAM file", required=True)
    required.add_argument("-p", "--threads", help="Number of threads", default=os.environ.get('NSLOTS', '1'))

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    ncpu = int(args.threads)
    f = os.path.abspath(args.input)

    flogger.info("Running on input file %s with %d threads.", f, ncpu)

    outdir = output.unique_output_dir("rrbs_coverage_sampling", reuse_empty=True)

    flogger.info("Output to directory %s", outdir)

    region_pad = 500  # number of BP to pad each region by
    n_perm = 10  # number of times to draw each region size at random

    bed_fn = os.path.join(GIT_LFS_DATA_DIR, 'mouse_cpg_island', 'grcm38_cpgisland.bed')
    tsv_fn = os.path.join(GIT_LFS_DATA_DIR, 'mouse_cpg_island', 'grcm38_cpgisland.tsv')
    indir = os.path.join(DATA_DIR_NON_GIT, 'rrbseq', 'GC-CV-7163', 'mouse', 'bismark')
    subdir = "GC-CV-7163-{i}_S{i}"

    chrom_lengths = genomics.reference_genome_chrom_lengths(tax_id=10090)
    # discard unplaced scaffolds, MT, X, Y
    chrom_lengths = chrom_lengths.loc[chrom_lengths.index.str.contains(r'^[0-9]')]

    # BAMs must be sorted and indexed
    # flist = [
    #     os.path.join(indir, subdir.format(i=i), 'GC-CV-7163-{i}_S{i}_pe.sorted.bam'.format(i=i))
    #     for i in range(1, 7)
    # ]

    # for f in flist:
    fstem = os.path.split(f)[-1]
    fstem = os.path.splitext(fstem)[0]
    fstem = fstem.replace('.sorted', '')
    fstem = re.sub(r'_pe$', '', fstem)

    fn = os.path.join(indir, f)
    flogger.info("Starting analysis of %s.", fn)
    s = pysam.AlignmentFile(fn, 'rb')

    # load tsv
    cpg_regions = pd.read_csv(tsv_fn, sep='\t', header=0)

    if ncpu > 1:
        pool = mp.Pool(ncpu)
        jobs = {}

    cov_cpg_islands = []
    cov_perms = []

    for i, row in cpg_regions.iterrows():
        region = (row.chrom, row.chromStart, row.chromEnd)
        kwds = {'region_pad': region_pad, 'n_perm': n_perm}
        if ncpu > 1:
            jobs[i] = pool.apply_async(
                get_one_coverage,
                args=(fn, region),
                kwds=kwds
            )
        else:
            try:
                res = get_one_coverage(fn, region, **kwds)
                cov_cpg_islands.append(res[0])
                cov_perms.append(res[1])
            except Exception:
                flogger.exception("Failed to extract region %s:%d-%d.", row.chrom, row.chromStart, row.chromEnd)

    if ncpu > 1:
        pool.close()
        for i, row in cpg_regions.iterrows():
            if (i % 1000) == 0:
                flogger.info("Region %d / %d", i, cpg_regions.shape[0])
            try:
                res = jobs[i].get(1e6)
                cov_cpg_islands.append(res[0])
                cov_perms.append(res[1])
            except Exception:
                flogger.exception("Failed to extract region %s:%d-%d.", row.chrom, row.chromStart, row.chromEnd)

    # save results
    to_write = {
        'cpg_islands': cov_cpg_islands,
        'permutations': cov_perms,
    }

    outfn = os.path.join(outdir, "%s.cpg_coverage.json" % fstem)
    flogger.info("Writing JSON results to %s", outfn)
    with open(outfn, 'wb') as fout:
        json.dump(to_write, fout)
    flogger.info("Completed file %s.", fn)
