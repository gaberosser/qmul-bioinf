#!/usr/bin/env python2
import sys
import os
import shutil
import tempfile
import subprocess
import argparse

sys.path.append(os.path.dirname(__file__) + '/../')
from utils import genomics, log
logger = log.get_console_logger()


def main():
    parser = argparse.ArgumentParser(description="Estimate GC fraction (mean and stdev). Results returned to stdout.")

    parser.add_argument("bam_fn", help="Path to BAM file")
    parser.add_argument("-q","--min-qual", help="Minimum mapping quality", default=None, type=int)
    parser.add_argument("-p", "--frac-reads", help="Fraction of reads to use", default=None, type=float)
    parser.add_argument("-n", "--num-reads", help="Maximum number of reads to use", default=None, type=float)
    parser.add_argument("--proper-pair", help="Only include proper pairs", action="store_true", default=False)
    parser.add_argument("--include-unmapped", help="Include unmapped reads", action="store_true", default=False)

    args, extra = parser.parse_known_args()
    num_reads = args.num_reads
    frac_reads = args.frac_reads

    if not os.path.isfile(args.bam_fn):
        logger.error("Invalid path to BAM file %s", args.bam_fn)
        parser.print_help()
        sys.exit(1)

    if num_reads is None:
        if frac_reads is None:
            logger.info("No fraction or number of reads supplied. Using default: n=1000.")
            num_reads = 1000

    if frac_reads is None:
        if num_reads is not None:
            logger.info("Using estimated BAM size to get a representative sample of %d reads.", num_reads)
            est_len = genomics.estimate_number_of_bam_reads(args.bam_fn)
            frac_reads = 10. / float(est_len) * num_reads
            logger.info("Using p=%.4e", frac_reads)

    it = genomics.estimate_gc_content_from_bam(
        args.bam_fn,
        frac_reads=frac_reads,
        proper_pair=args.proper_pair,
        min_qual=args.min_qual,
        n_reads=num_reads,
        include_unmapped=args.include_unmapped
    )
    for r in it:
        sys.stdout.write("%.3f\n" % r)

    return


if __name__ == "__main__":
    main()