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
    parser.add_argument("-n", "--num-reads", help="Approximate number of reads to return", default=None, type=float)
    parser.add_argument("--proper-pair", help="Only include proper pairs", action="store_true", default=False)
    parser.add_argument("--include-unmapped", help="Include unmapped reads", action="store_true", default=False)
    parser.add_argument("--pe", help="Paired end sequencing", action="store_true", default=False)

    args, extra = parser.parse_known_args()
    num_reads = args.num_reads
    frac_reads = args.frac_reads

    if args.proper_pair and not args.pe:
        raise ValueError("It doesn't make sense to specify SE reads but require reads mapped in a proper pair.")

    if not os.path.isfile(args.bam_fn):
        logger.error("Invalid path to BAM file %s", args.bam_fn)
        parser.print_help()
        sys.exit(1)

    if num_reads is not None and frac_reads is not None:
        raise ValueError("Must specify EITHER num_reads OR frac_reads")

    if num_reads is None and frac_reads is None:
        logger.info("No fraction or number of reads supplied. Using default: n=1000.")
        num_reads = 1000

    it = genomics.estimate_gc_content_from_bam(
        args.bam_fn,
        frac_reads=frac_reads,
        proper_pair=args.proper_pair,
        min_qual=args.min_qual,
        est_n_reads=num_reads,
        include_unmapped=args.include_unmapped,
        extra_args=tuple(extra)
    )

    for r in it:
        sys.stdout.write("%.3f\n" % r)

    return


if __name__ == "__main__":
    main()