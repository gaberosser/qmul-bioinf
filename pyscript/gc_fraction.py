#!/usr/bin/env python2
import sys
import os
import shutil
import tempfile
import subprocess
import argparse

sys.path.append(os.path.dirname(__file__) + '/../')
from utils import genomics


def main():
    parser = argparse.ArgumentParser(description="Estimate GC fraction (mean and stdev). Results returned to stdout.")

    parser.add_argument("bam_fn", help="Path to BAM file")
    parser.add_argument("-q","--min-qual", help="Minimum mapping quality", default=None, type=int)
    parser.add_argument("-p", "--frac-reads", help="Fraction of reads to use", default=0.001, type=float)
    parser.add_argument("-n", "--num-reads", help="Maximum number of reads to use", default=None, type=float)
    parser.add_argument("--proper-pair", help="Only include proper pairs", action="store_true", default=False)
    parser.add_argument("--include-unmapped", help="Include unmapped reads", action="store_true", default=False)

    args, extra = parser.parse_known_args()

    if args.bam_fn:
        res = genomics.estimate_gc_content_from_bam(
            args.bam_fn,
            frac_reads=args.frac_reads,
            proper_pair=args.proper_pair,
            min_qual=args.min_qual,
            n_reads=args.num_reads,
            include_unmapped=args.include_unmapped
        )
        for r in res:
            sys.stdout.write("%.3f\n" % r)
    else:
        parser.print_help()
        sys.exit(1)
    return


if __name__ == "__main__":
    main()