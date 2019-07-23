#!/usr/bin/env python

import argparse
import os
import vcf
import pandas as pd
import numpy as np


def get_calls_per_sample(fn, chrom=None, start_pos=None, end_pos=None):
    if chrom is None:
        if start_pos is not None or end_pos is not None:
            raise AttributeError("If start_pos or end_pos are specified, must also specify chrom.")

    rd = vcf.Reader(filename=fn)
    if chrom is not None:
        rd = rd.fetch(chrom, start=start_pos, end=end_pos)

    res = {}

    for rec in rd:
        res[str(rec)] = dict([(c.sample, c['GT'] if c.called else None) for c in rec.samples])

    return pd.DataFrame(res).transpose().fillna('')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get per-sample statistics relating to VCF file.")
    parser.add_argument('vcf', type=str, help='Path to VCF')
    parser.add_argument('--contig', type=str, help='Contig (optional)', required=False, default=None)
    parser.add_argument('--start', type=int, help='Start position (optional)', required=False, default=None)
    parser.add_argument('--end', type=int, help='End position (optional)', required=False, default=None)
    parser.add_argument('--outdir', type=str, help='Output directory', required=False, default='.')

    args = parser.parse_args()

    outdir = args.outdir
    out_file = os.path.split(args.vcf)[-1]
    out_file = os.path.splitext(out_file)[0]
    if args.contig is not None:
        out_file = "%s.%s" % (out_file, args.contig)
    if args.start is not None:
        out_file = "%s.%d" % (out_file, args.start)
    if args.end is not None:
        # if no start was specified, include zero in the filename for clarity
        if args.start is None:
            out_file = "%s.0" % out_file
        out_file = "%s-%d" % (out_file, args.end)

    out_fn = os.path.join(outdir, "%s.calls" % out_file)

    res = get_calls_per_sample(
        args.vcf,
        chrom=args.contig,
        start_pos=args.start,
        end_pos=args.end
    )

    if res.shape[0] > 0:
        res.to_csv(out_fn)
        print "Saved %d variant calls to %s" % (
            res.shape[0],
            out_fn
        )
    else:
        print "No variant calls found in the specified VCF file and range."