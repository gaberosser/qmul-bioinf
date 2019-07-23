#!/usr/bin/env python

import argparse
import re
import os
import vcf
import pandas as pd
import numpy as np


def get_gbm_specific_variants(
        fn,
        chrom=None,
        start_pos=None,
        end_pos=None
):
    """
    Extract variants that are found in the GBM sample but NOT the matching iNSC (or vice versa) in the specified region.
    We assume that the naming convention is GBM<pid> and iNSC<pid>
    TODO: make this more flexible?
    :param fn:
    :param chrom:
    :param start_pos:
    :param end_pos:
    :return:
    """
    if chrom is None:
        if start_pos is not None or end_pos is not None:
            raise AttributeError("If start_pos or end_pos are specified, must also specify chrom.")

    rd = vcf.Reader(filename=fn)
    samples = rd.samples

    sample_pairs = {}
    res = []

    for s in samples:
        regex = re.search(r'(?P<pid>[0-9]*)$', s)
        pid = regex.group('pid')
        sample_pairs.setdefault(pid, []).append(s)

    for k in sample_pairs:
        sample_pairs[k] = sorted(sample_pairs[k])
        if len(sample_pairs[k]) != 2:
            raise AttributeError("Sample with ID %s has %d matching entries. Expected 2." % (k, len(sample_pairs[k])))
        # res[k] = []

    if chrom is not None:
        rd = rd.fetch(chrom, start=start_pos, end=end_pos)

    for rec in rd:
        for k, pair in sample_pairs.items():
            res0 = rec.genotype(pair[0])
            res1 = rec.genotype(pair[1])
            if res0['GT'] != res1['GT']:
                res.append(rec)
    return res


def write_results_to_vcf(in_vcf_fn, out_vcf_fn, res, chrom=None):
    rd = vcf.Reader(filename=in_vcf_fn)
    if chrom is not None:
        # manually override the contigs list so it only contains the relevant chromosome
        rd.contigs = {chrom: rd.contigs[chrom]}
    with open(out_vcf_fn, 'wb') as f:
        writer = vcf.Writer(f, rd)
        for rec in res:
            writer.write_record(rec)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get per-sample statistics relating to VCF file.")
    parser.add_argument('vcf', type=str, help='Path to VCF')
    parser.add_argument('--contig', type=str, help='Contig (optional)', required=False, default=None)
    parser.add_argument('--start', type=int, help='Start position (optional)', required=False, default=None)
    parser.add_argument('--end', type=int, help='End position (optional)', required=False, default=None)
    parser.add_argument('--outdir', type=str, help='Output directory', required=False, default='.')

    args = parser.parse_args()

    res = get_gbm_specific_variants(
        args.vcf,
        chrom=args.contig,
        start_pos=args.start,
        end_pos=args.end
    )

    n = len(res)
    if n == 0:
        print "No variant calls found in the specified VCF file and range."
    else:
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

        out_fn = os.path.join(outdir, "%s.vcf" % out_file)
        write_results_to_vcf(args.vcf, out_fn, res, chrom=args.contig)
        print "Wrote %d variant calls to %s" % (n, out_fn)
