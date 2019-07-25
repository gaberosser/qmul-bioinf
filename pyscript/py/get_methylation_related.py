#!/usr/bin/env python

import argparse
import re
import sys
import os
import vcf
import pandas as pd
import re

sys.path.append(os.path.dirname(__file__) + '/../../')
from utils import string_manipulation


def get_methylation_related_variants(
    fn,
    search_term_arr,
    chrom=None,
    start_pos=None,
    end_pos=None
):
    """
    Extract variants linked to the terms in search_term_arr.
    We use the 'ANN' attribute in INFO to do this, running a simple text search. This may generate some false positives.
    :param fn:
    :param search_term_arr: Iterable containing search strings.
    :param chrom:
    :param start_pos:
    :param end_pos:
    :return:
    """
    if chrom is None:
        if start_pos is not None or end_pos is not None:
            raise AttributeError("If start_pos or end_pos are specified, must also specify chrom.")

    rd = vcf.Reader(filename=fn)
    if chrom is not None:
        rd = rd.fetch(chrom, start=start_pos, end=end_pos)

    the_regex = string_manipulation.construct_multiple_or_regex(search_term_arr, flags=re.IGNORECASE)
    res = []

    for rec in rd:
        the_ann = ''.join(rec.INFO['ANN'])
        srch = re.search(the_regex, the_ann)
        if srch is not None:
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


METH_GENES = [
    "A0A096LPK6", "AICDA", "ALKBH1", "ALKBH2", "ALKBH3", "APEX1", "APOBEC1", "APOBEC2", "APOBEC3A", "APOBEC3B",
    "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H", "ASZ1", "ATF7IP", "ATRX", "BAZ2A", "BEND3", "BRCA1",
    "CTCF", "CTCFL", "DDX4", "DMAP1", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "DPPA3", "EHMT1", "EHMT2", "EZH2",
    "FAM129B", "FKBP6", "FOS", "FTO", "GATA3", "GATAD2A", "GNAS", "GRHL2", "GSK3A", "GSK3B", "H1FOO", "HELLS", "HEMK1",
    "KCNQ1OT1", "KDM1B", "KMT2A", "KMT2E", "MAEL", "MBD1", "MBD2", "MBD3", "MECP2", "METTL4", "MGMT", "MIS18A", "MORC1",
    "MOV10L1", "MPHOSPH8", "MTA2", "MTRR", "MYC", "N6AMT1", "OTUD4", "PARP1", "PICK1", "PIK3CA", "PIKC3A", "PIWIL2",
    "PIWIL4", "PLD6", "PPM1D", "PRDM14", "PRMT5", "PRMT7", "RLF", "SPI1", "STPG4", "TDG", "TDRD1", "TDRD12", "TDRD5",
    "TDRD9", "TDRKH", "TET1", "TET2", "TET3", "TRIM28", "UHRF1", "UHRF2", "USP7", "USP9X", "WT1", "ZFP57", "ZMPSTE24"
]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get per-sample statistics relating to VCF file.")
    parser.add_argument('vcf', type=str, help='Path to VCF')
    parser.add_argument('--contig', type=str, help='Contig (optional)', required=False, default=None)
    parser.add_argument('--start', type=int, help='Start position (optional)', required=False, default=None)
    parser.add_argument('--end', type=int, help='End position (optional)', required=False, default=None)
    parser.add_argument('--outdir', type=str, help='Output directory', required=False, default='.')
    parser.add_argument('--search_terms', type=str, help='Comma-separated list of search terms', required=False)

    args = parser.parse_args()
    if args.search_terms is not None:
        args.search_terms = ','.split(args.search_terms)
    else:
        args.search_terms = METH_GENES

    res = get_methylation_related_variants(
        args.vcf,
        args.search_terms,
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
