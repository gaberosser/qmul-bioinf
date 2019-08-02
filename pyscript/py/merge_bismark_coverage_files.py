#!/usr/bin/env python

import argparse
import re
import csv
import sys
import os
import collections
import gzip

sys.path.append(os.path.dirname(__file__) + '/../../')


def merge_coverage_files(files, chroms=None):
    """
    Merge the specified bismark methylation coverage files. These are typically gzipped, but we'll check the file
    extension to be sure. The format is (tab delimited):
    CHROM START END PCT_METH M_COUNT U_COUNT
    where END = START (It's a single CpG site), M_COUNT is the number of reads supporting a methylated site (and vice versa
    for U_COUNT)
    :param files: Iterable of files to be merged
    :param chroms: Iterable specifying the chromosomes to include (None means include all chroms - this is the default
    operation)s
    :return:
    """
    if chroms is not None:
        chroms = set([str(t) for t in chroms])

    comb_m = collections.Counter()
    comb_u = collections.Counter()

    for fn in files:
        _, ext = os.path.splitext(fn)
        if ext.lower() == '.gz':
            opener = gzip.open
        else:
            opener = open
        with opener(fn, 'rb') as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if chroms is not None and row[0] not in chroms:
                    continue
                comb_m[(row[0], int(row[1]))] += int(row[4])
                comb_u[(row[0], int(row[1]))] += int(row[5])

    # sort
    keys = sorted(comb_m)

    # combine for result
    res = []
    for k in keys:
        this_m = comb_m[k]
        this_u = comb_u[k]
        res.append(
            [k[0], k[1], k[1], "%.2f" % (this_m / float(this_m + this_u) * 100.), this_m, this_u]
        )

    return res


def write_coverage_file(res, fstream):
    writer = csv.writer(fstream, delimiter='\t')
    writer.writerows(res)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge and sum multiple bismark coverage files from the same experiment.")
    parser.add_argument('files', type=argparse.FileType('rb'), nargs='+', help='Path to coverage files')
    parser.add_argument('-o', '--output', type=str, help='Output file', required=False, default=None)
    parser.add_argument('-c', '--chroms', type=str, help='Comma-separated list of chromosomes', required=False, default=None)

    args = parser.parse_args()
    files = [t.name for t in args.files]
    chroms = None
    if args.chroms is not None:
        chroms = args.chroms.split(',')

    out_fn = args.output
    out_stream = sys.stdout

    try:
        res = merge_coverage_files(files, chroms=chroms)
        if out_fn is not None:
            out_stream = gzip.open(out_fn, 'wb')
        write_coverage_file(res, out_stream)

    finally:

        out_stream.close()



