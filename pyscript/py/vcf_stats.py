import argparse
import sys
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

    samples = rd.samples
    res = {}

    for rec in rd:
        res[str(rec)] = dict([(c.sample, c.called) for c in rec.samples])

    return pd.DataFrame(res).transpose().astype(np.uint8)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get per-sample statistics relating to VCF file.")
    parser.add_argument('vcf', type=str, help='Path to VCF')
    parser.add_argument('--contig', type=str, help='Contig (optional)', required=False, default=None)
    parser.add_argument('--start', type=str, help='Start position (optional)', required=False, default=None)
    parser.add_argument('--end', type=str, help='End position (optional)', required=False, default=None)

    args = parser.parse_args()

    # TODO
    # res = get_calls_per_sample()
