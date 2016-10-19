import os
import sys
import optparse
import re
import dill
import csv
import pandas as pd
import collections

if __name__ == '__main__':
    """
    Input arguments:
    SAM file path - this is the output from htseq-count
    Output file path - the results are saved in dill format to this
    """

    optParser = optparse.OptionParser(
        usage="%prog [options] SAM_file output_file",
        description="""
        This script takes the SAM output file generated by htseq-count and assembles a Pandas dataframe containing
        the counts.
        """,
        epilog="""
        Written by Gabriel Rosser.
        """)

    if len(sys.argv) == 1:
        optParser.print_help()
        sys.exit(1)

    (opts, args) = optParser.parse_args()

    if len(args) != 2:
        sys.stderr.write(sys.argv[0] + ": Error: Please provide two input arguments.\n")
        sys.stderr.write("  Call with '-h' to get usage information.\n")
        sys.exit(1)

    ct = collections.Counter()

    ens_regex = re.compile(r'\tXF:Z:ENS')
    amb_regex = re.compile(r'\tXF:Z:__ambiguous')
    nof_regex = re.compile(r'\tXF:Z:__no_feature')
    anu_regex = re.compile(r'\tXF:Z:__alignment_not_unique')

    with open(args[0], 'rb') as f:
        for i, l in enumerate(f):
            if (i > 0) and ((i % 100000) == 0):
                sys.stderr.write("%d\n" % i)
            if re.match(ens_regex, l):
                g = l.strip('\tXF:Z:')
                ct[g] += 1
            elif re.match(amb_regex, l):
                ct['_ambiguous'] += 1
            elif re.match(nof_regex, l):
                ct['_no_feature'] += 1
            elif re.match(anu_regex, l):
                ct['_alignment_not_unique'] += 1

    df = pd.Series(ct)

    with open(args[1], 'wb') as f:
        dill.dump(df, f)
