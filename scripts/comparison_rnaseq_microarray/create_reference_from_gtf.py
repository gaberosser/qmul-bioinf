import os
import sys
import optparse
import HTSeq
import dill

# enumerate all permissible chromosome names
CHROMS = set([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT'])


def load_exons_from_gtf(gtf_fn, stranded=True):
    gtf_file = HTSeq.GFF_Reader(gtf_fn, end_included=True)
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=stranded)
    for feat in gtf_file:
        if feat.type == 'exon' and feat.iv.chrom in CHROMS:
            if feat.iv.chrom == 'MT':
                feat.iv.chrom = 'chrM'
            else:
                feat.iv.chrom = 'chr' + feat.iv.chrom
            exons[feat.iv] += feat.name
    return exons


def save_to_pickle(exons, fn):
    with open(fn, 'wb') as f:
        dill.dump(exons, f)


if __name__ == '__main__':

    optParser = optparse.OptionParser(
        usage="%prog [options] g[ft]f_file",
        description="""
        This script takes a gtf or gff reference_data file and calculates the location of all the annotated exons.
        Assume chromosome names are simple strings, i.e. '1', 'X', etc.
        The chromosome names are modified to ensure they commence with 'chr'.
        The resulting GenomicArrayOfSets object is saved using dill for quick loading on the next run.
        """,
        epilog="""
        Written by Gabriel Rosser using the tutorial created by Simon Anders (EMBL).
        2016.  Running HTSeq version %s.
        """ % HTSeq.__version__)

    optParser.add_option("-s", "--stranded", action="store_true", dest="stranded", help="required if reference_data data are stranded")
    optParser.add_option("-o", type="string", dest="output_fn", default="exon_reference.dill", help="path to output file")

    if len(sys.argv) == 1:
        optParser.print_help()
        sys.exit(1)

    (opts, args) = optParser.parse_args()

    if len(args) != 1:
        sys.stderr.write(sys.argv[0] + ": Error: Please provide one input argument.\n")
        sys.stderr.write("  Call with '-h' to get usage information.\n")
        sys.exit(1)

    exons = load_exons_from_gtf(args[0], stranded=opts.stranded)
    save_to_pickle(exons, opts.output_fn)
