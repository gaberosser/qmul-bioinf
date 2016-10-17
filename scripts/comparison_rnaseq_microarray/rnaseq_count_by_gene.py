import os
import sys
import optparse
import HTSeq
import dill
import csv
import collections

REF_FN = '/home/gabriel/data/Human_Genome/gtf/Homo_sapiens.GRCh38.85.gtf.gz'
XZ_DIR = '/home/gabriel/data/Xinyu_Zhang_NEB_mRNASeq_GC-XZ-2499_270115-19825815/'
XZ_SAMPLES = (
    'XZ-1-21077637/alignments/XZ-1.alignments.sorted.bam',
    'XZ-1-21077637/alignments/XZ-1.alignments.bam',
    'XZ-1-21184456',
)
CHROMS = set([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT'])


def invert_strand(iv, copy=True):
    """
    Invert the strand labels on sequence iv
    :param iv:
    :param copy: If True, make a copy first. Otherwise modify in-place.
    :return:
    """
    if copy:
        iv2 = iv.copy()
        if iv2.strand == "+":
          iv2.strand = "-"
        elif iv2.strand == "-":
          iv2.strand = "+"
        else:
          raise ValueError, "Illegal strand"
        return iv2
    else:
        if iv.strand == "+":
          iv.strand = "-"
        elif iv.strand == "-":
          iv.strand = "+"
        else:
          raise ValueError, "Illegal strand"


def count_gene_reads(bam_fn, exon_arr, strand_reverse=False):
    print "Counting genes with %s strandedness." % 'reversed' if strand_reverse else 'normal'
    bam_file = HTSeq.BAM_Reader(bam_fn)
    paired_aln = HTSeq.pair_SAM_alignments(bam_file, bundle=True)
    gene_counts = collections.Counter()

    ct = 0
    print "Iterating over BAM file for counts..."
    for aln_pair in paired_aln:
        if len(aln_pair) != 1:
            continue  # Skip multiple alignments
        aln1, aln2 = aln_pair[0]  # extract pair
        # only proceed with mapped pairs
        if (not aln1) or (not aln2) or not (aln1.aligned and aln2.aligned):
            gene_counts["_unmapped"] += 1
            continue
        if ct % 50000 == 0:
            print ct
        ct += 1

        if strand_reverse:
            invert_strand(aln1.iv, copy=False)
            invert_strand(aln2.iv, copy=False)

        iv1 = aln1.iv
        iv2 = aln2.iv

        gene_ids = set()
        for iv, val in exon_arr[iv1].steps():
            gene_ids |= val
        for iv, val in exon_arr[iv2].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            gene_counts[gene_id] += 1
        elif len(gene_ids) == 0:
            gene_counts["_no_feature"] += 1
        else:
            gene_counts["_ambiguous"] += 1

    return gene_counts


if __name__ == '__main__':
    """
    Input arguments:
    BAM file path
    Dill/pickle file path to reference_data exon array
    [-s, --stranded = yes]

    Outputs two or three files
    - BAM file sorted by name
    - CSV file with gene counts
    """

    optParser = optparse.OptionParser(
        usage="%prog [options] BAM_file DILL_file",
        description="""
        This script takes an aligned BAM file and a pre-prepared reference_data file and performs gene counting.
        The dill reference_data file is created using the script create_reference_from_gtf.py
        """,
        epilog="""
        Written by Gabriel Rosser using the tutorial created by Simon Anders (EMBL).
        2016.  Running HTSeq version %s.
        """ % HTSeq.__version__)

    optParser.add_option("-s", "--stranded", type="choice", dest="stranded",
                         choices=("yes", "reverse"), default="yes",
                         help="the direction of the strandedness. Specify 'yes' (normal direction) " +
                              "or 'reverse' (default: yes). " +
                              "'reverse' means 'yes' with reversed strand interpretation")

    if len(sys.argv) == 1:
        optParser.print_help()
        sys.exit(1)

    (opts, args) = optParser.parse_args()

    if len(args) != 2:
        sys.stderr.write(sys.argv[0] + ": Error: Please provide two input arguments.\n")
        sys.stderr.write("  Call with '-h' to get usage information.\n")
        sys.exit(1)

    # load annotated exon positions
    with open(args[1], 'rb') as f:
        exons = dill.load(f)

    # perform counting
    strand_reverse = opts.stranded == 'reverse'
    gene_count = count_gene_reads(args[0], exons, strand_reverse=strand_reverse)

    # save output CSV file
    out_path, fn = os.path.split(args[0])
    new_fn = fn.replace('.bam', '.counts.csv')
    if new_fn == fn:
        new_fn = fn + '.counts.csv'
    new_ff = os.path.join(out_path, new_fn)

    with open(new_ff, 'wb') as f:
        c = csv.writer(f)
        c.writerow(['Gene', 'Count'])
        c.writerows(
            sorted(gene_count.items(), key=lambda x: -x[1])
        )
