import os
import sys
import HTSeq
import dill
import collections
import pandas as pd
from settings import GIT_LFS_DATA_DIR

REF_DIR = os.path.join(GIT_LFS_DATA_DIR, 'ensemble_human_genome', 'gr37')
EXONS_PRECOMPUTED_FN = os.path.join(REF_DIR, 'exons.pkl')
REF_FN = os.path.join(REF_DIR, 'Homo_sapiens.GRCh37.85.gtf.gz')

fn1 = os.path.join(
    '/home/gabriel/data',
    'Xinyu_Zhang_NEB_mRNASeq_GC-XZ-2499_270115-19825815/XZ-1-21077637/alignments/XZ-1.alignments.sorted.bam'
)

CHROMS = set([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT'])
STRAND_REVERSE = False


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


# create reference if it doesn't exist
if not os.path.exists(EXONS_PRECOMPUTED_FN):
    exon_arr = load_exons_from_gtf(REF_FN, stranded=True)
    with open(EXONS_PRECOMPUTED_FN, 'wb') as f:
        dill.dump(exon_arr, f)
else:
    with open(EXONS_PRECOMPUTED_FN, 'rb') as f:
        exon_arr = dill.load(f)

# a bit of playing with the reference
chr_len = collections.OrderedDict()
chr_gene_count = collections.OrderedDict()
chr_exon_len = collections.OrderedDict()

prev_start = 0
prev_chr = None
for s, gs in exon_arr.steps():
    if s.chrom != prev_chr:
        chr_gene_count[s.chrom] = len(gs)
        chr_exon_len[s.chrom] = s.length if bool(gs) else 0
        # start of next chrom - record length of previous and count
        # don't do this at the start
        if prev_chr is not None:
            chr_len[prev_chr] = prev_start
        # update previous chromosome end
        prev_chr = s.chrom
    else:
        chr_gene_count[s.chrom] += len(gs)
        chr_exon_len[s.chrom] += s.length if bool(gs) else 0
    prev_start = s.start
# add final chrom.
chr_len[prev_chr] = prev_start

chr_len = pd.Series(chr_len)
chr_exon_len = pd.Series(chr_exon_len)

bam_file = HTSeq.BAM_Reader(fn1)
paired_aln = HTSeq.pair_SAM_alignments(bam_file, bundle=True)
exon_counts = collections.Counter()

ct = 0
print "Iterating over BAM file for counts..."
for aln_pair in paired_aln:
    if len(aln_pair) != 1:
        continue  # Skip multiple alignments
    aln1, aln2 = aln_pair[0]  # extract pair
    # only proceed with mapped pairs
    if (not aln1) or (not aln2) or not (aln1.aligned and aln2.aligned):
        exon_counts["_unmapped"] += 1
        continue
    if ct % 50000 == 0:
        print ct
    ct += 1

    if STRAND_REVERSE:
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
        exon_counts[gene_id] += 1
    elif len(gene_ids) == 0:
        exon_counts["_no_feature"] += 1
    else:
        exon_counts["_ambiguous"] += 1

exon_counts = pd.Series(exon_counts)
