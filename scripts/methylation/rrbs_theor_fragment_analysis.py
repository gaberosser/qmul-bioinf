import os
import collections
import gzip
import numpy as np
import pickle
import sys
import re
import pysam
import pybedtools
import subprocess
from matplotlib import pyplot as plt
import pandas as pd

sys.path.append(os.path.dirname(__file__) + '/../../')
from settings import DATA_DIR_NON_GIT, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
from utils import log, genomics, output

logger = log.get_console_logger(__name__)


def get_motif_locations(fa_reader, motif, references):
    for c in references:
        this_ref = fa_reader[c]
        it = re.finditer(motif, this_ref)
        for t in it:
            yield (c, t.start())


def create_cpg_bed(fa_fn, outfn, references=None):
    fa_reader = pysam.FastaFile(fa_fn)
    if references is None:
        references = fa_reader.references
    # get location of every CpG
    bed_arr = []
    for c, st in get_motif_locations(fa_reader, r'CG', references):
        bed_arr.append(
            (c, st, st + 1)
        )
    x = pybedtools.BedTool(bed_arr)
    x.saveas(outfn)


def create_ccgg_fragment_bed(fa_fn, outfn, references=None):
    fa_reader = pysam.FastaFile(fa_fn)
    if references is None:
        references = fa_reader.references
    # get location of every CCGG
    it = get_motif_locations(fa_reader, r'CCGG', references)
    c0, t0 = it.next()
    # the actual split occurs after the first C
    t0 += 1
    bed_arr = []
    for c, st in it:
        if c0 == c:
            bed_arr.append(
                (c, t0, st)
            )
        c0 = c
        t0 = st + 1
    x = pybedtools.BedTool(bed_arr)
    x.saveas(outfn)


def create_ccgg_fragment_saf(fa_fn, outfn, references=None):
    fa_reader = pysam.FastaFile(fa_fn)
    if references is None:
        references = fa_reader.references
    # get location of every CCGG
    it = get_motif_locations(fa_reader, r'CCGG', references)
    c0, t0 = it.next()
    # the actual split occurs after the first C (+1)
    # GTF is 1-indexed (+1)
    t0 += 2
    dat = []
    for c, st in it:
        # GTF is 1-indexed (+1)
        st += 1
        if c0 == c:
            dat.append(
                (c, t0, st)
            )
        c0 = c
        # next fragment starts one base later (+1)
        t0 = st + 1
    x = pd.DataFrame(dat, columns=['Chr', 'Start', 'End'])
    x.index.name = "GeneID"  # this has to be GeneID or we get an error from featureCounts?
    x.insert(3, 'Strand', pd.Series('+', index=x.index))
    x.to_csv(outfn, sep='\t')


def coverage_ignore_pairs(s, chrom, ix0, ix1):
    reads_seen = set()
    for rd in s.fetch(chrom, ix0, ix1):
        # check that the mate of this read hasn't already been counted
        if s.mate(rd) not in reads_seen:
            reads_seen.add(rd)
    return len(reads_seen)


if __name__ == "__main__":
    """
    Usage: rrbs_theor_fragment_analysis.py <BAM_FN>
    BAM_FN must be a sorted bam file
    """
    bam_fn = sys.argv[1]
    if not os.path.isfile(bam_fn):
        raise ValueError("Unable to find BAM file %s" % bam_fn)

    bam_dir = os.path.split(os.path.abspath(bam_fn))[0]
    outdir = output.unique_output_dir("rrbs_fragment_analysis", root_output_dir=bam_dir)

    # fixed output directory for BED regions
    bed_outdir = os.path.join(output.OUTPUT_DIR, "rrbs_theor_fragments")
    if not os.path.exists(bed_outdir):
        logger.info("Created output dir %s", bed_outdir)
        os.makedirs(bed_outdir)

    # same output directory for remaining results
    outfile = re.sub(r'(\.sorted)\.bam', ".coverage.pkl", os.path.split(bam_fn)[-1])
    outfn = os.path.join(outdir, outfile)
    fcounts_outfn = os.path.join(
        outdir,
        re.sub(r'(\.sorted)\.bam', '.mspi_fragments.counts', os.path.split(bam_fn)[-1])
    )
    depth_outfn = os.path.join(
        outdir,
        re.sub(r'(\.sorted)\.bam', '.cpg_coverage.gz', os.path.split(bam_fn)[-1])
    )

    logger.info("Output pickle file %s", outfn)
    logger.info("Output featureCount file %s", fcounts_outfn)

    chroms = [str(t) for t in range(1, 20)]

    # reference fasta file - may not be needed if we have already got the required BED files
    fa_fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'mouse',
        'ensembl',
        'GRCm38.p5.r88',
        'fa',
        'Mus_musculus.GRCm38.dna.primary_assembly.fa'
    )

    # do the required BED files exist? If not, create them now
    cg_bed_fn = os.path.join(bed_outdir, "cpg_regions.bed")
    if not os.path.exists(cg_bed_fn):
        logger.info("Unable to find BED file %s. Creating now.", cg_bed_fn)
        create_cpg_bed(fa_fn, cg_bed_fn, references=chroms)
        logger.info("Done")
    else:
        logger.info("Using existing BED file %s", cg_bed_fn)

    ccgg_saf_fn = os.path.join(bed_outdir, "ccgg_fragments.saf")
    if not os.path.exists(ccgg_saf_fn):
        logger.info("Unable to find SAF file %s. Creating now.", ccgg_saf_fn)
        create_ccgg_fragment_bed(fa_fn, ccgg_saf_fn, references=chroms)
        logger.info("Done")
    else:
        logger.info("Using existing SAF file %s", ccgg_saf_fn)

    # check whether BAM is sorted
    cmd = ""
    if genomics.bam_is_sorted(bam_fn):
        logger.info("BAM file %s is unsorted; we'll sort before running samtools depth.", bam_fn)
        cmd = "samtools depth -a -b {cg_bed_fn} {bam_fn}"
    else:
        cmd += "samtools sort {bam_fn} | samtools depth -a -b {cg_bed_fn}"
    cmd += " | gzip > {outfn}"
    cmd = cmd.format(
        cg_bed_fn=cg_bed_fn,
        bam_fn=bam_fn,
        outfn=depth_outfn
    )

    # CpG coverage
    logger.info("Calling samtools depth to get CpG coverage")
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()
    logger.info("Stdout: %s", stdout)
    logger.info("Stderr: %s", stderr)
    ccgg_coverage = None
    if p.returncode != 0:
        logger.error("samtools depth gave a non-zero return code (%s)", str(p.returncode))
    logger.info("Done.")

    # MspI theoretical fragment coverage
    cmd = "featureCounts -p -a {saf_fn} -T 12 -F SAF {bam_fn} -o {outfn}".format(
        saf_fn=ccgg_saf_fn,
        bam_fn=bam_fn,
        outfn=fcounts_outfn
    )
    logger.info("Calling featureCounts with the following command:")
    logger.info(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()
    logger.info("Stdout: %s", stdout)
    logger.info("Stderr: %s", stderr)
    ccgg_coverage = None
    if p.returncode != 0:
        logger.error("featureCounts gave a non-zero return code (%s)", str(p.returncode))

    # ccgg_coverage = pd.read_csv(fcounts_outfn, sep='\t', comment='#', header=0, index_col=0)
