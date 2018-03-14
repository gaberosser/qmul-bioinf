import os
import collections
import gzip
import numpy as np
import pickle
import sys
import re
import pysam
from matplotlib import pyplot as plt
import pandas as pd

sys.path.append(os.path.dirname(__file__) + '/../../')
from settings import DATA_DIR_NON_GIT, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
from utils import log, genomics, output

logger = log.get_console_logger(__name__)


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

    outdir = output.unique_output_dir("rrbs_theor_fragments", reuse_empty=True)
    outfile = os.path.split(bam_fn)[-1].replace('.sorted.bam', '.coverage.pkl')
    outfn = os.path.join(outdir, outfile)

    logger.info("Output file %s", outfn)

    chroms = [str(t) for t in range(1, 20)]

    basedir = os.path.join(DATA_DIR_NON_GIT, 'rrbseq', 'GC-CV-7163')

    fa_fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'mouse',
        'ensembl',
        'GRCm38.p5.r88',
        'fa',
        'Mus_musculus.GRCm38.dna.primary_assembly.fa'
    )
    fa_reader = pysam.FastaFile(fa_fn)

    # indir = os.path.join(basedir, 'trim_galore/mouse/bismark')
    # bam_files = [os.path.join(indir, 'GC-CV-7163-%d_S%d.sorted.bam' % (i, i)) for i in range(1, 7)]

    # get location of every CpG
    cpg_coords = {}
    n_cpg = 0
    for c in chroms:
        this_ref = fa_reader[c]
        it = re.finditer(r'CG', this_ref)
        cpg_coords[c] = [t.start() for t in it]
        n_cpg += len(cpg_coords[c])

    # get location of every restriction site
    ccgg_coords = {}
    n_ccgg = 0
    for c in chroms:
        this_ref = fa_reader[c]
        it = re.finditer(r'CCGG', this_ref)
        ccgg_coords[c] = [t.start() for t in it]
        n_ccgg += len(ccgg_coords[c])

    cpg_coverage = {}
    fragment_coverage = {}

    # for bam_fn in bam_files:
    s = pysam.AlignmentFile(bam_fn, 'rb')

    # for each theoretical fragment, compute the size and coverage
    ## FIXME: ARGH, this double counts pairs!!
    mspi_fragments = {}
    for c in chroms:
        print "CCGG fragment coverage. Chrom %s" % c
        this_coords = ccgg_coords[c]
        mspi_fragments[c] = []
        this_res = mspi_fragments[c]
        for i in range(1, len(ccgg_coords[c])):
            t0 = this_coords[i - 1]
            t1 = this_coords[i]
            this_cov = coverage_ignore_pairs(s, c, t0, t1)
            this_res.append([t1 - t0, this_cov])

    # get coverage of every CpG
    cpg_cov = {}
    for c in chroms:
        print "CpG coverage chromosome %s" % c
        cpg_cov[c] = [coverage_ignore_pairs(s, c, t0, t0 + 1) for t0 in cpg_coords[c]]
        # cov = s.count_coverage(c)
        # cova = np.array([cov[0][i] for i in cpg_coords[c]])
        # covc = np.array([cov[1][i] for i in cpg_coords[c]])
        # covg = np.array([cov[2][i] for i in cpg_coords[c]])
        # covt = np.array([cov[3][i] for i in cpg_coords[c]])
        # sa = cova.sum()
        # sc = covc.sum()
        # sg = covg.sum()
        # st = covt.sum()
        # pr = (sa + sg) / float(sa + sg + sc + st)
        # if pr > 0.01:
        #     print "Warning: chromosome %s has a high proportion of A and G bases where we expect C or T (%.2f)" % (
        #         c, pr
        #     )
        # cpg_cov[c] = covc + covt

    res_out = {
        'cpg_cov': cpg_cov,
        'fragment_coverage': mspi_fragments,
        'cpg_coords': cpg_coords,
        'fragment_coords': ccgg_coords
    }

    with open(outfn, 'wb') as fout:
        pickle.dump(res_out, fout)
