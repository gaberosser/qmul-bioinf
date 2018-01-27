import sys
import os
import re
import pandas as pd
import multiprocessing as mp
import json
import pysam
from StringIO import StringIO

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from settings import DATA_DIR_NON_GIT, GIT_LFS_DATA_DIR
from utils import genomics, output, log

clogger = log.get_console_logger(__name__)


def get_one_coverage(bam_fn, region, region_pad=500, n_perm=100):
    chrom, start, end = region
    start = max(start - region_pad, 0)
    end += region_pad  # no problem if we run off the end

    perm_depths = [None] * n_perm

    # get the coverage over this region
    region_depth = pysam.depth(bam_fn, "-a", "-r", "%s:%d-%d" % (chrom, start, end))
    # region_depth = pd.read_csv(StringIO(this_st_depth), sep='\t', header=None)
    # region_depth.columns = ['chrom', 'coord', 'cov']

    # this_st_depth = pysam.depth(bam_fn, "-r", "%s:%d-%d" % (chrom, start, end))
    # if len(this_st_depth) > 0:
    #     # more parsing is needed
    #     region_depth = pd.read_csv(StringIO(this_st_depth), sep='\t', header=None)
    #     region_depth.columns = ['chrom', 'coord', 'cov']

    for j in range(n_perm):
        chrom, start, end = genomics.random_genomic_interval(chrom_lengths, n_bp=end - start + 1)
        perm_depths[j] = pysam.depth(bam_fn, "-a", "-r", "%s:%d-%d" % (chrom, start, end))
        # this_depth = pd.read_csv(StringIO(this_st_depth), sep='\t', header=None)
        # this_depth.columns = ['chrom', 'coord', 'cov']
        # perm_depths[j] = this_depth

        # this_st_depth = pysam.depth(bam_fn, "-r", "%s:%d-%d" % (chrom, start, end))
        # if len(this_st_depth) > 0:
        #     this_depth = pd.read_csv(StringIO(this_st_depth), sep='\t', header=None)
        #     this_depth.columns = ['chrom', 'coord', 'cov']
        #     perm_depths[j] = this_depth

    return region_depth, perm_depths


if __name__ == "__main__":

    outdir = output.unique_output_dir("rrbs_coverage_sampling", reuse_empty=True)

    region_pad = 500  # number of BP to pad each region by
    n_perm = 5  # number of times to draw each region size at random

    bed_fn = os.path.join(GIT_LFS_DATA_DIR, 'mouse_cpg_island', 'grcm38_cpgisland.bed')
    tsv_fn = os.path.join(GIT_LFS_DATA_DIR, 'mouse_cpg_island', 'grcm38_cpgisland.tsv')
    indir = os.path.join(DATA_DIR_NON_GIT, 'rrbseq', 'GC-CV-7163', 'mouse', 'bismark')
    subdir = "GC-CV-7163-{i}_S{i}"

    chrom_lengths = genomics.reference_genome_chrom_lengths(tax_id=10090)
    # discard unplaced scaffolds, MT, X, Y
    chrom_lengths = chrom_lengths.loc[chrom_lengths.index.str.contains(r'^[0-9]')]

    # BAMs must be sorted and indexed
    flist = [
        os.path.join(indir, subdir.format(i=i), 'GC-CV-7163-{i}_S{i}_pe.sorted.bam'.format(i=i))
        for i in range(1, 7)
    ]

    for f in flist:
        fstem = os.path.split(f)[-1]
        fstem = os.path.splitext(fstem)[0]
        fstem = fstem.replace('.sorted', '')
        fstem = re.sub(r'_pe$', '', fstem)

        fn = os.path.join(indir, f)
        clogger.info("Starting analysis of %s.", fn)
        s = pysam.AlignmentFile(fn, 'rb')

        # load tsv
        cpg_regions = pd.read_csv(tsv_fn, sep='\t', header=0)

        # MP
        if 'NSLOTS' in os.environ:
            ncpu = int(os.environ['NSLOTS'])
        else:
            ncpu = mp.cpu_count()
        clogger.info("Using %d parallel processes", ncpu)

        pool = mp.Pool(ncpu)
        jobs = {}

        cov_cpg_islands = []
        cov_perms = []

        for i, row in cpg_regions.iterrows():
            region = (row.chrom, row.chromStart, row.chromEnd)
            jobs[i] = pool.apply_async(
                get_one_coverage,
                args=(fn, region),
                kwds={'region_pad': region_pad, 'n_perm': n_perm}
            )

        pool.close()
        # pool.join()
        for i, row in cpg_regions.iterrows():
            ## DEBUG!
            if i > 200:
                break
            if (i % 1000) == 0:
                clogger.info("Region %d / %d", i, cpg_regions.shape[0])
            try:
                res = jobs[i].get(1e6)
                cov_cpg_islands.append(res[0])
                cov_perms.append(res[1])
            except Exception:
                clogger.exception("Failed to extract region %s:%d-%d.", row.chrom, row.chromStart, row.chromEnd)

        # save results
        to_write = {
            'cpg_islands': cov_cpg_islands,
            'permutations': cov_perms,
        }

        outfn = os.path.join(outdir, "%s.cpg_coverage.json" % fstem)
        clogger.info("Writing JSON results to %s", outfn)
        with open(outfn, 'wb') as fout:
            json.dump(to_write, fout)
        clogger.info("Completed file %s.", fn)


            # for each region, get reads
    # for i, row in cpg_regions.iterrows():
        # chrom = row.chrom
        # start = max( - region_pad, 0)
        # end = row.chromEnd + region_pad  # no problem if we run off the end
        # # get the coverage over this region
        # try:
        #     this_st_depth = pysam.depth(fn, "-r", "%s:%d-%d" % (chrom, start, end))
        #     if len(this_st_depth) > 0:
        #         # more parsing is needed
        #         this_depth = pd.read_csv(StringIO(this_st_depth), sep='\t', header=None)
        #         this_depth.columns = ['chrom', 'coord', 'cov']
        #         cov_cpg_islands.append(this_depth)
        #     else:
        #         cov_cpg_islands.append(0)
        #
        #     this_perms = []
        #     for j in range(n_perm):
        #         chrom, start, end = genomics.random_genomic_interval(chrom_lengths, n_bp=end - start + 1)
        #         this_st_depth = pysam.depth(fn, "-r", "%s:%d-%d" % (chrom, start, end))
        #         if len(this_st_depth) > 0:
        #             this_depth = pd.read_csv(StringIO(this_st_depth), sep='\t', header=None)
        #             this_depth.columns = ['chrom', 'coord', 'cov']
        #             this_perms.append(this_depth)
        #         else:
        #             this_perms.append(0)
        #     cov_perms.append(this_perms)
        #
        # except Exception as exc:
        #     print "Failed to extract region %s:%d-%d - %s" % (
        #         chrom, start, end, repr(exc)
        #     )

