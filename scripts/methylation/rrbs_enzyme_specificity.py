import os
import collections
import gzip
import numpy as np
from scipy import stats
import re
from settings import DATA_DIR, LOCAL_DATA_DIR, GIT_LFS_DATA_DIR
import pysam
from matplotlib import pyplot as plt
import pandas as pd
import multiprocessing as mp
from utils import log, genomics, output
logger = log.get_console_logger(__name__)


if __name__ == "__main__":
    outdir = output.unique_output_dir("rrbs_enzyme_specificity", reuse_empty=True)

    basedir = os.path.join(DATA_DIR, 'rrbseq', 'GC-CV-7163')

    indir = os.path.join(basedir, 'trim_galore_mouse/bismark')
    bam_fn = os.path.join(indir, 'GC-CV-7163-6_S6_pe.sorted.bam')
    cov_fn = os.path.join(indir, 'GC-CV-7163-6_S6_bismark.cov.gz')
    s = pysam.AlignmentFile(bam_fn, 'rb')
    chroms = [str(t) for t in range(1, 20)]

    # theoretical (binomial) distribution of inferred methylation by coverage
    Ns = [10, 20, 50, 100]
    cs = ['k', 'r', 'b', 'g']
    ps = [0.1, 0.25, 0.5]
    for p in ps:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for N, c in zip(Ns, cs):
            x = np.arange(N)
            y = stats.binom.pmf(x, N, p)
            ax.plot(x / float(N), y * len(x), '%s-o' % c, label='N = %d' % N)
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.set_xlabel("Inferred methylation level")
        ax.set_ylabel("Normalised probability density")
        ax.legend(loc='upper right')
        fig.savefig(os.path.join(outdir, "binomial_p%.2f.png" % p), dpi=200)

    reads = {}
    for c in chroms:
        reads[c] = []
        it = s.fetch(c)
        for seg in it:
            if seg.flag == 99 or seg.flag == 163:
                # 99: properly mapped in pair, first in pair, forward strand
                # 163: properly mapped in pair, second in pair, forward strand
                if seg.tlen > 0:
                    reads[c].append(seg)

    # read.get_forward_sequence() returns the ORIGINAL sequence, as it appears in the fq file.
    # read.seq returns the ALIGNED sequence
    # this will be the reverse complement if the read mapped to the reverse strand
    # this will be the same if the read mapped to the forward strand

    # Here's an educational example
    c = '19'
    rd_a = reads[c][0]
    rd_b = s.mate(rd_a)
    if rd_a.is_read1:
        rd1 = rd_a
        rd2 = rd_b
    else:
        rd1 = rd_b
        rd2 = rd_a

    print "Chromosome %s" % c
    print "Read 1: %s" % rd1.qname
    print "Read 2: %s" % rd2.qname

    if rd1.is_reverse and not rd2.is_reverse:
        print "R1 is reversed, R2 is forward (F2R1)."
    elif not rd1.is_reverse and rd2.is_reverse:
        print "R1 is forward, R2 is reverse (F2R1)."
    else:
        print "R1 and R2 have the same direction (%s) (probably bad?)" % ('reverse' if rd1.is_reverse else 'forward')

    # get the sequences from the reference
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

    # now we just want to check the equivalence of an aligned read and the reference sequence
    # TODO

    # for each fragment: does it contain a CCGG site?
    contains_ccgg = []
    contains_ccgg_template = []
    for c in chroms:
        this_ref = fa_reader[c]
        print "Chromosome %s" % c
        for rd in reads[c]:
            if rd.is_proper_pair and not rd.is_reverse:
                ref_seq = this_ref[rd.reference_start:rd.reference_start + rd.template_length]
                if 'CCGG' in ref_seq:
                    contains_ccgg.append(rd)
                    contains_ccgg_template.append(ref_seq)
                # print "F read seq: %s" % rd.seq
                # print "Ref seq: %s" % this_ref[rd.reference_start:rd.reference_start + rd.reference_length]
                # print "*********"

    nc = sum([len(v) for v in reads.values()])
    print "%d / %d fragments contain CCGG (%.2f%%)" % (
        len(contains_ccgg),
        nc,
        len(contains_ccgg) / float(nc) * 100
    )

    chrom_ct = collections.Counter()
    coords = {}
    for rd in contains_ccgg:
        chrom_ct[rd.reference_name] += 1
        coords.setdefault(rd.reference_name, []).append(rd.reference_start)
    all_coords = reduce(lambda x, y: x + y, coords.values())

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

    # for each theoretical fragment, compute the size and coverage
    mspi_fragments = {}
    for c in chroms:
        print "CCGG fragment coverage. Chrom %s" % c
        this_coords = ccgg_coords[c]
        mspi_fragments[c] = []
        this_res = mspi_fragments[c]
        for i in range(1, len(ccgg_coords[c])):
            t0 = this_coords[i - 1]
            t1 = this_coords[i]
            this_res.append([t1 - t0, s.count(c, t0, t1)])

    size_cov = np.concatenate(mspi_fragments.values(), axis=0)




    # get coverage of every CpG
    cpg_cov = {}
    for c in chroms:
        print "Chromosome %s" % c
        cov = s.count_coverage(c)
        cova = [cov[0][i] for i in cpg_coords[c]]
        covc = [cov[1][i] for i in cpg_coords[c]]
        covg = [cov[2][i] for i in cpg_coords[c]]
        covt = [cov[3][i] for i in cpg_coords[c]]
        sa = sum(cova)
        sc = sum(covc)
        sg = sum(covg)
        st = sum(covt)
        pr = (sa + sg) / float(sa + sg + sc + st)
        if pr > 0.01:
            print "Warning: chromosome %s has a high proportion of A and G bases where we expect C or T (%.2f)" % (
                c, pr
            )
        cpg_cov[c] = covc + covt

    cpg_cov_all_nz = np.array(reduce(lambda x, y: x + y, [[t for t in x if t > 0] for x in cpg_cov.values()]))

    # make an inverse CDF (is that called a CEDF?)
    # cc = np.sort(cpg_cov_all_nz)[::-1]
    cov = []
    ecdf = []
    for x in np.unique(cpg_cov_all_nz):
        cov.append(x)
        ecdf.append((cpg_cov_all_nz >= x).sum())
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(cov[:25], ecdf[:25])
    ax.set_xticks(cov[:25])
    ax.set_xlabel('Minimum coverage')
    ax.set_ylabel('Number of CpG sites')

    # inverse CDF of low coverage region

    fig = plt.figure(figsize=(8.5, 5))
    ax1 = fig.add_subplot(111)
    ax1.bar(cov[:25], np.array(ecdf[:25]) / float(n_cpg) * 100)
    ax1.set_xticks(cov[:25])
    ax1.set_xlabel('Minimum coverage')
    ax1.set_ylabel('% CpG sites')
    ax2 = ax1.twinx()
    h = ax2.plot(cov[:25], np.array(ecdf[:25]) / 1e6, 'x')
    ax2.set_ylim(np.array(ax1.get_ylim()) / 100 * n_cpg / 1e6)
    ax2.set_ylabel("Number of CpG sites (millions)")
    h[0].set_visible(False)
    fig.tight_layout()

    fig.savefig(os.path.join(outdir, "rrbs_cpg_coverage_low.png"), dpi=200)

    # inverse CDF of higher coverage
    cov_vals = np.array([20, 30, 40, 50, 60, 70, 80, 90, 100, 200])
    ecdf_vals = np.array([ecdf[cov.index(t)] for t in cov_vals])

    fig = plt.figure(figsize=(8.5, 5))
    ax1 = fig.add_subplot(111)
    ax1.bar(range(len(cov_vals)), ecdf_vals / float(n_cpg) * 100)
    ax1.set_xticks(range(len(cov_vals)))
    ax1.set_xticklabels(cov_vals)
    ax1.set_xlabel('Minimum coverage')
    ax1.set_ylabel('% CpG sites')
    ax2 = ax1.twinx()
    h = ax2.plot(range(len(cov_vals)), ecdf_vals / 1e3, 'x')
    ax2.set_ylim(np.array(ax1.get_ylim()) / 100 * n_cpg / 1e3)
    ax2.set_ylabel("Number of CpG sites (thousands)")
    h[0].set_visible(False)
    fig.tight_layout()

    fig.savefig(os.path.join(outdir, "rrbs_cpg_coverage_high.png"), dpi=200)

    cpg_tsv_fn = os.path.join(GIT_LFS_DATA_DIR, 'mouse_cpg_island', 'grcm38_cpgisland.tsv')

    # load tsv
    cpg_regions = pd.read_csv(cpg_tsv_fn, sep='\t', header=0)
    region_pad = 2000

    # get methylation levels in CpG sites in promoters and CpG islands
    promoters = genomics.get_promoter_regions(tax_id=10090)
    promoters = [p for p in promoters if p['chr'] in chroms]
    promoters_by_region = {}
    for p in promoters:
        promoters_by_region.setdefault((p['chr'], p['strand']), collections.defaultdict(list))
        promoters_by_region[(p['chr'], p['strand'])][(p['promoter_region_start'], p['promoter_region_end'])].append(p)
        # promoters_by_region.setdefault(p['chr'], collections.defaultdict(list))
        # promoters_by_region[p['chr']][(p['promoter_region_start'], p['promoter_region_end'])].append(p)

    # get CpG sites associated with the unique promoter regions
    cpg_coords_promoter_regions = {}
    for c in chroms:
        # strand doesn't matter here, but we need to check both
        the_cpg_coords = np.array(cpg_coords[c])
        cpg_coords_promoter_regions.setdefault(c, [])
        for strand in ['+', '-']:
            if (c, strand) in promoters_by_region:
                for region in promoters_by_region[(c, strand)]:
                    reg_min = min(region)
                    reg_max = max(region)
                    in_region_idx = (the_cpg_coords >= reg_min) & (the_cpg_coords <= reg_max)
                    cpg_coords_promoter_regions[c].extend(the_cpg_coords[in_region_idx])

    # read the coverage file, splitting into chromosomes
    cov = {}
    raw = pd.read_csv(cov_fn, sep='\t', header=None, index_col=None)
    raw.columns = ['chrom', 'coord', 'coord1', 'pct_meth', 'n_meth', 'n_unmeth']
    raw.loc[:, 'chrom'] = raw.loc[:, 'chrom'].astype(str)
    for c in chroms:
        this_raw = raw.loc[raw.chrom == c]
        cov[c] = this_raw.drop('chrom', axis=1)

    # methylation and coverage of CpG island regions
    methylation_in_cpg_islands = []
    methylation_in_promoters = []

    for c in chroms:
        this_data = cpg_regions.loc[cpg_regions.chrom == c, ['chromStart', 'chromEnd', 'length', 'cpgNum']]
        for _, r in this_data.iterrows():
            idx = (cov[c].loc[:, 'coord'] >= r.chromStart) & (cov[c].loc[:, 'coord'] <= r.chromEnd)
            new_dat = cov[c].loc[idx, ['coord', 'n_meth', 'n_unmeth', 'pct_meth']]
            new_dat.insert(0, 'coverage', new_dat.n_meth + new_dat.n_unmeth)
            new_dat.insert(0, 'beta', new_dat.n_meth / new_dat.coverage)
            new_dat.drop('pct_meth', axis=1, inplace=True)
            new_dat.insert(0, 'chrom', c)
            new_dat.insert(0, 'length', r.length)
            new_dat.insert(0, 'n_cpg', r.cpgNum)
            methylation_in_cpg_islands.append(new_dat)

        this_data = promoters_by_region[(c, '+')].keys() + promoters_by_region[(c, '-')].keys()
        for start, end in this_data:
            idx = (cov[c].loc[:, 'coord'] >= start) & (cov[c].loc[:, 'coord'] <= end)
            new_dat = cov[c].loc[idx, ['coord', 'n_meth', 'n_unmeth', 'pct_meth']]
            new_dat.insert(0, 'coverage', new_dat.n_meth + new_dat.n_unmeth)
            new_dat.insert(0, 'beta', new_dat.n_meth / new_dat.coverage)
            new_dat.drop('pct_meth', axis=1, inplace=True)
            new_dat.insert(0, 'chrom', c)
            methylation_in_promoters.append(new_dat)

    methylation_in_cpg_islands = pd.concat(methylation_in_cpg_islands, axis=0)
    methylation_in_promoters = pd.concat(methylation_in_promoters, axis=0)
