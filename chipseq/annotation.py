import pandas as pd
import numpy as np
import multiprocessing as mp
import tabix
from . import feature_enrichment
from . import loader
from settings import LOCAL_DATA_DIR
import collections
import os


SOURCES = {'ensembl', 'havana', 'ensembl_havana'}


def assign_peaks_to_basic_features(peak_dat, gtf_fn, tss_pad=500):
    """
    Given the peak data (from MACS), get basic peak assignment: TSS, exon, intron, intergenic
    (in that order of priority)
    :param peak_dat: Pandas DataFrame containing ChIP peaks of interest. Must have the columns start, end and chrom.
    :param gtf_fn: Path to a sorted, BGzipped, tabix-indexed GTF file.
    :return:
    """
    tb = tabix.open(gtf_fn)
    peak_assignment = collections.defaultdict(float)
    for i, row in peak_dat.iterrows():
        qry = tb.query(row.chrom, row.start, row.end)
        hit_gene = False
        hit_exon = False
        hit_tss = False
        for t in qry:
            _, typ, ftr, i0, i1, _, strand, _, attr = t
            if typ in SOURCES:
                if ftr in {'gene', 'transcript'}:
                    hit_gene = True
                if ftr == 'exon':
                    if 'exon_number "1"' in attr:
                        tss_loc = int(i0) if strand == '+' else int(i1)
                        if (row.start - tss_pad) <= tss_loc <= (row.end + tss_pad):
                            hit_tss = True
                        else:
                            hit_exon = True
                    else:
                        hit_exon = True
        if hit_tss:
            peak_assignment['tss'] += 1.
        elif hit_exon:
            peak_assignment['exon'] += 1.
        elif hit_gene:
            peak_assignment['intron'] += 1.
        else:
            peak_assignment['intergenic'] += 1.
    return peak_assignment


if __name__ == '__main__':
    run_type = 'default'
    # recreate the essential operation of Homer suite's annotatePeaks.pl
    peaks = loader.load_macs2_by_patient('all', run_type='default')
    dat = peaks.data.loc[:, 'GBM050H3K27ac'].dropna()

    # 1) overlap of peaks with simple genomic features (intron, exon, tss, intergenic)
    # define TSS as being exon 1 start +/- a padding distance
    # order of priority:
    # 1) TSS (or close to a TSS)
    # 2) Exon
    # 3) Intron
    # 4) Intergenic (none of the above)
    tss_pad = 500

    # in order to use tabix (fast genomic region indexer / looker-upper), we need a sorted, bgzipped GTF file
    fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human',
        'ensembl',
        'GRCh38.release87',
        'gtf',
        'Homo_sapiens.GRCh38.87.sorted.gtf.bgz'
    )
    tb = tabix.open(fn)

    peak_assignment = collections.defaultdict(float)

    chroms = [str(t) for t in range(1, 23)] + ['X', 'Y']
    for c in dat.chrom.unique():
        print "Computing peak assignment for chromosome %s" % c
        this_dat = dat.loc[dat.chrom == c]
        for i, row in this_dat.iterrows():
            qry = tb.query(c, row.start, row.end)
            hit_gene = False
            hit_exon = False
            hit_tss = False
            for t in qry:
                _, typ, ftr, i0, i1, _, strand, _, attr = t
                if typ in SOURCES:
                    if ftr in {'gene', 'transcript'}:
                        hit_gene = True
                    if ftr == 'exon':
                        if 'exon_number "1"' in attr:
                            tss_loc = int(i0) if strand == '+' else int(i1)
                            if (row.start - tss_pad) <= tss_loc <= (row.end + tss_pad):
                                hit_tss = True
                            else:
                                hit_exon = True
                        else:
                            hit_exon = True
            if hit_tss:
                peak_assignment['tss'] += 1.
            elif hit_exon:
                peak_assignment['exon'] += 1.
            elif hit_gene:
                peak_assignment['intron'] += 1.
            else:
                peak_assignment['intergenic'] += 1.



    reg, names = feature_enrichment.get_gene_tss_from_gtf(fn, distance=0)
    closest_tss = {}
    for c in chroms:
        this_tss = [t for t in names if t[0] == c]
        this_tss_loc = np.array([t[1] for t in this_tss])
        closest_tss[c] = {}
        this_dat = dat.loc[dat.chrom == c]
        search_ix = np.searchsorted(this_tss_loc, this_dat.start + this_dat.rel_peak_pos, side='left')
        ix = search_ix
        ix[ix == this_tss_loc.shape[0]] -= 1
        tss_coords = this_tss_loc[ix]
