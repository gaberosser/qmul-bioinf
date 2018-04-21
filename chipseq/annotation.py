import pandas as pd
import numpy as np
import tabix
import collections
import os
import subprocess


def prepare_tabix_indexed_gtf(gtf_fn):
    """
    It is necessary to run this procedure to index a GTF file for use with the tabix application.
    This gtf must be sorted first
    The GTF must first be SORTED, then BGZIPPED, then INDEXED
    """
    filestem, ext = os.path.splitext(gtf_fn)
    if ext.lower() == '.gz':
        cmd = "zcat {fn}"
        outfile = filestem + ".bgz"
    else:
        cmd = "cat {fn}"
        outfile = gtf_fn + ".bgz"
    cmd += " | bedtools sort | bgzip > {outfn}"
    cmd = cmd.format(fn=gtf_fn, outfn=outfile)
    subprocess.call(cmd, shell=True)
    cmd = "tabix {fn}".format(fn=outfile)
    subprocess.call(cmd, shell=True)


def assign_peaks_to_basic_features(peak_dat, gtf_fn, tss_pad=500, sources=None):
    """
    Given the peak data (from MACS), get basic peak assignment: TSS, exon, intron, intergenic
    (in that order of priority)
    :param peak_dat: Pandas DataFrame containing ChIP peaks of interest. Must have the columns start, end and chrom.
    :param gtf_fn: Path to a sorted, BGzipped, tabix-indexed GTF file.
    :param tss_pad: Number of bases (in both directions) considered a hit for the TSS.
    :param sources: If supplied, limit the GTF sources to these
    :return:
    """
    tb = tabix.open(gtf_fn)
    peak_assignment = pd.Series(index=peak_dat.index)
    assignment_count = collections.defaultdict(float)
    for i, row in peak_dat.iterrows():
        qry = tb.query(str(row.chrom), row.start, row.end)
        hit_gene = False
        hit_exon = False
        hit_tss = False
        for t in qry:
            _, typ, ftr, i0, i1, _, strand, _, attr = t
            if (sources is None) or (typ in sources):
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
            assignment_count['tss'] += 1.
            peak_assignment[i] = 'tss'
        elif hit_exon:
            assignment_count['exon'] += 1.
            peak_assignment[i] = 'exon'
        elif hit_gene:
            assignment_count['intron'] += 1.
            peak_assignment[i] = 'intron'
        else:
            assignment_count['intergenic'] += 1.
            peak_assignment[i] = 'intergenic'
    return peak_assignment


def get_closest_tss(peak_loc, tss_loc, tss_names):
    """
    For every peak in peak_loc, find the nearest TSS in tss_loc. Provide the name (from tss_names) and distance.
    Distance is defined so that negative values are upstream of the TSS.
    NB not considering chromosomes here - so split by those in an outer loop.
    :param peak_loc: Numpy array containing the coordinates of the peaks (all in the same chromosome).
    :param tss_loc: Numpy array containing the coordinates of each TSS (all in the matching chromosome).
    :param tss_names: Numpy array containing the name(s) of each TSS (all in the matching chromosome).
    :return: Pandas DataFrame with two columns: 'gene' and 'distance_to_tss'.
    """
    N_pk = len(peak_loc)
    N = len(tss_loc)
    ix = np.searchsorted(tss_loc, peak_loc, side='left')

    closest_tss = pd.DataFrame(index=range(N_pk), columns=['gene', 'distance_to_tss'])

    if N == 0:
        # no peaks to search, so return an empty dataframe
        return closest_tss

    # ix == 0 and ix == N are special cases (closest TSS is unambiguous)
    closest_tss.loc[ix == 0, 'gene'] = tss_names[0]
    closest_tss.loc[ix == N, 'gene'] = tss_names[N - 1]
    closest_tss.loc[ix == 0, 'distance_to_tss'] = peak_loc[ix == 0] - tss_loc[0]
    closest_tss.loc[ix == N, 'distance_to_tss'] = peak_loc[ix == N] - tss_loc[N - 1]

    # everything else: we need to check "one up and one down" for the closest match
    for i in range(1, N):
        n_in_segment = (ix == i).sum()
        dist_down = peak_loc[ix == i] - tss_loc[i - 1]
        dist_up = tss_loc[i] - peak_loc[ix == i]
        winning_ix_rel = np.array([0, 1])[(dist_up < dist_down).astype(int)]
        # winning_ix_abs = np.array([i - 1, i])[(dist_up < dist_down).astype(int)]
        winning_ix_abs = winning_ix_rel + i - 1
        closest_tss.loc[ix == i, 'gene'] = tss_names[winning_ix_abs]
        closest_tss.loc[ix == i, 'distance_to_tss'] = np.vstack((dist_down, -dist_up)).transpose()[
            range(n_in_segment),
            winning_ix_rel
        ]

    return closest_tss
