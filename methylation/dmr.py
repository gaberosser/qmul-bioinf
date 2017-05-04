from load_data import methylation_array
from methylation.process import m_from_beta, merge_illumina_probe_gene_classes
from scipy import ndimage
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
from utils.output import unique_output_dir
import numpy as np
import logging
import multiprocessing as mp
logger = logging.getLogger(__name__)
if len(logger.handlers) == 0:
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.DEBUG)

NORM_METHOD = 'swan'
OUTDIR = unique_output_dir('dmr', reuse_empty=True)
CLASSES = {
    'gene',
    'tss',
    'island'
}


def get_chr(dat, chr, col='CHR'):
    return dat.loc[dat.loc[:, col] == chr].sort_values(by='MAPINFO', axis=0)


def get_coords(dat, gene_class, col='merged_class'):
    return dat.loc[dat.loc[:, col].str.contains(gene_class), 'MAPINFO']


def identify_region(coords, n_min, d_max):
    # genes = {}
    probes = {}
    # coords = get_coords(dat, cl)
    di = np.diff(coords.values.flat) <= d_max
    ll, nn = ndimage.measurements.label(di)
    j = 0
    for i in range(1, nn + 1):
        this_idx = np.where(ll == i)[0]
        if this_idx.size + 1 >= n_min:
            p = coords.iloc[np.arange(this_idx[0], this_idx[-1] + 2)].index
            # locs = dat.loc[probes, 'MAPINFO']
            # g = dat.loc[p, 'UCSC_RefGene_Name']
            probes[j] = p.tolist()
            # genes[j] = g
            j += 1
    # return probes, genes
    return probes


class PyDMR(object):
    def __init__(self, anno, beta, m, n_jobs=-1):
        self.anno = anno
        self.b = beta
        self.m = m
        if n_jobs == -1:
            n_jobs = mp.cpu_count()

        self.n_jobs = n_jobs
        if n_jobs > 1:
            logger.info("Creating a pool of %d workers to do my bidding.", n_jobs)
            self.pool = mp.Pool(self.n_jobs)

        logger.info("Merging classes.")
        self.add_merged_probe_classes()

    def identify_regions(self, n_min=4, d_max=200):

        dmr_probes = {}

        if self.n_jobs > 1:
            jobs = {}

        for chr in pd.factorize(self.anno.CHR)[1]:
            print chr
            p1 = {}
            dat = get_chr(self.anno, chr)
            for cl in CLASSES:
                coords = get_coords(dat, cl)
                if self.n_jobs > 1:
                    jobs[(chr, cl)] = self.pool.apply_async(identify_region, args=(coords, n_min, d_max))
                else:
                    p2 = identify_region(coords, n_min, d_max)
                    p1[cl] = p2
            dmr_probes[chr] = p1

            if self.n_jobs > 1:
                # fill in the dict from the deferred results
                for (chr, cl), j in jobs.items():
                    try:
                        p2 = j.get(1e3)
                        dmr_probes[chr][cl] = p2
                    except Exception:
                        logger.exception("Failed to compute DMR for chr %s class %s", chr, cl)

        return dmr_probes

    def add_merged_probe_classes(self):
        self.anno.loc[:, 'merged_class'] = merge_illumina_probe_gene_classes(
            self.anno.loc[:, 'UCSC_RefGene_Group'], self.anno.loc[:, 'Relation_to_UCSC_CpG_Island']
        )

    def parameter_sweep(self, n_min_arr, d_max_arr):
        pass

    def test_regions(self, dmr_probes):
        pass


if __name__ == "__main__":

    anno = methylation_array.load_illumina_methylationepic_annotation()
    b = methylation_array.gbm_nsc_methylationepic('swan')
    m = m_from_beta(b)
    anno.loc[:, 'merged_class'] = merge_illumina_probe_gene_classes(
        anno.loc[:, 'UCSC_RefGene_Group'], anno.loc[:, 'Relation_to_UCSC_CpG_Island']
    )

    # split gene symbols and store as a set
    anno.loc[:, 'UCSC_RefGene_Name'] = \
        anno.UCSC_RefGene_Name.str.split(';').apply(lambda x: set(x) if isinstance(x, list) else None)

    # plot them
    chr1 = anno.loc[anno.CHR == '1'].sort_values(by='MAPINFO', axis=0)

    loc_from = 1000000
    loc_to = loc_from + 500000
    idx = (loc_from <= chr1.MAPINFO) & (chr1.MAPINFO < loc_to)

    strand_y = chr1.Strand.where(chr1.Strand == 'F', 0)
    strand_y.loc[strand_y == 'F'] = 1
    strand_y = strand_y.astype(int)

    col_mapper = {
        'tss': [1, 0, 0],
        'gene': [0, 1, 0],
        'island': [0, 0, 1]
    }
    colour_no_role = [0.5, 0.5, 0.5]

    fig = plt.figure(figsize=(10, 3))
    ax = fig.add_subplot(111)

    t = chr1.loc[idx]
    for grp, c in col_mapper.items():
        x = t.MAPINFO[idx & t.merged_class.str.contains(grp)].values
        y = strand_y[idx & t.merged_class.str.contains(grp)].values
        ax.scatter(x, y, c=c, marker='|', label=grp, alpha=0.7, zorder=3)

    x = t.MAPINFO[idx & (t.merged_class == '')].values
    y = strand_y[idx & (t.merged_class == '')].values
    ax.scatter(x, y, c=colour_no_role, marker='|', alpha=0.3, label='none', zorder=2)
    ax.plot([loc_from, loc_to], [0, 0], 'k-', zorder=1, alpha=0.3)
    ax.plot([loc_from, loc_to], [1, 1], 'k-', zorder=1, alpha=0.3)
    ax.legend(loc='center right', frameon=True, facecolor='w')
    ax.yaxis.set_ticks([0, 1])
    ax.yaxis.set_ticklabels(['R', 'F'])
    ax.set_xlabel('Chromosomal coordinate')
    ax.set_title('Chromosome 1: {:,} - {:,}'.format(loc_from, loc_to))
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "chr1_probe_demo.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "chr1_probe_demo.pdf"))

    # identify connected regions with same class, ignoring strand
    dm = 1.4
    n_min = 4
    d_max = 200

    from sklearn.cluster import DBSCAN

    d = DBSCAN(eps=d_max, min_samples=n_min, n_jobs=-1, metric="cityblock")
    lbl = {}
    n_lbl = {}
    for cl in CLASSES:
        coords = chr1.loc[chr1.merged_class.str.contains(cl), 'MAPINFO']
        coords = coords.values.reshape((coords.size, 1))
        fit = d.fit(coords)
        lbl[cl] = fit.labels_
        n_lbl[cl] = len(set(fit.labels_).difference({-1}))


    from bx import intervals
    import HTSeq

    gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)

    dmr_genes = {}
    dmr_probes = {}

    for chr in pd.factorize(anno.CHR)[1]:
        print chr
        g1 = {}
        p1 = {}
        dat = get_chr(anno, chr)
        for cl in CLASSES:
            g2 = {}
            p2 = {}
            coords = get_coords(dat, cl)
            di = np.diff(coords.values.flat) <= d_max
            ll, nn = ndimage.measurements.label(di)
            j = 0
            for i in range(1, nn + 1):
                this_idx = np.where(ll == i)[0]
                if this_idx.size + 1 >= n_min:
                    probes = coords.iloc[np.arange(this_idx[0], this_idx[-1] + 2)].index
                    locs = dat.loc[probes, 'MAPINFO']
                    genes = dat.loc[probes, 'UCSC_RefGene_Name']
                    # ivl = HTSeq.GenomicInterval(chr, locs.min(), locs.max())
                    # gas[ivl] += tuple(genes)
                    p2[j] = probes.tolist()
                    g2[j] = genes
                    j += 1
            p1[cl] = p2
            g1[cl] = g2
        dmr_genes[chr] = g1
        dmr_probes[chr] = p1

