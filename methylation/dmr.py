from load_data import methylation_array
from methylation.process import m_from_beta, merge_illumina_probe_gene_classes
from methylation.plots import illustrate_probes
from scipy import ndimage
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import operator
import itertools
from functools import partial
from collections import defaultdict
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


class GenomicRegionCollection(object):
    def __init__(self, anno, onetoone_fields=None):
        self.anno = anno

        self.beta = beta
        self.m = m

        self.probes_by_chr = {}
        self.probes_by_class = {}

        # counters

    def add_region_by_probelist(self, chr, cls, probes):
        self.probes_by_chr.setdefault(chr, {}).setdefault(cls, {})

class PyDMR(object):
    def __init__(self, anno, beta, m, n_jobs=-1):
        self.anno = anno
        self.b = beta
        self.m = m
        if n_jobs == -1:
            n_jobs = mp.cpu_count()

        self.n_jobs = n_jobs
        self.start_pool()

        logger.info("Merging classes.")
        self.add_merged_probe_classes()

    def start_pool(self):
        if self.n_jobs > 1:
            logger.info("Creating a pool of %d workers to do my bidding.", self.n_jobs)
            self.pool = mp.Pool(self.n_jobs)

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

    def test_regions(self, dmr_probes):
        pass


def region_count(dmr_probes):
    n = defaultdict(int)
    for chr, v1 in dmr_probes.iteritems():
        for cls, v2 in v1.iteritems():
            n[cls] += len(v2)
    return dict(n)


def probe_count(dmr_probes):
    p = defaultdict(set)
    for chr, v1 in dmr_probes.iteritems():
        for cls, v2 in v1.iteritems():
            p[cls].update(reduce(operator.add, v2.values()))
    return dict([(cls, len(p[cls])) for cls in p])


def dmr_sweep_wrapper(n, d, dmr_obj=None):
    return dmr_obj.identify_regions(n, d)


def dmr_region_parameter_sweep(dmr_obj, n_min_arr, d_max_arr, n_jobs=None):
    n_regions = defaultdict(lambda: np.zeros((len(n_min_arr), len(d_max_arr))))
    n_probes = defaultdict(lambda: np.zeros((len(n_min_arr), len(d_max_arr))))

    n_map = dict([(n, i) for i, n in enumerate(n_min_arr)])
    d_map = dict([(d, i) for i, d in enumerate(d_max_arr)])

    def record_output(n, d, dmr_probes):
        n_reg = region_count(dmr_probes)
        n_pro = probe_count(dmr_probes)
        for cls, r in n_reg.items():
            p = n_pro[cls]
            n_regions[cls][n_map[n], d_map[d]] = r
            n_probes[cls][n_map[n], d_map[d]] = p

    n_jobs_orig = dmr_obj.n_jobs
    if n_jobs is None:
        n_jobs = dmr_obj.n_jobs
    elif n_jobs == -1:
        n_jobs = mp.cpu_count()

    try:
        if n_jobs > 1:
            jobs = {}
            dmr_obj.n_jobs = 1
            dmr_obj.pool.close()
            dmr_obj.pool = None
            pool = mp.Pool(n_jobs)
            func = partial(dmr_sweep_wrapper, dmr_obj=dmr_obj)

        for n, d in itertools.product(n_min_arr, d_max_arr):
            if n_jobs > 1:
                jobs[(n, d)] = pool.apply_async(func, args=(n, d))
            else:
                dmr_probes = dmr_obj.identify_regions(n, d)
                record_output(n, d, dmr_probes)

        if n_jobs > 1:
            for (n, d), j in jobs.items():
                dmr_probes = j.get(1e3)
                record_output(n, d, dmr_probes)

    finally:
        dmr_obj.n_jobs = n_jobs_orig
        dmr_obj.start_pool()

    return n_regions, n_probes


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

    loc_from = 1000000
    loc_to = loc_from + 500000

    ax = illustrate_probes(anno, anno.merged_class, '1', loc_from, loc_to, alpha_classed=1.)
    fig = ax.figure
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "chr1_probe_demo.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "chr1_probe_demo.pdf"))

    # identify connected regions with same class, ignoring strand
    dm = 1.4
    n_min = 4
    d_max = 200

    # from sklearn.cluster import DBSCAN
    #
    # d = DBSCAN(eps=d_max, min_samples=n_min, n_jobs=-1, metric="cityblock")
    # lbl = {}
    # n_lbl = {}
    # for cl in CLASSES:
    #     coords = chr1.loc[chr1.merged_class.str.contains(cl), 'MAPINFO']
    #     coords = coords.values.reshape((coords.size, 1))
    #     fit = d.fit(coords)
    #     lbl[cl] = fit.labels_
    #     n_lbl[cl] = len(set(fit.labels_).difference({-1}))

    obj = PyDMR(anno, b, m)
    # reg = obj.identify_regions(n_min=n_min, d_max=d_max)