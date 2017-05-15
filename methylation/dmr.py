from load_data import methylation_array
from methylation.process import m_from_beta, merge_illumina_probe_gene_classes
from methylation import plots
from scipy import ndimage, stats
from stats import nht
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
            probes[j] = p.tolist()
            j += 1
    # return probes, genes
    return probes


def get_bounded_region(reg_coll, anno, chr, loc_from, loc_to):
    """
    :param reg_coll: The output of PyDMR.identify_regions. This is a dict, indexed by chr, with nested dict elements
    indexed by probe class.
    """
    portion = anno.loc[(anno.MAPINFO > loc_from) & (anno.MAPINFO < loc_to) & (anno.CHR == '1')]
    this_regs = {}
    for cls, d in reg_coll[chr].iteritems():
        for k, prs in d.items():
            if len(portion.index.intersection(prs) == len(prs)):
                this_regs.setdefault(cls, []).append(prs)
    return this_regs


def test_region_1vs1(probes, dat=None, min_median_foldchange=1.4):
    """
    For the supplied probe list (which defines a region), run a statistical test to determine whether the median
    difference observed in the data (e.g. M values or beta values) is meaningful. This function is designed to
    compare only the case where there are two samples.
    :param probes: Iterable containing the probe IDs to be tested
    :param dat: The data corresponding to the probes. Columns are samples, rows are probes.
    :param min_median_foldchange:
    :return: P value, unless comparison fails due to insufficient data, in which case None
    """
    y = dat.loc[probes]
    y1 = y.iloc[:, 0].dropna()
    y2 = y.iloc[:, 1].dropna()
    if y1.size != y2.size:
        logger.error("Data are not matched among the two groups")
        return
    if y1.size < 4:
        logger.error("Unable to compute statistics for <4 observations")
        return
    return nht.mannwhitneyu_test(y1, y2)


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

    def test_regions(self, dmr_probes, use_data='beta'):
        if use_data == 'beta':
            data = self.b
        elif use_data == 'm':
            data = self.m
        else:
            raise AttributeError("Unrecognised use_data value.")
        pvals = {}

        if self.n_jobs > 1:
            jobs = {}

        for chr, v1 in dmr_probes.iteritems():
            print chr
            p1 = {}
            dat = get_chr(self.anno, chr)
            for typ, v2 in v1.iteritems():
                p2 = []
                for _, these_probes in v2.iteritems():
                    # TODO
                    pass
        #         if self.n_jobs > 1:
        #             jobs[(chr, cl)] = self.pool.apply_async(identify_region, args=(coords, n_min, d_max))
        #         else:
        #             p2 = identify_region(coords, n_min, d_max)
        #             p1[cl] = p2
        #     dmr_probes[chr] = p1
        #
        #     if self.n_jobs > 1:
        #         # fill in the dict from the deferred results
        #         for (chr, cl), j in jobs.items():
        #             try:
        #                 p2 = j.get(1e3)
        #                 dmr_probes[chr][cl] = p2
        #             except Exception:
        #                 logger.exception("Failed to compute DMR for chr %s class %s", chr, cl)
        #
        # return dmr_probes



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
            p[cls].update(reduce(operator.add, v2.values(), []))
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
                try:
                    dmr_probes = dmr_obj.identify_regions(n, d)
                    record_output(n, d, dmr_probes)
                except Exception:
                    logger.exception("Failed with n_min=%d, d_max=%.1f", n, d)

        if n_jobs > 1:
            for (n, d), j in jobs.items():
                try:
                    dmr_probes = j.get(1e3)
                    record_output(n, d, dmr_probes)
                except Exception:
                    logger.exception("Failed with n_min=%d, d_max=%.1f", n, d)

    finally:
        dmr_obj.n_jobs = n_jobs_orig
        dmr_obj.start_pool()

    return n_regions, n_probes


def plot_n_region_heatmap(dat, n_arr, d_arr, ax=None, **kwargs):
    # heatmap plots for individual classes and combined
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    sns.heatmap(
        dat,
        cmap='RdBu_r',
        xticklabels=d_arr,
        yticklabels=n_arr,
        ax=ax,
        **kwargs
    )
    ax.set_xlabel("Maximum distance between probes")
    ax.set_ylabel("Minimum number of probes")
    ax.figure.tight_layout()

    return ax

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

    if False:

        # carry out a parameter sweep to see how the n_min and d_max parameters affect the number of valid regions
        d_max_arr = np.arange(100, 1100, 100)
        n_min_arr = np.arange(4, 11)
        nreg, npro = dmr_region_parameter_sweep(obj, d_max_arr=d_max_arr, n_min_arr=n_min_arr)

        # combine results for the distinct classes
        nreg_all = reduce(operator.add, nreg.values())
        npro_all = reduce(operator.add, npro.values())
        n_pro_tot = float(anno.shape[0])

        # heatmap plots
        # number of regions

        ax = plot_n_region_heatmap(nreg_all, n_min_arr, d_max_arr)
        cax = [a for a in ax.figure.get_axes() if a is not ax][0]
        cax.set_ylabel('Number of regions considered', rotation=270, labelpad=14)

        ax.figure.savefig(os.path.join(OUTDIR, 'parameter_sweep_regions_total.png'), dpi=200)
        ax.figure.savefig(os.path.join(OUTDIR, 'parameter_sweep_regions_total.pdf'))

        fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True)
        for i, k in enumerate(nreg.keys()):
            plot_n_region_heatmap(nreg[k], n_min_arr, d_max_arr, ax=axs[i], cbar=False)
            axs[i].set_title(k)
            plt.setp(axs[i].xaxis.get_ticklabels(), rotation=90)
            if i == 0:
                plt.setp(axs[i].yaxis.get_ticklabels(), rotation=0)
            else:
                axs[i].yaxis.label.set_visible(False)
        fig.tight_layout()

        fig.savefig(os.path.join(OUTDIR, 'parameter_sweep_regions_by_class.png'), dpi=200)
        fig.savefig(os.path.join(OUTDIR, 'parameter_sweep_regions_by_class.pdf'))

        ax = plot_n_region_heatmap(npro_all / n_pro_tot * 100, n_min_arr, d_max_arr, vmin=0, vmax=100)
        cax = [a for a in ax.figure.get_axes() if a is not ax][0]
        cax.set_ylabel('% probes retained', rotation=270, labelpad=14)

        ax.figure.savefig(os.path.join(OUTDIR, 'parameter_sweep_probes_total.png'), dpi=200)
        ax.figure.savefig(os.path.join(OUTDIR, 'parameter_sweep_probes_total.pdf'))

        fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True)
        for i, k in enumerate(npro.keys()):
            nt = anno.merged_class.str.contains(k).sum()
            plot_n_region_heatmap(npro[k] / nt * 100, n_min_arr, d_max_arr, ax=axs[i], cbar=False, vmin=0, vmax=100)
            axs[i].set_title(k)
            plt.setp(axs[i].xaxis.get_ticklabels(), rotation=90)
            if i == 0:
                plt.setp(axs[i].yaxis.get_ticklabels(), rotation=0)
            else:
                axs[i].yaxis.label.set_visible(False)
        fig.tight_layout()

        fig.savefig(os.path.join(OUTDIR, 'parameter_sweep_probes_by_class.png'), dpi=200)
        fig.savefig(os.path.join(OUTDIR, 'parameter_sweep_probes_by_class.pdf'))

    reg = obj.identify_regions(n_min=n_min, d_max=d_max)

    # plot illustrating probe locations/classes and regions of interest

    loc_from = 1000000
    loc_to = loc_from + 200000

    ax = plots.illustrate_probes(anno, anno.merged_class, '1', loc_from, loc_to, alpha_classed=1.)
    fig = ax.figure
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "chr1_probe_demo.png"), dpi=200)
    fig.savefig(os.path.join(OUTDIR, "chr1_probe_demo.pdf"))

    # add regions
    bounded_reg = get_bounded_region(reg, anno, '1', loc_from, loc_to)
    plots.illustrate_regions(anno, bounded_reg, CLASSES, '1', loc_from, loc_to, ax=ax)


