from load_data import methylation_array
from methylation.process import m_from_beta, beta_from_m, merge_illumina_probe_gene_classes
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
from statsmodels.sandbox.stats import multicomp

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


def init_pool(_data):
    global pdata
    pdata = _data


def add_merged_probe_classes(anno):
    anno.loc[:, 'merged_class'] = merge_illumina_probe_gene_classes(
        anno.loc[:, 'UCSC_RefGene_Group'], anno.loc[:, 'Relation_to_UCSC_CpG_Island']
    )


def identify_cluster(coords, n_min, d_max):
    probes = {}
    di = np.diff(coords.values.flat) <= d_max
    ll, nn = ndimage.measurements.label(di)
    j = 0
    for i in range(1, nn + 1):
        this_idx = np.where(ll == i)[0]
        if this_idx.size + 1 >= n_min:
            p = coords.iloc[np.arange(this_idx[0], this_idx[-1] + 2)].index
            probes[j] = p.tolist()
            j += 1
    return probes


def get_bounded_region(reg_coll, anno, chr, loc_from, loc_to):
    """
    :param reg_coll: The output of identify_clusters. This is a nested dict.
    """
    portion = anno.loc[(anno.MAPINFO > loc_from) & (anno.MAPINFO < loc_to) & (anno.CHR == chr)]
    this_regs = {}
    for cls, d in reg_coll[chr].iteritems():
        for k, prs in d.items():
            # check whether all probes are present in the selected region of the genome
            if len(portion.index.intersection(prs) == len(prs)):
                this_regs.setdefault(cls, []).append(prs)
    return this_regs


def identify_clusters(anno, n_min=4, d_max=200, n_jobs=1):
    clusters = {}

    if n_jobs > 1:
        jobs = {}
        pool = mp.Pool(processes=n_jobs)

    for chr in pd.factorize(anno.CHR)[1]:
        print chr
        p1 = {}
        dat = get_chr(anno, chr)
        for cl in CLASSES:
            coords = get_coords(dat, cl)
            if n_jobs > 1:
                jobs[(chr, cl)] = pool.apply_async(identify_cluster, args=(coords, n_min, d_max))
            else:
                p2 = identify_cluster(coords, n_min, d_max)
                p1[cl] = p2
        clusters[chr] = p1

        if n_jobs > 1:
            # fill in the dict from the deferred results
            for (chr, cl), j in jobs.items():
                try:
                    p2 = j.get(1e3)
                    clusters[chr][cl] = p2
                except Exception:
                    logger.exception("Failed to compute DMR for chr %s class %s", chr, cl)
    return clusters


def median_change(y1, y2):
    return np.abs(np.median(y1 - y2))


def test_cluster_data_values(y1, y2, min_median_change=1.4):
    """
    Runs tests for relevance (min median change) and significance (Mann Whitney U) for the two datasets.
    :param y1:
    :param y2:
    :param min_median_change: The minimum absolute median change required to call a change significant
    :return:
    """
    if len(y1) != len(y2):
        raise ValueError("Data are not matched among the two groups")

    if len(y1) < 4:
        raise ValueError("Unable to compute statistics for <4 observations")

    res = {'median_change': median_change(y1, y2)}
    if res['median_change'] > min_median_change:
        res['pval'] = nht.mannwhitneyu_test(y1, y2)

    return res


def test_cluster_data_values_parallel(probes, **kwargs):
    """
    Parallel version of `test_cluster_data_values`.
    Assumes that the data variable is available (it should be a global property in the pool)
    """
    y1 = pdata[0].loc[probes].dropna()
    y2 = pdata[1].loc[probes].dropna()
    return test_cluster_data_values(y1, y2, **kwargs)


def test_clusters(clusters, data, samples, min_median_change=1.4, n_jobs=1):
    """
    Compare beta or M values between two samples.
    Each cluster is marked as either 'not of interest', 'relevant and non-significant', or
    'relevant AND significant'. Only relevant regions are tested, reducing the number of hypothesis tests and
    increasing the power.
    :param dmr_probes: Output of identify regions
    :param samples: Iterable of length two, containing strings referencing the sample names.
    :param use_data: String specifying which data to use for comparison,
    :return: Dict with same structure as dmr_probes. Each cluster has an associated dictionary:
        {'relevant': True/False, 'pval': None if not relevant, float if relevant}
    """
    if len(samples) != 2:
        raise AttributeError("samples must have length two and reference columns in the data")

    res = {}

    y1 = data.loc[:, samples[0]]
    y2 = data.loc[:, samples[1]]

    if n_jobs > 1:
        pool = mp.Pool(n_jobs, initializer=init_pool, initargs=((y1, y2),))
        jobs = {}

    for chr, d in clusters.iteritems():
        # chromosome loop
        logger.info("Chromosome %s", chr)
        res[chr] = {}
        for typ, probedict in d.iteritems():
            # cluster type loop
            res[chr][typ] = {}
            for j, these_probes in probedict.iteritems():
                if n_jobs > 1:
                    jobs[(chr, typ, j)] = pool.apply_async(
                        test_cluster_data_values_parallel,
                        args=(these_probes,),
                        kwds={'min_median_change': min_median_change}
                    )
                else:
                    try:
                        this_y1 = y1.loc[these_probes].dropna()
                        this_y2 = y2.loc[these_probes].dropna()
                        res[chr][typ][j] = test_cluster_data_values(this_y1, this_y2,
                                                                    min_median_change=min_median_change)
                    except Exception:
                        logger.error("Failed on chr %s type %s number %s", chr, typ, j)

    if n_jobs > 1:
        # close pool and wait for execution to complete
        for (chr, typ, j), task in jobs.iteritems():
            try:
                res[chr][typ][j] = task.get(1e3)
            except Exception:
                logger.exception("Failed on chr %s type %s number %s", chr, typ, j)
    pool.close()
    return res



class ProbeCluster(object):
    def __init__(self, pids, anno, cls=None, chr=None, beta=None, m=None):
        """
        Represents a collection of methylation probes
        :param pids: The IDs of the probes belonging to this cluster
        :param anno: The annotation data for this probe, from the Illumina manifest. Probes are in rows.
        May include other columns (e.g. `merged_class`).
        :param beta, m: If supplied, these are the beta and M data corresponding to these probes, in a pandas DataFrame.
        Columns represent samples, rows represent probes.
        """
        self.anno = anno
        self.pids = pids
        self.beta = beta
        self.m = m
        self.n_probes = len(pids)

        if chr is None:
            self.chr = self.anno.CHR[0]
            if (self.anno.CHR != self.chr).any():
                raise AttributeError("Probes are located on different chromosomes")
        else:
            self.chr = chr

        self.cls = cls

        self.pval = None
        self.padj = None

    @property
    def coord_list(self):
        return self.anno.MAPINFO

    @property
    def coord_range(self):
        cs = self.coord_list
        return (min(cs), max(cs))

    def check_inputs(self, samples, use_data):
        if use_data == 'beta':
            data = self.beta
        elif use_data == 'm':
            data = self.m
        else:
            raise AttributeError("Unsupported data")
        if data is None:
            raise AttributeError("Cannot continue without any data")

        if samples is None:
            if data.shape[1] > 2:
                raise ValueError("If the number of samples is >2, must specify which two samples to compare")
            else:
                samples = data.columns
        return samples, data

    def relevant(self, min_median_diff, samples=None, use_data='beta'):
        samples, data = self.check_inputs(samples, use_data)
        y1 = data.loc[:, samples[0]]; y2 = data.loc[:, samples[1]]

        return np.abs(np.median(y1 - y2)) > min_median_diff

    def compute_pvalue(self, samples=None, use_data='beta'):
        if self.n_probes < 4:
            logger.error("Unable to compute statistics for <4 observations")
            return
        samples, data = self.check_inputs(samples, use_data)
        y1 = data.loc[:, samples[0]]; y2 = data.loc[:, samples[1]]
        self.pval = nht.mannwhitneyu_test(y1, y2)
        return self.pval



def dict_iterator(d, parents=None, level=1, n_level=None):
    if parents is None:
        parents = []
    if (n_level is not None) and (level > n_level):
        yield (parents, d)
    else:
        for k in d.iterkeys():
            if isinstance(d[k], dict):
               for x in dict_iterator(d[k], parents + [k], level=level+1, n_level=n_level):
                   yield x
            else:
                yield (parents + [k], d[k])


def dict_by_sublevel(d, level, key, n_level=None):
    """
    Given the sublevel with specified key, rebuild the dictionary
    :param d:
    :param level:
    :param key:
    :param n_level:
    :return:
    """
    # translate level for zero indexing
    j = level - 1
    g = dict_iterator(d, n_level=n_level)
    res = {}
    for k, val in g:
        if len(k) <= level:
            raise ValueError("Requested level is too low for this dictionary")
        if k[j] == key:
            remaining_keys = k[:j] + k[(j+1):]
            par = res
            for rk in remaining_keys[:-1]:
                par = par.setdefault(rk, {})
            par[remaining_keys[-1]] = val
    return res


def filter_dictionary(d, filt, n_level=None):
    g = dict_iterator(d, n_level=n_level)
    res = {}
    for k, val in g:
        if not filt(val):
            continue
        par = res
        for rk in k[:-1]:
            par = par.setdefault(rk, {})
        par[k[-1]] = val
    return res


class ClusterCollection(object):
    """
    Class for managing a collection of clusters.
    This allows indexing by chromosome or class type and adds methods for efficient bulk computation of p values.
    """
    def __init__(self, cluster_dict=None, pval_dict=None):
        """
        If supplied, clusters are stored in a dictionary.
        First level: chromosome
        Second level: probe class
        Third level: list of clusters
        The same hierarchy is used in `self._clusters` and `self._pvals`
        :param cluster_dict:
        :param pval_dict:
        """
        self._clusters = cluster_dict or {}


        self.clusters_by_chromosome = {}
        self.clusters_by_class = {}

        self._pvals = pval_dict or {}
        self.pvals_by_chromosome = {}
        self.pvals_by_class = {}

    # TODO?


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


def dmr_sweep_wrapper(n, d, anno=None):
    return identify_clusters(anno, n_min=n, d_max=d, n_jobs=1)


def dmr_region_parameter_sweep(anno, n_min_arr, d_max_arr, n_jobs=1):
    # FIXME: the parallel version could be MUCH faster if we declare global anno data and use that rather than
    # passing anno each time
    n_clusters = defaultdict(lambda: np.zeros((len(n_min_arr), len(d_max_arr))))
    n_probes = defaultdict(lambda: np.zeros((len(n_min_arr), len(d_max_arr))))

    n_map = dict([(n, i) for i, n in enumerate(n_min_arr)])
    d_map = dict([(d, i) for i, d in enumerate(d_max_arr)])

    def record_output(n, d, dmr_probes):
        n_reg = region_count(dmr_probes)
        n_pro = probe_count(dmr_probes)
        for cls, r in n_reg.items():
            p = n_pro[cls]
            n_clusters[cls][n_map[n], d_map[d]] = r
            n_probes[cls][n_map[n], d_map[d]] = p

    if n_jobs > 1:
        jobs = {}
        pool = mp.Pool(n_jobs)
        func = partial(dmr_sweep_wrapper, anno=anno)

    for n, d in itertools.product(n_min_arr, d_max_arr):
        if n_jobs > 1:
            jobs[(n, d)] = pool.apply_async(func, args=(n, d))
        else:
            try:
                dmr_probes = func(n, d)
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

    return n_clusters, n_probes


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

    # reduce anno down to probes in the data
    anno = anno.loc[anno.index.intersection(b.index)]

    # add merged class column to annotation
    add_merged_probe_classes(anno)

    # split gene symbols and store as a set
    anno.loc[:, 'UCSC_RefGene_Name'] = \
        anno.UCSC_RefGene_Name.str.split(';').apply(lambda x: set(x) if isinstance(x, list) else None)

    # identify connected regions with same class, ignoring strand
    dm = 1.4
    n_min = 4
    d_max = 200
    n_jobs = 12

    # slower way to assess clustering: but we might want to use this in future?
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


    clusters = identify_clusters(anno, n_min=n_min, d_max=d_max, n_jobs=n_jobs)
    test_results = test_clusters(clusters, m, samples=('GBM019', 'Dura019'), min_median_change=dm, n_jobs=n_jobs)
    test_results_only_relevant = filter_dictionary(test_results, filt=lambda x: 'pval' in x, n_level=3)

    # correct for multiple hypothesis testing
    # TODO: use statsmodels to test


    if False:

        # carry out a parameter sweep to see how the n_min and d_max parameters affect the number of valid regions
        d_max_arr = np.arange(100, 1100, 100)
        n_min_arr = np.arange(4, 11)
        nreg, npro = dmr_region_parameter_sweep(anno, d_max_arr=d_max_arr, n_min_arr=n_min_arr, n_jobs=n_jobs)

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

    # plot illustrating probe locations/classes and regions of interest
    if False:
        loc_from = 1000000
        loc_to = loc_from + 200000

        ax = plots.illustrate_probes(anno, anno.merged_class, '1', loc_from, loc_to, alpha_classed=1.)
        fig = ax.figure
        fig.tight_layout()
        fig.savefig(os.path.join(OUTDIR, "chr1_probe_demo.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "chr1_probe_demo.pdf"))

        # add regions
        bounded_reg = get_bounded_region(clusters, anno, '1', loc_from, loc_to)
        plots.illustrate_regions(anno, bounded_reg, CLASSES, '1', loc_from, loc_to, ax=ax)


