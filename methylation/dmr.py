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
import copy
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

CLASSES = {
    'gene',
    'tss',
    'island'
}

TEST_METHODS = {
    'none', 'mwu', 'mwu_permute', 'wsrt', 'wsrt_permute'
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


def get_clusters_by_location(reg_coll, anno, chr, loc_from, loc_to):
    """
    :param reg_coll: The output of identify_clusters. This is a nested dict.
    """
    portion = anno.loc[(anno.MAPINFO > loc_from) & (anno.MAPINFO < loc_to) & (anno.CHR == chr)]
    this_regs = {}
    for cls, d in reg_coll[chr].iteritems():
        for k, prs in d.items():
            # check whether all probes are present in the selected region of the genome
            if len(portion.index.intersection(prs) == len(prs)):
                this_regs.setdefault(cls, {})[k] = prs
    return this_regs


def identify_clusters(anno, n_min=4, d_max=200, n_jobs=1):
    clusters = {}

    if n_jobs > 1:
        jobs = {}
        pool = mp.Pool(processes=n_jobs)

    for chr in pd.factorize(anno.CHR)[1]:
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

    # run back over each chromosome
    # for each cluster, add a list of classes it belongs to
    classes = {}
    for chr, dat in clusters.items():
        t = {}
        for cl in CLASSES:
            for x in dat[cl].values():
                t.setdefault(tuple(x), set()).add(cl)
        classes[chr] = t

    # reform the dictionary, but with meaningful cluster IDs and cluster classes added
    clusters_new = {}
    for chr, dat in clusters.items():
        track_multiples = {}
        clusters_new[chr] = {}
        j = 0
        for cl in CLASSES:
            clusters_new[chr][cl] = {}
            for x in dat[cl].values():
                probe_tup = tuple(x)
                probe_classes = classes[chr][probe_tup]
                if len(probe_classes) == 1:
                    clusters_new[chr][cl][j] = {
                        'probes': x,
                        'classes': probe_classes,
                    }
                    j += 1
                elif probe_tup in track_multiples:
                    # we've already added this cluster under a different class, so we'll just reference it here
                    cid, pr = track_multiples[probe_tup]
                    clusters_new[chr][cl][cid] = pr
                else:
                    # we need to add this cluster and track it
                    clusters_new[chr][cl][j] = {
                        'probes': x,
                        'classes': probe_classes,
                    }
                    track_multiples[probe_tup] = (j, clusters_new[chr][cl][j])
                    j += 1

    return clusters_new


def median_change(y1, y2):
    return np.abs(np.median(y1 - y2))


def reshape_data(dat):
    if len(dat.shape) == 1:
        dat = dat.reshape((dat.size, 1))
    elif len(dat.shape) == 2:
        pass
    else:
        raise ValueError("Invalid data supplied with %d dimensions (must be either 1D or 2D)" % len(dat.shape))
    return dat


def test_cluster_data_values(y1, y2, min_median_change=1.4, method='mwu_permute', test_kwargs=None):
    """
    Runs tests for relevance (min median change) and significance (Mann Whitney U) for the two datasets.
    Mann Whitney is equivalent to the Wilcoxon rank-sum test for independent samples.
    :param y1: Data for group 1. If replicates are included, one column per replicate.
    :param y2: Data for group 2.
    :param min_median_change: The minimum absolute median change required to call a change significant
    :param method: The statistical test applied to determine significant results. Options are:
        None: do not apply a statistical test (relevance check only)
        mwu: Mann-Whitney U statistic, non-parametric analogue to the t test. Supports replicates (fastest).
        wsrt: Wilcoxon Signed Rank Test, non-parametric analogue to the paired t test. Only supports n=1 comparisons.
        mwu_permute: Permutation test with the Mann-Whtiney U statistic. Supports replicates. (default)
        wsrt_permute: Permutation test with the Wilcoxon Signed Rank T statistic. Only supports n=1 comparisons.
    :param test_kwargs: If supplied, these kwargs are passed to the testing function.
    :return:
    """

    if method is not None and method not in TEST_METHODS:
        raise AttributeError("Requested method %s not supported. Supported methods are %s" % (
            method,
            ', '.join(TEST_METHODS)
        ))

    if test_kwargs is None:
        test_kwargs = {}

    if len(y1) != len(y2):
        raise ValueError("Data are not matched among the two groups")

    if len(y1) < 4:
        raise ValueError("Unable to compute statistics for <4 observations")

    y1 = reshape_data(np.asarray(y1))
    n1 = y1.shape[1]
    y2 = reshape_data(np.asarray(y2))
    n2 = y2.shape[1]

    if n1 == 1:
        m1 = y1
    elif n1 > 1:
        m1 = np.median(y1, axis=1)
    else:
        raise ValueError("y1 has an empty dimension (no samples?)")

    if n2 == 1:
        m2 = y2
    elif n2 > 1:
        m2 = np.median(y2, axis=1)
    else:
        raise ValueError("y1 has an empty dimension (no samples?)")

    res = dict(
        median_change=np.median(m1 - m2),
        median1=np.median(m1),
        median2=np.median(m2),
    )

    if method is None or method == 'none':
        # don't run the statistical test
        pass
    else:
        if np.abs(res['median_change']) > min_median_change:
            if method == 'mwu':
                res['pval'] = nht.mannwhitneyu_test(y1.flat, y2.flat, **test_kwargs)
            elif method == 'wsrt':
                if (n1 > 1) or (n2 > 1):
                    raise AttributeError("Wilcoxon signed rank test is not suitable for replicates.")
                res['pval'] = nht.wilcoxon_signed_rank_test(y1, y2, **test_kwargs)
            elif method == 'mwu_permute':
                res['pval'] = wilcoxon_rank_sum_permutation(y1, y2, **test_kwargs)
            elif method == 'wsrt_permute':
                raise NotImplementedError

    return res


def test_cluster_data_values_parallel(probes, **kwargs):
    """
    Parallel version of `test_cluster_data_values`.
    Assumes that the data variable is available (it should be a global property in the pool)
    """
    y1 = pdata[0].loc[probes].dropna()
    y2 = pdata[1].loc[probes].dropna()
    return test_cluster_data_values(y1, y2, **kwargs)


def test_clusters(clusters, data, samples, min_median_change=1.4, n_jobs=1, **kwargs):
    """
    Compare beta or M values between two samples.
    Each cluster is marked as either 'not of interest', 'relevant and non-significant', or
    'relevant AND significant'. Only relevant regions are tested, reducing the number of hypothesis tests and
    increasing the power.
    :param dmr_probes: Output of identify regions
    :param samples: Iterable of length two, containing strings referencing the sample names.
    :param kwargs: Passed to test_cluster_data_values
    :return: Dict with same structure as dmr_probes. Each cluster has an associated dictionary:
        {'relevant': True/False, 'pval': None if not relevant, float if relevant}
    """
    if len(samples) != 2:
        raise AttributeError("samples must have length two and reference columns in the data")

    test_kwargs = dict(kwargs)
    test_kwargs.update({'min_median_change': min_median_change})

    res = {}

    y1 = data.loc[:, samples[0]]
    y2 = data.loc[:, samples[1]]

    if n_jobs > 1:
        pool = mp.Pool(n_jobs, initializer=init_pool, initargs=((y1, y2),))
        jobs = {}

    for chr, d in clusters.iteritems():
        # chromosome loop
        res[chr] = {}
        for typ, cldict in d.iteritems():
            # cluster type loop
            res[chr][typ] = {}
            for j, probedict in cldict.iteritems():
                res[chr][typ][j] = dict(probedict)
                if n_jobs > 1:
                    jobs[(chr, typ, j)] = pool.apply_async(
                        test_cluster_data_values_parallel,
                        args=(probedict['probes'],),
                        kwds=test_kwargs
                    )
                else:
                    try:
                        this_y1 = y1.loc[probedict['probes']].dropna()
                        this_y2 = y2.loc[probedict['probes']].dropna()
                        res[chr][typ][j].update(
                            test_cluster_data_values(this_y1, this_y2, **test_kwargs)
                        )
                    except Exception:
                        logger.error("Failed on chr %s type %s number %s", chr, typ, j)

    if n_jobs > 1:
        # close pool and wait for execution to complete
        for (chr, typ, j), task in jobs.iteritems():
            try:
                res[chr][typ][j].update(task.get(1e3))
            except Exception:
                logger.exception("Failed on chr %s type %s number %s", chr, typ, j)
    pool.close()
    return res


def mht_correction(test_results, alpha=0.05, method='fdr_bh', copy=False):
    """
    :param test_results: As generated by test_clusters.
    :param copy: If True, return a copy of the results. Otherwise (default), modify results in place and return a
    new dictionary referencing them.
    """
    # Filter to include only relevant clusters
    # NB this doesn't copy values, so any modifications are represented in the original dictionary
    test_results_only_relevant = filter_dictionary(test_results, filt=lambda x: 'pval' in x, n_level=3, copy=copy)
    # gather pvalues
    keys, attrs = zip(*dict_iterator(test_results_only_relevant, n_level=3))
    pvals = [t['pval'] for t in attrs]
    rej, padj, _, _ = multicomp.multipletests(pvals, alpha=alpha, method=method)

    # add results to dictionary
    for r, pa, d in zip(rej, padj, attrs):
        d['padj'] = pa
        d['rej_h0'] = r

    return test_results_only_relevant


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
        # if len(k) <= level:
        if len(k) < level:
            raise ValueError("Requested level is too low for this dictionary")
        if k[j] == key:
            remaining_keys = k[:j] + k[(j+1):]
            par = res
            for rk in remaining_keys[:-1]:
                par = par.setdefault(rk, {})
            par[remaining_keys[-1]] = val
    return res


def filter_dictionary(d, filt, n_level=None, copy=False):
    g = dict_iterator(d, n_level=n_level)
    res = {}
    for k, val in g:
        if not filt(val):
            continue
        par = res
        for rk in k[:-1]:
            par = par.setdefault(rk, {})
        par[k[-1]] = val
    if copy:
        res = copy.deepcopy(res)
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
            p[cls].update(reduce(operator.add, [t['probes'] for t in v2.values()], []))
    return dict([(cls, len(p[cls])) for cls in p])


def dmr_sweep_parallel_wrapper(n, d, anno=None):
    return identify_clusters(anno, n_min=n, d_max=d, n_jobs=1)


def dmr_region_parameter_sweep(anno, n_min_arr, d_max_arr, n_jobs=1):
    n_clusters = defaultdict(lambda: np.zeros((len(n_min_arr), len(d_max_arr))))
    n_probes = defaultdict(lambda: np.zeros((len(n_min_arr), len(d_max_arr))))
    clusters = defaultdict(dict)

    n_map = dict([(n, i) for i, n in enumerate(n_min_arr)])
    d_map = dict([(d, i) for i, d in enumerate(d_max_arr)])

    def record_output(n, d, dmr_probes):
        n_reg = region_count(dmr_probes)
        n_pro = probe_count(dmr_probes)
        for cls, r in n_reg.items():
            p = n_pro[cls]
            n_clusters[cls][n_map[n], d_map[d]] = r
            n_probes[cls][n_map[n], d_map[d]] = p
            clusters[cls][(n, d)] = dmr_probes

    if n_jobs > 1:
        jobs = {}
        func = partial(dmr_sweep_parallel_wrapper, anno=anno)
        pool = mp.Pool(n_jobs)
        # pool = mp.Pool(n_jobs, initializer=init_pool, initargs=(anno,))

    for n, d in itertools.product(n_min_arr, d_max_arr):
        if n_jobs > 1:
            # jobs[(n, d)] = pool.apply_async(dmr_sweep_parallel_wrapper, args=(n, d))
            jobs[(n, d)] = pool.apply_async(func, args=(n, d))
        else:
            try:
                dmr_probes = identify_clusters(anno, n_min=n, d_max=d)
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

    return clusters, n_clusters, n_probes


def cluster_permutation_test(dat1, dat2, probes=None, n_probe=None, n_perm=1000):
    """
    Permute observed data and determine the null distribution of the observed median value change
    Permutation is carried out so that alike probes are being considered.
    :param dat1:
    :param dat2:
    :param probes:
    :param n_perm:
    :return:
    """
    if n_probe is None and probes is None:
        raise AttributeError("Must supply either probes or n_probes")
    if n_probe is not None and probes is not None:
        raise AttributeError("Must supply either probes or n_probes")
    if n_probe is None:
        n_probe = len(probes)
    m = dat1.size
    medians = []
    for i in range(n_perm):
        idx = np.random.choice(m, n_probe, replace=True)
        t1 = dat1[idx]; t2 = dat2[idx]
        medians.append((t2 - t1).median())
    if probes is not None:
        obs = (dat2[probes] - dat1[probes]).median()
        return obs, np.array(sorted(medians))
    else:
        return np.array(sorted(medians))


def cluster_confined_permutation_test(dat1, dat2, anno, probes=None, n_probe=None, n_perm=1000):
    """
    Permute observed data and determine the null distribution of the observed median value change
    Permutation is carried out so that probes must be consecutive
    :param dat1:
    :param dat2:
    :param anno: Required to sort probes by location along the genome
    :param probes:
    :param n_perm:
    :return:
    """
    if n_probe is None and probes is None:
        raise AttributeError("Must supply either probes or n_probes")
    if n_probe is not None and probes is not None:
        raise AttributeError("Must supply either probes or n_probes")
    if n_probe is None:
        n_probe = len(probes)

    m = dat1.size

    # sort by chromosome then genomic coordinate
    # we don't need chromosomes to be in any particular order
    anno_sorted = anno.sort_values(by=['CHR', 'MAPINFO'])
    sort_idx = anno_sorted.index
    dat1 = dat1[sort_idx]
    dat2 = dat2[sort_idx]

    # find the indices of the changeover points for chromosomes. We need to discard any results that span these.
    # the index gives the last entry for each chromosome
    last_idx = np.where(np.diff(anno_sorted.CHR.factorize()[0]) != 0)[0]

    medians = []
    i_max = m - n_probe
    count = 0
    while count < n_perm:
    # for i in range(n_perm):
        i0 = np.random.randint(i_max)
        i1 = i0 + n_probe
        # check for overlap with breakpoint
        if ((i0 <= last_idx) & (i1 > last_idx)).any():
            continue
        t1 = dat1[i0:i1]; t2 = dat2[i0:i1]
        med = (t2 - t1).median()
        if np.isnan(med).any():
            import ipdb; ipdb.set_trace()
        medians.append((t2 - t1).median())
        count += 1

    if probes is not None:
        obs = (dat2[probes] - dat1[probes]).median()
        return obs, np.array(sorted(medians))
    else:
        return np.array(sorted(medians))


def wilcoxon_signed_rank_permutation(dat1, dat2, n_max=9999, return_stats=False):
    assert dat1.size == dat2.size, "Incompatible data sizes"

    d_diff = dat1 - dat2

    n = dat1.size
    n_it = 2 ** n
    b_sample = n_it > n_max
    if b_sample:
        n_it = n_max
        stat = np.zeros(n_it)
        print "Sampling strategy with %d iterations" % n_it
        multipliers = stats.binom.rvs(1, 0.5, size=(n, n_it)) * -2 + 1
        perm_d = np.tile(d_diff, (n_it, 1)).transpose() * multipliers
        for i in range(n_it):
            stat[i] = nht.wilcoxon_signed_rank_statistic(perm_d[:, i])
    else:
        stat = np.zeros(n_it)
        print "Exact strategy with %d iterations" % n_it
        for i in range(n_it):
            str_fmt = ("{0:0%db}" % n).format(i)
            cc = (np.array(list(str_fmt)) == '1') * -2 + 1
            stat[i] = nht.wilcoxon_signed_rank_statistic(d_diff * cc)

    t = stats.wilcoxon(dat1, dat2)
    p = (stat <= t.statistic).sum() / float(n_it)
    if return_stats:
        return (p, stat)
    return p


def wilcoxon_rank_sum_permutation(x, y, n_max=1999, return_stats=False):
    """
    :param n_max: The maximum number of iterations to run. If the exact test requires more than this, we revert to
    sampling. Setting this to None, 0 or negative forces exact sampling, but this might be very slow and expensive.
    """
    force_exact = (n_max is None) or (n_max <= 0)
    x = reshape_data(np.asarray(x))
    y = reshape_data(np.asarray(y))

    nx = x.shape[1]
    ny = y.shape[1]

    k = x.shape[0]

    if k != y.shape[0]:
        raise ValueError("The number of probes in the two samples does not match.")

    n_it1 = 2 ** k
    n_it2 = nx * ny
    n_it = n_it1 * n_it2

    if not force_exact and n_it > n_max:
        n_it = n_max
        stat = np.zeros(n_it)

        multipliers = stats.binom.rvs(1, 0.5, size=(k, n_it))
        jxs = np.random.randint(nx, size=n_it)
        jys = np.random.randint(ny, size=n_it)
        for i in range(n_it):
            perm_x = x.copy()
            perm_y = y.copy()
            idx1 = multipliers[:, i] == 1
            jx = jxs[i]
            jy = jys[i]
            # perform the data swap
            perm_x[idx1, jx] = y[idx1, jy]
            perm_y[idx1, jy] = x[idx1, jx]
            stat[i] = nht.mannwhitneyu_statistic(perm_x.flat, perm_y.flat)

    else:
        stat = np.zeros(n_it)
        count = 0
        for i in range(n_it1):
            str_fmt = ("{0:0%db}" % k).format(i)
            idx1 = np.array(list(str_fmt)) == '1'
            for jx in range(nx):
                for jy in range(ny):
                    perm_x = x.copy()
                    perm_y = y.copy()
                    # perform the data swap
                    perm_x[idx1, jx] = y[idx1, jy]
                    perm_y[idx1, jy] = x[idx1, jx]
                    stat[count] = nht.mannwhitneyu_statistic(perm_x.flat, perm_y.flat)
                    count += 1

    u = nht.mannwhitneyu_statistic(x.flat, y.flat)
    p = (stat <= u).sum() / float(n_it)

    if return_stats:
        return (p, stat)
    return p

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


if __name__ == "__main__":
    OUTDIR = unique_output_dir('dmr', reuse_empty=True)

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

    fdr = 0.05
    dm = 1.4
    n_min = 6
    d_max = 400
    n_jobs = mp.cpu_count()

    clusters = identify_clusters(anno, n_min=n_min, d_max=d_max, n_jobs=n_jobs)
    n_clusters = len(list(
        dict_iterator(clusters)
    ))

    test_results = {}
    test_results_relevant = {}
    test_results_significant = {}
    for sid in ['018', '019', '031']:
        samples = ('GBM%s' % sid, 'Dura%s' % sid)
        test_results[sid] = test_clusters(clusters, m, samples=samples, min_median_change=dm, n_jobs=n_jobs)
        test_results_relevant[sid] = mht_correction(test_results[sid], alpha=fdr)

    for sid in ['018', '019', '031']:
        test_results_significant[sid] = filter_dictionary(
            test_results_relevant[sid],
            filt=lambda x: x['rej_h0'],
            n_level=3
        )
        n_relevant = len(list(
            dict_iterator(test_results_relevant[sid], n_level=3)
        ))
        n_significant = len(list(
            dict_iterator(
                test_results_significant[sid],
                n_level=3
            )
        ))
        print "GBM%s vs Dura%s: %d clusters, %d relevant (%f %%), %d significant (%f %%)" % (
            sid, sid,
            n_clusters,
            n_relevant,
            n_relevant / float(n_clusters) * 100.,
            n_significant,
            n_significant / float(n_clusters) * 100.
        )

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
        bounded_reg = get_clusters_by_location(clusters, anno, '1', loc_from, loc_to)
        plots.illustrate_regions(anno, bounded_reg, CLASSES, '1', loc_from, loc_to, ax=ax)

    # multipanel plot showing probe locations/classes, regions of potential interest and methylation data for the
    # same probes
    if False:
        the_chr = '1'
        loc_from = 1000000
        loc_to = loc_from + 200000

        for sid in ['018', '019', '031']:
            fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(12, 6))

            plots.illustrate_probes(anno, anno.merged_class, the_chr, loc_from, loc_to, alpha_classed=1., ax=axs[0])
            bounded_reg = get_clusters_by_location(clusters, anno, the_chr, loc_from, loc_to)
            plots.illustrate_regions(anno, bounded_reg, CLASSES, ax=axs[0])

            portion = anno.loc[(anno.MAPINFO > loc_from) & (anno.MAPINFO < loc_to) & (anno.CHR == the_chr)]
            samples = ['Dura%s' % sid, 'GBM%s' % sid]

            the_data = m.loc[portion.index, samples]
            the_locs = portion.MAPINFO.values

            axs[1].scatter(the_locs, the_data.loc[:, samples[0]], facecolor='none', edgecolor='b', marker='o', label=samples[0])
            axs[1].scatter(the_locs, the_data.loc[:, samples[1]], facecolor='none', edgecolor='r', marker='o', label=samples[1])
            axs[1].legend(loc='upper left', frameon=True, facecolor='w')

            # link the test results back to the probes for plotting
            # do this separately for relevant and significant clusters
            ymin = the_data.min().min()
            ymax = the_data.max().max()
            height = ymax - ymin

            # pick out significant and relevant clusters for plotting
            sign_for_plot = {}
            rele_for_plot = {}
            for (kcls, kid), pids in dict_iterator(bounded_reg):
                the_res = test_results[sid][the_chr][kcls][kid]
                if 'pval' not in the_res:
                    continue
                elif the_res['rej_h0']:
                    sign_for_plot.setdefault(kcls, {})[kid] = pids
                else:
                    rele_for_plot.setdefault(kcls, {})[kid] = pids

            plots.illustrate_regions(anno, sign_for_plot, CLASSES, ax=axs[-1], ylim=(ymin, ymax), linestyle='-')
            plots.illustrate_regions(anno, rele_for_plot, CLASSES, ax=axs[-1], ylim=(ymin, ymax), linestyle='--')

            axs[0].xaxis.label.set_visible(False)
            axs[0].set_ylabel('Strand')
            axs[1].set_xlabel('Chromosomal coordinate')
            axs[1].set_ylabel('M value')
            fig.tight_layout()

            fig.savefig(os.path.join(OUTDIR, "chr1_%s_demo_with_Mvalues.png" % sid), dpi=200)
            fig.savefig(os.path.join(OUTDIR, "chr1_%s_demo_with_Mvalues.pdf" % sid))
