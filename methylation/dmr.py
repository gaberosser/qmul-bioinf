from load_data import methylation_array
from methylation.process import m_from_beta, beta_from_m, merge_illumina_probe_gene_classes
from methylation import plots
from scipy import ndimage, stats
from stats import nht
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import pickle
import types
import operator
import itertools
from functools import partial
from collections import defaultdict
import os
import json
import copy
from utils.output import unique_output_dir
import numpy as np
import logging
import multiprocessing as mp
import types
from statsmodels.sandbox.stats import multicomp
from utils.dictionary import dict_by_sublevel, dict_iterator, filter_dictionary, nested_dict_to_flat, flat_dict_to_nested

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


class BasicLogicException(Exception):
    pass


class TestResultEncoder(json.JSONEncoder):
    """
    Handles saving results that may contain Bool or sets
    """
    def default(self, o):
        if isinstance(o, set):
            return tuple(o)
        elif isinstance(o, bool) or isinstance(o, np.bool_):
            return int(o)
        return super(TestResultEncoder, self).default(o)


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
            if len(portion.index.intersection(prs['probes'])) == len(prs['probes']):
                this_regs.setdefault(cls, {})[k] = prs['probes']
    return this_regs


def identify_clusters(anno, n_min=4, d_max=200, n_jobs=1, **kwargs):
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
    # we use a tuple list of probe IDs as the hash
    classes = {}
    for chr, dat in clusters.items():
        t = {}
        for cl in CLASSES:
            for x in dat[cl].values():
                t.setdefault(tuple(x), set()).add(cl)
        classes[chr] = t

    # reform the dictionary, but with meaningful cluster IDs and cluster classes added
    # since each cluster may appear multiple times (once per class), keep track of what we've seen so we can
    # skip multiples
    res = []
    clusters_seen = set()
    # keep a running count of cluster ID
    cl_id = 1
    for chr, dat in clusters.items():
        for cl in CLASSES:
            for x in dat[cl].values():
                probe_tup = tuple(x)
                if probe_tup in clusters_seen:
                    # nothing to do
                    continue
                res.append(ProbeCluster(
                    probe_tup, anno, cls=classes[chr][probe_tup], chr=chr, cluster_id=cl_id
                ))
                cl_id += 1
                clusters_seen.add(probe_tup)

    return res


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
    The test is run as y1 - y2.
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


def test_clusters_in_place(clusters, data, samples, min_median_change=1.4, n_jobs=1, add_nested_classes=True, **kwargs):
    """
    Compare beta or M values between two samples.
    Each cluster is marked as either 'not of interest', 'relevant and non-significant', or
    'relevant AND significant'. Only relevant regions are tested, reducing the number of hypothesis tests and
    increasing the power.
    :param dmr_probes: Output of identify regions
    :param samples: Iterable of length two, containing strings referencing the sample names.
    :param kwargs: Passed to test_cluster_data_values
    :param add_nested_classes: If True (default), we add dict keys for convenience to access individual classes.
    :return: Dict with same structure as dmr_probes. Each cluster has an associated dictionary:
        {'relevant': True/False, 'pval': None if not relevant, float if relevant}
    """
    if len(samples) != 2:
        raise AttributeError("samples must have length two and reference columns in the data")

    test_kwargs = dict(kwargs)
    test_kwargs.update({'min_median_change': min_median_change})

    res = {}
    res_all = {}

    y1 = data.loc[:, samples[0]]
    y2 = data.loc[:, samples[1]]

    if n_jobs > 1:
        pool = mp.Pool(n_jobs, initializer=init_pool, initargs=((y1, y2),))
        jobs = {}

    for cl in clusters:
        pids = cl.pids
        if n_jobs > 1:
            jobs[cl] = pool.apply_async(
                test_cluster_data_values_parallel,
                args=(pids,),
                kwds=test_kwargs
            )
            res_all[cl.cluster_id] = cl.summary_dict
        else:
            try:
                this_y1 = y1.loc[pids].dropna()
                this_y2 = y2.loc[pids].dropna()
                res_all[cl.cluster_id] = test_cluster_data_values(this_y1, this_y2, **test_kwargs)
                res_all[cl.cluster_id].update(cl.summary_dict)
            except Exception:
                logger.exception("Failed on chr %s class %s number %d", str(cl.chr), str(cl.cls), cl.cluster_id)

    if n_jobs > 1:
        # close pool and wait for execution to complete
        pool.close()
        for cl, task in jobs.iteritems():
            try:
                res_all[cl.cluster_id].update(task.get(1e3))
            except Exception as exc:
                logger.exception(
                    "Failed on chr %s class %s number %d: %s",
                    str(cl.chr),
                    str(cl.cls),
                    cl.cluster_id,
                )

    if add_nested_classes:
    # for convenience, provide nested dictionaries allowing access to single class type
        return cluster_class_nests(res_all)
    else:
        return res_all


def cluster_class_nests(res_all):
    """
    Given the dictionary of DMR results, containing all classes, create a new nested structure that allows access to
    all results (key: 'all') but also to nested dictionaries containing only a given class (key: 'island', etc.)
    """
    classes = set()
    for x in res_all.values():
        classes.update(x['classes'])

    res = {'all': res_all}
    for cls in classes:
        res[cls] = dict([
            (k, v) for k, v in res_all.iteritems() if cls in v['classes']
        ])

    return res


def test_clusters(clusters, data, samples, min_median_change=1.4, n_jobs=1, mht='fdr_bh', alpha=0.05, **kwargs):
    """
    Compare beta or M values between two samples.
    Each cluster is marked as either 'not of interest', 'relevant and non-significant', or
    'relevant AND significant'. Only relevant regions are tested, reducing the number of hypothesis tests and
    increasing the power.
    :param clusters: Dictionary keyed by cluster ID, values are ProbeCluster objects, as generated by identify_clusters.
    :param samples: Iterable of length two, containing strings referencing the sample names.
    :param min_median_change: Minimum median change in data values between the groups to be declared relevant.
     Statistical tests are not carried out if this threshold is not met.
    :param mht: If supplied, apply MHT correction after pvalues have been computed. Disable by setting None.
    :param alpha: Only used if mht is True. `alpha` value used as a cutoff on the adjusted p value.
    :param kwargs: Passed to test_cluster_data_values
    :return: Dict with same keys supplied in clusters dict.dmr_probes. Each cluster has an associated dictionary:
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

    for cl_id, cl in clusters.iteritems():
        pids = cl.pids
        if n_jobs > 1:
            jobs[cl] = pool.apply_async(
                test_cluster_data_values_parallel,
                args=(pids,),
                kwds=test_kwargs
            )
        else:
            try:
                this_y1 = y1.loc[pids].dropna()
                this_y2 = y2.loc[pids].dropna()
                res[cl.cluster_id] = test_cluster_data_values(this_y1, this_y2, **test_kwargs)
            except Exception:
                logger.exception("Failed on chr %s class %s number %d", str(cl.chr), str(cl.cls), cl.cluster_id)

    if n_jobs > 1:
        # close pool and wait for execution to complete
        pool.close()
        pool.join()
        for cl, task in jobs.iteritems():
            try:
                res[cl.cluster_id] = task.get(1e3)
            except Exception:
                logger.exception("Failed on chr %s class %s number %d", str(cl.chr), str(cl.cls), cl.cluster_id)


    if mht:
        mht_correction(res, alpha=alpha, method=mht)

    return res


def mht_correction(test_results, alpha=0.05, method='fdr_bh'):
    """
    Apply a multiple hypothesis testing correction to the test_results. All modifying is done in-place
    :param test_results: As generated by test_clusters.
    """
    # we'll work with the full set of results as other representations (per class) are just pointers
    # the_dict = test_results['all']
    # assume that the supplied dictionary represents all classes
    the_dict = test_results

    # gather pvalues
    # keep a list of keys to ensure order is maintained
    keys = []
    pvals = []
    for k, v in the_dict.iteritems():
        if 'pval' in v:
            pvals.append(v['pval'])
            keys.append(k)

    if len(pvals) == 0:
        logger.warn("No DMR results have associated pvalues. Possibly no relevant comparisons?")
        return

    rej, padj, _, _ = multicomp.multipletests(pvals, alpha=alpha, method=method)

    # add results back to results dictionary
    for k, pa, r in zip(keys, padj, rej):
        the_dict[k].update({
            'padj': pa,
            'rej_h0': r
        })

    # we've done all the modifying in-place, so no return value


class ProbeCluster(object):
    def __init__(
            self,
            pids,
            anno,
            cluster_id=None,
            cls=None,
            chr=None,
            anno_requires_lookup=True,
            compute_genes=True
    ):
        """
        Represents a collection of methylation probes
        :param pids: The IDs of the probes belonging to this cluster
        :param anno: The full annotation data from the Illumina manifest. Probes are in rows.
        May include other columns (e.g. `merged_class`).
        We start by reducing this to just the relevant rows, unless anno_requires_lookup is False.
        :param cls: May be either a string (for a single class), an iterable (for multiple classes) or None (no class).
        :param beta, m: If supplied, these are the beta and M data corresponding to these probes, in a pandas DataFrame.
        :param anno_requires_lookup: If True (default), we start by performing a lookup for relevant rows in anno.
        :param compute_genes: If True (default), compute genes at the start.
        Columns represent samples, rows represent probes.
        """
        self.pids = list(pids)
        self.n_probes = len(pids)
        self.chr = chr
        self.cluster_id = cluster_id

        if hasattr(cls, '__iter__'):
            self.cls = set(cls)
        elif cls is not None:
            self.cls = {cls}
        else:
            self.cls = cls

        if anno_requires_lookup:
            self.anno = anno.loc[self.pids]

        else:
            self.anno = anno

        if chr is None:
            # no chromosome supplied: try to get this fromthe probe IDs
            matching_chr = list(set(self.anno.loc[:, 'CHR'].values))
            if len(matching_chr) > 1:
                raise AttributeError("Probes are located on different chromosomes")
            self.chr = matching_chr[0]

        self._genes = None

        if compute_genes:
            _ = self.genes

    @classmethod
    def from_dict(cls, summary_dict, anno):
        return cls(
            summary_dict['pids'],
            anno,
            cluster_id=summary_dict['cid'],
            cls=summary_dict['classes'],
            chr=summary_dict['chr']
        )

    @property
    def summary_dict(self):
        return {
            'pids': self.pids,
            'genes': self.genes,
            'cls': self.cls,
            'chr': self.chr,
            'cid': self.cluster_id,
        }

    def __getitem__(self, item):
        return self.summary_dict[item]

    def __hash__(self):
        return self.hash().__hash__()

    def hash(self):
        return tuple(self.pids)

    @property
    def genes(self):
        if self._genes is not None:
            return self._genes
        genes = self.anno.UCSC_RefGene_Name.dropna()
        if len(genes) == 0:
            self._genes = set()
        else:
            self._genes = reduce(set.union, genes.tolist())
        return self._genes

    @property
    def coord_list(self):
        return self.anno.MAPINFO

    @property
    def coord_range(self):
        cs = self.coord_list
        return (min(cs), max(cs))


def test_results_to_table(dat, exclude_clusters_with_no_genes=True):
    """
    Convert the nested dictionary structure to a flat table.
    :param dat: Nested dictionary, as output by test_clusters. Assumed structure:
    pid, cluster class (or 'all'), cluster id
    We only need to use 'all', as this has all the cluster classes annotated anyway
    """
    classes = list(CLASSES)
    cols = ['cluster_id', 'chr', 'genes'] + \
           ['class_%s' % t for t in classes] + \
           ['median_1', 'median_2', 'median_delta', 'padj']
    dtypes = ['uint16', 'object', 'object'] + \
             ['bool'] * len(classes) + \
             ['float', 'float', 'float', 'float']
    rows = {}
    seen = {}

    the_dict = dict_by_sublevel(dat, 2, 'all')
    for (pid, cluster_id), attrs in dict_iterator(the_dict, n_level=2):

        if (len(attrs['genes']) == 0) and exclude_clusters_with_no_genes:
            # don't add this DMR with no associated genes
            continue

        # check whether we have already processed this result
        seen.setdefault(pid, set())
        this_hsh = (cluster_id, chr)
        if this_hsh in seen[pid]:
            continue
        else:
            seen[pid].add(this_hsh)
        this_row = [cluster_id, chr, tuple(attrs['genes'])] + \
                   [t in attrs['classes'] for t in classes] + \
                   [attrs['median1'], attrs['median2'], attrs['median_change'], attrs['padj']]
        rows.setdefault(pid, []).append(this_row)

    # run through the pids and generate a pandas DataFrame.
    # additionally run a quick sanity check that no rows are duplicated
    tbl = {}
    for pid, r in rows.iteritems():
        tbl[pid] = pd.DataFrame(r, columns=cols).astype(dict(zip(cols, dtypes)))
        if tbl[pid].duplicated().any():
            raise BasicLogicException("Duplicate rows found with PID %s" % pid)

    return tbl


class DmrResults(object):
    """
    Class for managing a collection of DMR results.
    This allows indexing by chromosome or class type and adds methods for efficient bulk computation of p values.
    """
    def __init__(self, clusters=None, anno=None):
        """
        :param anno: Dataframe containing annotation data for the methylation array.
        :param clusters: List of proposed clusters, as computed by `identify_clusters()`. Either a list, in which case
        cluster IDs are automatically generated, or a dictionary with keys as IDs.
        """
        self.clusters = None
        self.set_clusters(clusters)

        self.results = None

        self.anno = anno
        self.min_median_change = None

        self._clusters_by_class = None
        self._results_by_class = None
        self._results_relevant = None
        self._results_relevant_by_class = None
        self._results_significant = None
        self._results_significant_by_class = None

        self._classes = None

    def copy(self):
        # copy attributes across; only results need to be deep copied, others can be shared
        # skip hidden attributes as these can be recomputed quickly
        new = self.__class__()
        new.anno = self.anno
        new.clusters = self.clusters
        for k, v in self.__dict__.iteritems():
            if k in {'anno', 'clusters'}:
                continue
            if k[0] == '_':
                continue
            setattr(new, k, copy.deepcopy(v))
        return new

    def set_clusters(self, clusters):
        if isinstance(clusters, dict):
            self.clusters = clusters
        elif hasattr(clusters, '__iter__'):
            self.clusters = dict([
                (cl.cluster_id, cl) for cl in clusters
            ])
        elif clusters is None:
            self.clusters = None
        else:
            raise AttributeError ("Unrecognised data format for clusters.")

    def identify_clusters(self, d_max, n_min, n_jobs=1, anno=None, **kwargs):
        if anno is None:
            if self.anno is None:
                raise ValueError("Must supply annotation data either at object construction or to this function.")
            else:
                anno = self.anno
        else:
            if self.anno is None:
                self.anno = anno
        cl = identify_clusters(anno, n_min=n_min, d_max=d_max, n_jobs=n_jobs)
        self.set_clusters(cl)

    def test_clusters(self, data, samples, min_median_change=1.4, n_jobs=1, **kwargs):
        self.min_median_change = min_median_change
        self.results = test_clusters(
            self.clusters,
            data,
            samples,
            min_median_change=min_median_change,
            n_jobs=n_jobs,
            **kwargs
        )

        self._results_by_class = None
        self._results_relevant = None
        self._results_relevant_by_class = None
        self._results_significant = None
        self._results_significant_by_class = None

    def to_table(self, include='all', skip_geneless=True):
        """
        Convert results and cluster attributes to a flat table
        :param include: Controls which clusters are used {'all', 'relevant', 'significant'}
        :param skip_geneless: If True (default), do not include probe clusters that have no associated genes
        :return:
        """
        # convert classes to a list to ensure reproducible iteration order
        classes = []
        if self._classes is not None:
            classes = list(self.classes)
        cols = ['cluster_id', 'chr', 'genes'] + \
               ['class_%s' % t for t in classes] + \
               ['median_1', 'median_2', 'median_delta', 'padj']
        dtypes = ['uint16', 'object', 'object'] + \
                 ['bool'] * len(classes) + \
                 ['float', 'float', 'float', 'float']
        rows = []

        # set the cluster dictionary used for iteration based on `include`
        if include == 'all':
            cl_dict = self.clusters
        elif include == 'relevant':
            cl_dict = self.clusters_relevant
        elif include == 'significant':
            cl_dict = self.clusters_significant
        else:
            raise AttributeError("Unrecognised include option. Supported: 'all', 'relevant', 'significant'.")

        for cluster_id, cl in cl_dict.items():
            # don't add if this DMR has no associated genes
            if len(cl.genes) == 0 and skip_geneless:
                continue
            the_res = self.results[cluster_id]
            this_row = [cluster_id, cl.chr, tuple(cl.genes)] + \
                       [t in cl.cls for t in classes] + \
                       [the_res['median1'], the_res['median2'], the_res['median_change'], the_res.get('padj', None)]
            rows.append(this_row)

        # run through the pids and generate a pandas DataFrame.
        # additionally run a quick sanity check that no rows are duplicated
        tbl = pd.DataFrame(rows, columns=cols).astype(dict(zip(cols, dtypes))).set_index('cluster_id')
        if tbl.duplicated().any():
            raise BasicLogicException("Duplicate rows found.")

        return tbl

    @property
    def classes(self):
        """
        Compute a set of all classes in the clusters
        """
        if self._classes is None:
            classes = set()
            for x in self.clusters.values():
                classes.update(x.cls)
            self._classes = classes
        return self._classes

    def separate_by_class(self, lookup, convert=None):
        if convert is None:
            convert = lookup
        dest = {}
        for cls in self.classes:
            dest[cls] = dict([t for t in lookup.items() if cls in convert[t[0]].cls])
        return dest

    def clusters_by_class(self, cls=None):
        if self._clusters_by_class is None:
            self._clusters_by_class = self.separate_by_class(self.clusters)
        if cls is None:
            return self._clusters_by_class
        return self._clusters_by_class[cls]

    def results_by_class(self, cls=None):
        if self._results_by_class is None:
            self._results_by_class = self.separate_by_class(self.results, convert=self.clusters)
        if cls is None:
            return self._results_by_class
        return self._results_by_class[cls]

    @property
    def clusters_relevant(self):
        return dict([
            (k, self.clusters[k]) for k in self.results_relevant
        ])

    @property
    def clusters_significant(self):
        return dict([
            (k, self.clusters[k]) for k in self.results_significant
        ])

    @property
    def results_relevant(self):
        if self._results_relevant is None:
            self._results_relevant = filter_dictionary(self.results, lambda x: 'pval' in x, n_level=1)
        return self._results_relevant

    @property
    def results_significant(self):
        if self._results_significant is None:
            self._results_significant = filter_dictionary(self.results, lambda x: x.get('rej_h0', False), n_level=1)
        return self._results_significant

    def results_relevant_by_class(self, cls=None):
        if self._results_relevant_by_class is None:
            self._results_relevant_by_class = self.separate_by_class(self.results_relevant, convert=self.clusters)
        if cls is None:
            return self._results_relevant_by_class
        return self._results_relevant_by_class[cls]

    def results_significant_by_class(self, cls=None):
        if self._results_significant_by_class is None:
            self._results_significant_by_class = self.separate_by_class(self.results_significant, convert=self.clusters)
        if cls is None:
            return self._results_significant_by_class
        return self._results_significant_by_class[cls]


class DmrResultCollection(object):
    """
    Represents a collection of DmrResult objects. All objects must share the same clusters and anno attributes,
    but can differ in the results dictionary.
    """

    def __init__(self, **objs):
        """
        :param objs: The DmrResults objects. These can be bundled into a nested dictionary structure
        """
        if len(objs) == 0:
            raise ValueError("Must supply some `DmrResults` objects.")

        # flat representation of supplied (possibly nested) objects
        self._flat = nested_dict_to_flat(objs)

        # check objects for compatibility
        first = True
        clusters = None
        anno = None

        for x in self._flat.values():
            if first:
                clusters = x.clusters
                anno = x.anno
                first = False
            else:
                if x.clusters is not clusters:
                    raise AttributeError("All DmrResults objects must share the same `clusters` attribute.")
                if x.anno is not anno:
                    raise AttributeError("All DmrResults objects must share the same `anno` attribute.")

        self.objects = objs
        self.clusters = clusters
        self.anno = anno

    def __getitem__(self, item):
        """
        Act as a direct lookup into objects
        """
        return self.objects[item]

    def __iter__(self):
        return self.objects.__iter__()

    def keys(self):
        return self.objects.keys()

    def iterkeys(self):
        return self.objects.iterkeys()

    def apply(self, func, *args, **kwargs):
        # TODO: test thoroughly
        if isinstance(func, str):
            s = func
            def func(x):
                if hasattr(x.__class__, s):
                    inst = getattr(x.__class__, s)
                    # not data, but could be a property
                    if isinstance(inst, property):
                        return getattr(x, s)
                    elif isinstance(inst, types.MethodType):
                        return getattr(x, s)(*args, **kwargs)
                    else:
                        # what is this? a static data object?
                        raise NotImplementedError()
                else:
                    # attribute or function
                    inst = getattr(x, s)
                    if not hasattr(inst, '__call__'):
                        return inst
                    else:
                        return inst(*args, **kwargs)

        flat_apply = dict([
            (k, func(v)) for k, v in self._flat.items()
        ])
        return flat_dict_to_nested(flat_apply)

    def clusters_by_class(self, cls=None):
        return self.apply(lambda x: x.clusters_by_class(cls=cls))

    def results_by_class(self, cls=None):
        return self.apply(lambda x: x.results_by_class(cls=cls))

    @property
    def clusters_relevant(self):
        return self.apply('clusters_relevant')

    @property
    def clusters_significant(self):
        return self.apply('clusters_significant')

    @property
    def results_relevant(self):
        return self.apply('results_relevant')

    @property
    def results_significant(self):
        return self.apply('results_significant')

    def results_relevant_by_class(self, cls=None):
        return self.apply(lambda x: x.results_relevant_by_class(cls=cls))

    def results_significant_by_class(self, cls=None):
        return self.apply(lambda x: x.results_significant_by_class(cls=cls))

    # FIXME: JSON components disabled for now, because using tuples as a dict key breaks JSON support

    # def to_json(self, fn):
    #     """
    #     Save this collection to a JSON file. We don't store the annotation, as it can be reloaded from elsewhere.
    #     """
    #     # extract results from flat dict
    #     results = dict([(k, v.results) for k, v in self._flat.items()])
    #     # convert clusters to a serializable form
    #     clusters = dict([
    #         (k, v.summary_dict) for k, v in self.clusters.items()
    #     ])
    #     x = dict(clusters=clusters, results=results)
    #     with open(fn, 'wb') as f:
    #         json.dump(x, f, cls=TestResultEncoder)
    #
    # @classmethod
    # def from_json(cls, fn, anno):
    #     """
    #     Load from a json file previously saved with to_json
    #     """
    #     with open(fn, 'rb') as f:
    #         x = json.load(f)
    #     # recreate clusters
    #     clusters = dict([
    #         (k, ProbeCluster.from_dict(v, anno)) for k, v in x['clusters'].items()
    #     ])
    #
    #     # recreate objects
    #     objects = {}
    #     for k, v in x['results'].items():
    #         objects[k] = DmrResults(clusters=clusters, anno=anno)
    #         objects[k].results = v
    #
    #     return cls(**flat_dict_to_nested(objects))

    def to_pickle(self, fn, include_annotation=True):
        """
        Save this collection to pickle. We can store the annotation there, too
        """
        # extract results from flat dict
        results = dict([(k, v.results) for k, v in self._flat.items()])
        x = dict(clusters=self.clusters, results=results)
        if include_annotation:
            x['anno'] = self.anno
        with open(fn, 'wb') as f:
            pickle.dump(x, f)

    @classmethod
    def from_pickle(cls, fn, anno=None):
        """
        Load from pickle.
        :param anno: must be included if it is not bundled
        """
        with open(fn, 'rb') as f:
            x = pickle.load(f)

        if 'anno' in x:
            if anno is not None:
                logger.warning("Supplied annotation data will override bundled annotation data.")
            else:
                anno = x['anno']
        elif anno is None:
            raise AttributeError("Must supply annotation data if not bundled")

        # recreate objects
        objects = {}
        for k, v in x['results'].items():
            objects[k] = DmrResults(clusters=x['clusters'], anno=anno)
            objects[k].results = v

        return cls(**flat_dict_to_nested(objects))

    def summarise_dmr_count(self, kind='significant'):
        """
        Generate a dict summarising the number of DMRs in each result set
        :param kind: Specify the type of results to include. Valid options: significant, relevant, all
        :return:
        """
        if kind == 'significant':
            func = lambda x: len(x.results_significant)
        elif kind == 'relevant':
            func = lambda x: len(x.results_relevant)
        elif kind == 'all':
            func = lambda x: len(x.results)
        else:
            raise NotImplementedError("Unrecognised kind of result to include %s" % kind)
        return self.apply(func)


# TODO: work in progress!!
# class _DmrResultCollection(type):
#     """
#     Class factory that programmatically adds properties to the base _DmrResultsCollection class
#     This is useful so that we can have an arbitrary number of recursive methods that can be applied using `apply`.
#     """
#
#     def __init__(cls, name, bases, namespace):
#         super(_DmrResultCollection, cls).__init__(name, bases, namespace)
#         recursive_props = (
#             'clusters_relevant',
#             'clusters_significant',
#             'results_relevant',
#             'results_significant',
#         )
#
#         recursive_methods_by_cls = (
#             'clusters_relevant_by_class',
#             'clusters_significant_by_class',
#             'results_relevant_by_class',
#             'results_significant_by_class',
#         )
#
#         for rp in recursive_props:
#             setattr(cls, rp, property(lambda x: x.apply(rp)))
#
#         for rm in recursive_methods_by_cls:
#             def the_func(obj, cls=None):
#                 return obj.apply(lambda x: x.results_relevant_by_class(cls=cls))
#             setattr(cls, rm, the_func)


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


def compute_cross_dmr(
        me_data,
        me_meta,
        anno,
        pids,
        dmr_params,
        external_references=(('GIBCO', 'NSC'),),
):
    """
    Compute every possible DMR permutation in the available data.
    This function is quite specific to our setup, so might not be easily generalised (TODO)
    It computes every possible combination of GBM vs iNSC and GBM vs ref NSC.
    :return DmrResultCollection
    """

    obj = DmrResults(anno=anno)
    obj.identify_clusters(**dmr_params)
    res = {}

    # loop over GBM groups
    for pid1 in pids:
        res.setdefault(pid1, {})
        the_idx1 = me_meta.index.str.contains(pid1) & (me_meta.loc[:, 'type'] == 'GBM')
        # loop over iNSC groups
        for pid2 in pids:
            the_idx2 = me_meta.index.str.contains(pid2) & (me_meta.loc[:, 'type'] == 'iNSC')
            the_idx = the_idx1 | the_idx2
            the_groups = me_meta.loc[the_idx, 'type'].values
            the_samples = me_meta.index[the_idx].groupby(the_groups).values()
            the_obj = obj.copy()
            the_obj.test_clusters(me_data,
                                  samples=the_samples,
                                  n_jobs=dmr_params['n_jobs'],
                                  min_median_change=dmr_params['delta_m_min'],
                                  method=dmr_params['dmr_test_method'],
                                  alpha=dmr_params['alpha'],
                                  **dmr_params['test_kwargs']
                                  )
            res[pid1][pid2] = the_obj

        # loop over external reference NSC groups
        for er, er_type in external_references:
            the_idx2 = me_meta.index.str.contains(er) & (me_meta.loc[:, 'type'] == er_type)
            the_idx = the_idx1 | the_idx2
            the_groups = me_meta.loc[the_idx, 'type'].values
            the_samples = me_meta.index[the_idx].groupby(the_groups).values()

            the_obj = obj.copy()
            the_obj.test_clusters(me_data,
                                  samples=the_samples,
                                  n_jobs=dmr_params['n_jobs'],
                                  min_median_change=dmr_params['delta_m_min'],
                                  method=dmr_params['dmr_test_method'],
                                  **dmr_params['test_kwargs']
                                  )
            res[pid1][er] = the_obj

    return DmrResultCollection(**res)


def count_dmr_genes(res):
    """

    :param res: test_results (or significant/relevant only) from compute_dmr. This is the final nested dictionary,
    keyed by cluster ID and with value dictionary containing the kw 'genes'
    :return: The number of genes within this set of results
    """
    the_genes = set()
    for x in res.values():
        the_genes.update(x['genes'])
    return len(the_genes)


def venn_set_to_wide_dataframe(
    data,
    venn_set,
    set_labels,
    include_sets=None,
    full_data=None,
    cols_to_include=('median_delta', 'padj'),
    direction_col='median_delta'
):
    """
    Given the input DMR data and Venn sets, generate a wide format dataframe containing all the data, one column
    per patient and one row per gene.
    Optionally filter the sets to include only a subset.
    Optionally include non-significant results too.
    :param data: Dict containing DMR results, keyed by the entries of set_labels.
    These are produced using DmrResults.to_table(), i.e. they are pd.DataFrames
    :param venn_set:
    :param set_labels:
    :param include_sets:
    :param full_data: If supplied, this has the same format as `data`, but the lists are complete so that even non-
    significant results can be accessed.
    :param cols_to_include: Iterable of columns to include in the output
    :param direction_col: The name of the column in the input data to use for determing direction of change.
    :return:
    """
    if include_sets is not None:
        venn_set = dict([
            (k, v) for k, v in venn_set.items() if k in include_sets
        ])

    res = []
    for k in venn_set:
        ids = venn_set[k]
        # populate with individual patient results
        blocks = []
        consistency_check = []
        for i, t in enumerate(k):
            pid = set_labels[i]
            cols = [pid] + ["%s_%s" % (pid, lbl) for lbl in cols_to_include]
            this_datum = pd.DataFrame(
                index=ids,
                columns=cols
            )
            if t == '1':
                # this member is included here
                this_datum.loc[ids, pid] = 'Y'
                for c in cols_to_include:
                    this_datum.loc[ids, "%s_%s" % (pid, c)] = data[pid].loc[ids, c]
                # consistency check
                cc = data[pid].loc[ids, direction_col]
                cc.name = pid
                consistency_check.append(cc)
            else:
                this_datum.loc[ids, pid] = 'N'
                if full_data is not None:
                    for c in cols_to_include:
                        this_datum.loc[ids, "%s_%s" % (pid, c)] = full_data[pid].loc[ids, c]

            blocks.append(this_datum)

        core_block = pd.concat(blocks, axis=1)
        # assess consistency of direction
        consist = pd.Series(index=ids)

        if len(consistency_check) > 0:
            consistency_check = pd.concat(consistency_check, axis=1)
            # figure out what kind of consistency check is required based on data type
            if isinstance(consistency_check.values[0, 0], str):
                idx = consistency_check.apply(lambda col: col == consistency_check.iloc[:, 0]).all(axis=1)
            else:
                idx = consistency_check.apply(
                    lambda col: np.sign(col) == np.sign(consistency_check.iloc[:, 0])
                ).all(axis=1)
            consist.loc[idx] = 'Y'
            consist.loc[~idx] = 'N'

        core_block.insert(core_block.shape[1], 'consistent', consist)
        res.append(core_block)

    # check: no features should be in more than one data entry
    for i, k in enumerate(venn_set):
        for j, k2 in enumerate(venn_set):
            if k == k2: continue
            bb = len(res[i].index.intersection(res[j].index))
            if bb > 0:
                raise AttributeError("Identified %d genes that are in BOTH %s and %s" % (bb, k, k2))

    res = pd.concat(res, axis=0)

    return res



if __name__ == "__main__":
    OUTDIR = unique_output_dir('dmr', reuse_empty=True)

    anno = methylation_array.load_illumina_methylationepic_annotation()
    b, meta = methylation_array.hgic_methylationepic('swan')
    m = m_from_beta(b)

    gbm_sample_dict = {
        '018': 'GBM018_P10',
        '019': 'GBM019_P4',
        '031': 'GBM031_P4',
    }
    dura_sample_dict = {
        '018': 'DURA018_NSC_N2_P6',
        '019': 'DURA019_NSC_N8C_P2',
        '031': 'DURA031_NSC_N44B_P2',
    }

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
        samples = (gbm_sample_dict[sid], dura_sample_dict[sid])
        test_results[sid] = test_clusters_in_place(clusters, m, samples=samples, min_median_change=dm, n_jobs=n_jobs)
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
            samples = [dura_sample_dict[sid], gbm_sample_dict[sid]]

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
