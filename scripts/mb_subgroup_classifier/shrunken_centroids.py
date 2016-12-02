import numpy as np
import pandas as pd
from load_data import microarray_data, allen_human_brain_atlas
import collections
import operator
import log
from matplotlib import pyplot as plt
import seaborn as sns

logger = log.get_console_logger(__name__)


def soft_threshold(x, delta):
    """
    Apply soft thresholding to the input array x
    :param x: Array of any dimension
    :param delta: The absolute quantity to subtract from each entry in x
    :return: A new array of the same shape as x
    """
    u = np.abs(x) - delta
    u[u < 0] = 0.
    u *= np.sign(x)
    return u


def tile_series(x, k, index=None, columns=None):
    """
    Tile the supplied series to create a DataFrame with k identical columns
    :param x: Input data
    :param k: Number of times to tile
    :return:
    """
    a = pd.DataFrame(
        np.tile(x, (k, 1)).transpose(),
        columns=columns,
        index=index
    )
    return a


class NearestCentroidClassifier(object):
    def __init__(self, data, labels, flat_prior=False):
        """
        Statistical model in which the data are described as a collection of centroids in p-dimensional space.
        :param data: The input data, array of dim P x N, where N is the number of samples and P is the number of features.
        :param labels: 1D array of length P containing the group labels
        :param flat_prior: If True, the prior probabilities for group membership are assumed to be uniform. Otherwise
         (the default) the group membership in the input data are assumed to reflect the underlying probability.
        """
        self.data = pd.DataFrame(data)
        self.p, self.n = self.data.shape
        self.labels = pd.Series(labels, index=self.data.columns)
        self.label_idx, self.labels_ = pd.Series(labels, index=self.data.columns).factorize()
        self.k = len(self.labels_)

        # compute global and group-specific centroids and standard errors
        self.global_centroid = self.data.mean(axis=1)
        # use .loc to ensure column order is maintained
        # not necessary, but tidy
        self.class_centroids = data.groupby(self.labels, axis=1).mean().loc[:, self.labels_]
        self.n_class = data.groupby(self.labels, axis=1).size().loc[self.labels_]
        self.m_k = np.sqrt(1 / self.n_class + 1 / float(self.n))

        if flat_prior:
            self.prior_prob = pd.Series(np.ones(self.k) / float(self.k), index=self.labels_)
        else:
            self.prior_prob = self.n_class / self.n_class.sum()

        # pooled within-class stdev
        y = self.data.copy()
        for i, l in enumerate(self.labels_):
            y.loc[:, self.labels == l] = y.loc[:, self.labels == l].subtract(self.class_centroids.loc[:, l], axis=0).values
        s_i_sq = (y ** 2).sum(axis=1) / (self.n - self.k)
        self.s_i = s_i_sq ** 0.5

        # s_0: median value of s over all genes
        self.s0 = self.s_i.median()

        # linear discriminant
        d_num = self.class_centroids.subtract(self.global_centroid, axis=0)

        # tiled dataframe (P x K) for convenience
        a = tile_series(self.s_i, self.k, columns=self.labels_, index=self.data.index)

        self.s_eff = a + self.s0
        self.d_den = self.s_eff.multiply(self.m_k)
        self.d_ik = d_num / self.d_den

        # placeholders for shrinkage operations
        self.d_ik_shrink = self.d_ik
        self.class_centroids_shrink = self.class_centroids
        self.delta = None

    def shrink_centroids(self, delta):
        """
        Shrink the class centroids by the absolute quantity delta using soft thresholding.
        This has the effect of setting self.class_centroids_shrink and self.d_ik_shrink
        :param delta: Absolute size of shrinkage, applied to linear discriminants
        :return: None
        """
        self.delta = delta
        self.d_ik_shrink = soft_threshold(self.d_ik, delta=delta)
        self.class_centroids_shrink = (self.d_den * self.d_ik_shrink).add(self.global_centroid, axis=0)
        if self.num_nonzero_features == 0:
            logger.warning("No features remain after shrinking. Attempting to classify will result in an error.")

    @property
    def num_nonzero_features(self):
        """
        :return: The number of features for which at least one of the subgroups has a non-zero linear discriminant.
        Features that are all zeros do not have any influence over the inferred classification.
        """
        return (self.d_ik_shrink != 0.).any(axis=1).sum()

    @property
    def remaining_d_ik(self):
        """
        :return: The entries of d_ik_shrink that have at least one non-zero component
        """
        return self.d_ik_shrink.loc[self.d_ik_shrink.any(axis=1)]

    @property
    def remaining_class_centroid_diff(self):
        """
        :return: The entries of class_centroid_shrink - global_centroid that have at least one non-zero component
        Sort by descending magnitude
        """
        rel = self.class_centroids_shrink.subtract(obj.global_centroid, axis=0)
        rel = rel.loc[rel.any(axis=1)]
        sort_idx = (rel ** 2).sum(axis=1).sort_values(ascending=False).index
        return rel.loc[sort_idx]

    def discriminant_score(self, vec):
        """
        Compute the discriminant score for an arbitrary vector of length P for each of the classes
        :param vec:
        :return: Array length K
        """
        if self.num_nonzero_features == 0:
            logger.error("No features remaining. Reduce shrinkage parameter.")
            raise ValueError("No features remaining. Reduce shrinkage parameter.")
        num = self.class_centroids_shrink.subtract(vec, axis=0)
        den = self.s_eff
        return ((num / den) ** 2).sum(axis=0) - 2 * np.log(self.prior_prob)

    def class_probabilities(self, vec):
        dk = self.discriminant_score(vec)
        # if the discriminant scores are high, all of the exponential functions could give a numerical 0.
        # therefore subtract the minimum value to guarantee at least one
        # this does not affect the RATIO of the values, which is what is measured.
        dk -= dk.min()
        gk = np.exp(-0.5 * dk)
        return gk / gk.sum()

    def classify(self, vec, int_index=False):
        """
        Compute the classification of the P-length array vec
        :param vec:
        :param int_index: If True, return the integer that indexes the label, else return the label string itself.
        :return:
        """
        dk = self.discriminant_score(vec)
        if int_index:
            dk.index = range(self.k)
        return dk.argmin()

    def assess_training_data(self):
        return self.assess_test_data(self.data, self.labels)

    def assess_test_data(self, mat, true_labels):
        """
        Run the classification on the samples in the input matrix
        :param mat: P x M DataFrame, where M is the number of test samples
        :param true_labels: The true labels (length M). These are used only after classification to assess the result.
        :return: inferred_labels, n_correct, n_incorrect
        """
        m = mat.columns.size
        res = pd.Series([self.classify(mat.iloc[:, i]) for i in range(m)], index=mat.columns)
        n_correct = (res == true_labels).sum()
        n_incorrect = m - n_correct
        return res, n_correct, n_incorrect


def run_validation(deltas, train_data, train_labels, test_data, test_labels, **kwargs):
    n_err_test = pd.Series(index=deltas)
    n_err_train = pd.Series(index=deltas)
    clas = pd.Series(index=deltas)
    obj = NearestCentroidClassifier(data=train_data, labels=train_labels, **kwargs)
    for d in deltas:
        obj.shrink_centroids(d)
        try:
            c, nc, ni = obj.assess_test_data(test_data, true_labels=test_labels)
            n_err_test.loc[d] = ni
            clas.loc[d] = c.values
            _, nc, ni = obj.assess_training_data()
            n_err_train.loc[d] = ni
        except ValueError:
            # no more features
            continue
    return n_err_train, n_err_test, clas


if __name__ == '__main__':
    # res, meta = microarray_data.load_annotated_microarray_gse37382(index_field='gene_symbol')
    from settings import DATA_DIR
    import os
    import multiprocessing as mp
    from scripts.comparison_rnaseq_microarray import consts

    # it's useful to maintain a list of known upregulated genes
    nano_genes = []
    for _, arr in consts.NANOSTRING_GENES:
        nano_genes.extend(arr)
    nano_genes.remove('EGFL11')
    nano_genes.append('EYS')

    meta_fn = os.path.join(DATA_DIR, 'microarray_GSE37382', 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=0, sep=',')
    infile = os.path.join(DATA_DIR, 'microarray_GSE37382', 'data.ann.txt.gz')
    res = pd.read_csv(infile, sep='\t', header=0, index_col=0)
    # fix sample names
    res.columns = res.columns.str.replace('.', '_')
    # aggregate by gene symbol
    res = res.drop(['ENTREZID', 'GENENAME', 'ACCNUM', 'ENSEMBL'], axis=1)
    res = res.loc[~res.SYMBOL.isnull()]
    res = res.groupby('SYMBOL', axis=0).median()

    # for statistical purposes, split the data into training/validation and test sets in the ratio 3:1

    # shuffle columns first
    n = res.columns.size
    cols = res.columns[np.random.permutation(n)]
    res = res.loc[:, cols]
    meta = meta.loc[cols]

    # now remove the top 25% for testing
    cut_idx = int(n * 0.25)
    X_test = res.iloc[:, :cut_idx]
    meta_test = meta.iloc[:cut_idx]
    X_train = res.iloc[:, cut_idx:]
    meta_train = meta.iloc[cut_idx:]

    idx, labels = meta_train.subgroup.factorize()

    # define shrinkage delta values
    deltas = np.linspace(0., 12., 60)
    FLAT_PRIOR = False

    # split test-training, maintaining proportion of classes
    k = 10
    training = collections.defaultdict(list)
    test = collections.defaultdict(list)

    # pre-compute the permuted indices for each class
    train_idx = []
    test_idx = []
    perm_ind = {}
    chunks = collections.defaultdict(list)
    for i, l in enumerate(labels):
        subgrp_ind = np.where(idx == i)[0]
        nl = (idx == i).sum()
        # random order
        pix = np.random.permutation(nl)
        perm_ind = subgrp_ind[pix]
        ch = np.array_split(perm_ind, indices_or_sections=k)
        for j in range(k):
            chunks[j].extend(ch[j])

    # Xv
    train_idx = []
    test_idx = []
    n_err_test = pd.DataFrame(index=range(k), columns=deltas)
    n_err_train = pd.DataFrame(index=range(k), columns=deltas)
    clas = pd.DataFrame(index=range(k), columns=deltas)

    p = mp.Pool()
    jobs = []


    for j in range(k):
        te = chunks[j]
        te_ind = X_train.columns[te]
        tr = reduce(operator.add, [chunks[z] for z in set(range(k)).difference([j])])
        tr_ind = X_train.columns[tr]
        test_idx.append(te_ind)
        train_idx.append(tr_ind)

        test_data = X_train.loc[:, te_ind]
        train_data = X_train.loc[:, tr_ind]
        test_labels = meta_train.subgroup.loc[te_ind]
        train_labels = meta_train.subgroup.loc[tr_ind]

        pargs = (deltas, train_data, train_labels, test_data, test_labels)
        pkwargs = dict(flat_prior=FLAT_PRIOR)
        jobs.append(
            p.apply_async(run_validation, args=pargs, kwds=pkwargs)
        )

    p.close()
    for j in range(k):
        etrain, etest, cl = jobs[j].get(1e6)
        n_err_test.loc[j] = etest
        n_err_train.loc[j] = etrain
        clas.loc[j] = cl
        logger.info("Complete xv run %d", j)

    # final run with all training data for assessment on test
    etrain, etest, _ = run_validation(deltas, X_train, meta_train.subgroup, X_test, meta_test.subgroup, flat_prior=FLAT_PRIOR)
    obj = NearestCentroidClassifier(data=X_train, labels=meta_train.subgroup, flat_prior=FLAT_PRIOR)

    xv_error_rate = n_err_test.sum(axis=0) / float(X_train.columns.size)
    train_error_rate_mean = n_err_train.sum(axis=0) / float(sum([t.size for t in train_idx]))
    train_error_rate = etrain / float(X_train.columns.size)
    test_error_rate = etest / float(X_test.columns.size)

    n_remain = pd.Series(index=deltas)
    for d in deltas:
        obj.shrink_centroids(d)
        n_remain.loc[d] = obj.num_nonzero_features

    fig = plt.figure()
    ax = fig.add_subplot(111)
    xv_error_rate.plot(ax=ax, label='XV error')
    train_error_rate_mean.plot(ax=ax, label='XV mean training error')
    train_error_rate.plot(ax=ax, label='Training error overall')
    test_error_rate.plot(ax=ax, label='Test error')
    ax.legend(loc='upper left')
    ax.set_xlabel('Shrinkage parameter')
    ax.set_ylabel('Error rate')

