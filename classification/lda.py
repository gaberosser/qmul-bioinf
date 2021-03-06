import numpy as np
import pandas as pd
import multiprocessing as mp
from utils import log
import operator

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
    # reinstate original sign
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
            # subtract mean class centroid from all samples in that class
            y.loc[:, self.labels == l] = y.loc[:, self.labels == l].subtract(self.class_centroids.loc[:, l], axis=0).values
        s_i_sq = (y ** 2).sum(axis=1) / (self.n - self.k)
        self.s_i = s_i_sq ** 0.5

        # s_0: median value of s over all genes
        self.s0 = self.s_i.median()

        # linear discriminant
        d_num = self.class_centroids.subtract(self.global_centroid, axis=0)

        # tiled dataframe (P x K) containing pooled w/c stdev
        a = tile_series(self.s_i, self.k, columns=self.labels_, index=self.data.index)

        self.s_eff = a + self.s0  # shape (P x K) but all K columns are identical
        # normalise so that each column gives the standard deviation of the matching column in d_num
        self.d_den = self.s_eff.multiply(self.m_k, axis=1)
        # this is the linear discriminant for each class
        # it's like a t-statistic in the sense that it is (sample mean - global mean) / sample std err
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
        # compute the new shrunken class centroids
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
        rel = self.class_centroids_shrink.subtract(self.global_centroid, axis=0)
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

    def confusion_matrix(self, mat, true_labels):
        """
        Generate the confusion matrix for the supplied data
        :param mat: P x M DataFrame, where M is the number of test samples
        :param true_labels: The true labels (length M).
        :return: pd.DataFrame
        """
        m = mat.columns.size
        label_idx, label_str = true_labels.factorize()
        # include all labels, even if some are absent from the supplied dataset
        label_str = label_str.union(self.labels_)

        # nl = len(self.labels_)
        nl = len(label_str)

        n_grp = pd.Series(0, index=label_str)
        cts = true_labels.groupby(true_labels).size()
        n_grp.loc[cts.index] = cts

        pred_labels = pd.Series([self.classify(mat.iloc[:, i]) for i in range(m)], index=mat.columns)

        cm = np.zeros((nl + 2, nl + 1))
        # for i, l in enumerate(self.labels_):
        for i, l in enumerate(label_str):
            # rows: inferred
            inf_idx = pred_labels[pred_labels == l].index
            if len(inf_idx) == 0:
                # skip this row - no inferred observations
                continue
            truth = true_labels[inf_idx]
            # for j, l2 in enumerate(self.labels_):
            for j, l2 in enumerate(label_str):
                # cols: truth
                cm[i, j] = (truth == l2).sum()
        # row sum
        cm[:nl, nl] = cm[:nl].sum(axis=1)
        # col sum
        cm[nl, :nl + 1] = cm[:nl, :nl + 1].sum(axis=0)
        # true positive rate
        a = cm.diagonal()[:nl]
        b = n_grp.astype(float).values
        c = np.zeros_like(a)
        c[b > 0] = a[b > 0] / b[b > 0]

        cm[nl + 1, :nl] = c

        cols = pd.MultiIndex.from_product([('True',), list(label_str) + ['Total',]])
        rows = pd.MultiIndex.from_product([('Inferred',), list(label_str) + ['Total', 'TPR']])
        cm = pd.DataFrame(data=cm, index=rows, columns=cols)
        return cm

    def confusion_matrix_training(self):
        return self.confusion_matrix(self.data, self.labels)


def run_validation(deltas, train_data, train_labels, test_data=None, test_labels=None, **kwargs):
    """
    :param test_data, test_labels: If provided, these are the data and labels for the validation set. Either both or
    neither must be provided. If None, only the training data are used for reporting.
    """
    b_test = True
    if test_data is None:
        if test_labels is not None:
            raise ValueError("Must provide either both or neither of test_data and test_labels")
        else:
            b_test = False
    elif test_labels is None:
        raise ValueError("Must provide either both or neither of test_data and test_labels")

    n_err_test = pd.Series(index=deltas)
    n_err_train = pd.Series(index=deltas)
    clas = pd.Series(index=deltas)
    obj = NearestCentroidClassifier(data=train_data, labels=train_labels, **kwargs)
    for d in deltas:
        obj.shrink_centroids(d)
        try:
            if b_test:
                c, nc, ni = obj.assess_test_data(test_data, true_labels=test_labels)
                n_err_test.loc[d] = ni
                clas.loc[d] = c.values
            _, nc, ni = obj.assess_training_data()
            n_err_train.loc[d] = ni
        except ValueError:
            # no more features
            continue
    if b_test:
        return n_err_train, n_err_test, clas
    else:
        return n_err_train


def cross_validation(deltas, data, labels, k=10, n_core=None, **kwargs):
    """

    :param deltas:
    :param data: Rows correspond to features, columns to samples
    :param labels:
    :param k: The -fold of the XV
    :param n_core: If not none, parallel processing is used
    :param kwargs:
    :return:
    """
    if n_core is not None:
        p = mp.Pool(n_core)
        jobs = []

    sample_names = np.random.permutation(data.columns)
    chunks = np.array_split(sample_names, indices_or_sections=k)
    true_clas = [labels.loc[ch] for ch in chunks]
    if any([len(t) == 0 for t in chunks]):
        raise ValueError("One or more groups contains no entries")
    if any([len(t) < 4 for t in chunks]):
        logger.warning("Small chunks (<4) detected. Consider reducing k?")

    train_idx = []
    test_idx = []

    n_err_test = pd.DataFrame(index=range(k), columns=deltas)
    n_err_train = pd.DataFrame(index=range(k), columns=deltas)
    clas = pd.DataFrame(index=range(k), columns=deltas)

    ind_set = set(range(k))

    for j in range(k):
        te_idx = chunks[j]
        te = data.loc[:, te_idx]
        tr_idx = np.concatenate([chunks[z] for z in ind_set.difference([j])])
        tr = data.loc[:, tr_idx]
        test_idx.append(te_idx)
        train_idx.append(tr_idx)

        te_labels = labels.loc[te_idx]
        tr_labels = labels.loc[tr_idx]

        pargs = (deltas, tr, tr_labels, te, te_labels)
        if n_core is not None:
            jobs.append(
                p.apply_async(run_validation, args=pargs, kwds=kwargs)
            )
        else:
            n_err_train.loc[j], n_err_test.loc[j], clas.loc[j] = run_validation(*pargs, **kwargs)
            logger.info("Complete xv run %d", j)

    if n_core is not None:
        p.close()
        for j in range(k):
            n_err_train.loc[j], n_err_test.loc[j], clas.loc[j] = jobs[j].get(1e6)
            logger.info("Complete xv run %d", j)

    return n_err_train, n_err_test, clas, true_clas