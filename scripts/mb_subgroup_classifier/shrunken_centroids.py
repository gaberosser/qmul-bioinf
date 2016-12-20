import numpy as np
import pandas as pd
from load_data import microarray_data, rnaseq_data, allen_human_brain_atlas
from scripts.mb_subgroup_classifier.load import load_xz_rnaseq, load_xiaonan_microarray
from microarray import process
from plotting import heatmap, utils
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
    from settings import DATA_DIR
    import os
    import multiprocessing as mp
    from scripts.comparison_rnaseq_microarray import consts

    # define shrinkage delta values
    n_core = 14
    deltas = np.linspace(0., 12., 60)
    # FLAT_PRIOR = False
    FLAT_PRIOR = True
    k = 10  # XV

    # it's useful to maintain a list of known upregulated genes
    nano_genes = []
    for grp, arr in consts.NANOSTRING_GENES:
        nano_genes.extend(arr)
    nano_genes.remove('EGFL11')
    nano_genes.append('EYS')

    # load Ncott data (285 non-WNT MB samples)

    ncott, ncott_meta = microarray_data.load_annotated_microarray_gse37382(
        aggr_field='SYMBOL',
        aggr_method='max'
    )
    sort_idx = ncott_meta.subgroup.sort_values().index
    ncott_meta = ncott_meta.loc[sort_idx]
    ncott = process.yugene_transform(ncott.loc[:, sort_idx])

    # X = ncott.copy()
    # m = ncott_meta.copy()

    # load Allen (healthy cerebellum)

    # he, he_meta = allen_human_brain_atlas.cerebellum_microarray_reference_data(agg_field='gene_symbol', agg_method='max')

    # combine

    # common_genes = ncott.index.intersection(he.index)
    # res = pd.DataFrame(index=common_genes, columns=ncott.columns.union(he.columns))
    # res.loc[common_genes, he.columns] = he.loc[common_genes].values
    # res.loc[common_genes, ncott.columns] = ncott.loc[common_genes].values
    # res = res.astype(float)
    #
    # he_meta = pd.DataFrame(data=[[None, None, 'control']] * he.columns.size, index=he.columns, columns=ncott_meta.columns)
    # meta = pd.concat((ncott_meta, he_meta), axis=0).loc[res.columns]
    # X = res.copy()
    # m = meta.copy()
    # X = process.yugene_transform(X)

    # load Kool dataset
    kool, kool_meta = microarray_data.load_annotated_microarray_gse10327(
        aggr_field='SYMBOL',
        aggr_method='max',
    )
    sort_idx = kool_meta.subgroup.sort_values().index
    kool_meta = kool_meta.loc[sort_idx]
    kool = process.yugene_transform(kool.loc[:, sort_idx])
    kool_meta.loc[:, 'subgroup'] = (
        kool_meta.loc[:, 'subgroup'].str
            .replace('A', 'WNT')
            .replace('B', 'SHH')
            .replace('E', 'Group 3')
            .replace('C', 'Group 4')
            .replace('D', 'Group 4')
    )

    # X = kool.copy()
    # m = kool_meta.copy()

    # load Robinson dataset
    robi, robi_meta = microarray_data.load_annotated_microarray_gse37418(aggr_field='SYMBOL', aggr_method='max')
    robi_meta = robi_meta.loc[~robi_meta.subgroup.isin(['U', 'SHH OUTLIER'])]
    sort_idx = robi_meta.subgroup.sort_values().index
    robi_meta = robi_meta.loc[sort_idx]
    robi = process.yugene_transform(robi.loc[:, sort_idx])
    robi_meta.loc[:, 'subgroup'] = robi_meta.subgroup.str.replace('G3', 'Group 3').replace('G4', 'Group 4')

    X = robi.copy()
    m = robi_meta.copy()

    # if we lump groups C and D together, the classification becomes almost perfect (as you might expect):
    # m.loc[:, 'subgroup'] = m.subgroup.replace('Group 4', 'Group C/D')
    # m.loc[:, 'subgroup'] = m.subgroup.replace('Group 3', 'Group C/D')

    # for statistical purposes, split the data into training/validation and test sets in the ratio 3:1

    # shuffle columns first
    n = X.columns.size
    cols = X.columns[np.random.permutation(n)]
    X = X.loc[:, cols]
    m = m.loc[cols]

    # now remove the top 25% for testing
    cut_idx = int(n * 0.25)
    X_test = X.iloc[:, :cut_idx]
    meta_test = m.iloc[:cut_idx]
    X_train = X.iloc[:, cut_idx:]
    meta_train = m.iloc[cut_idx:]

    idx, labels = meta_train.subgroup.factorize()

    # split test-training, maintaining proportion of classes
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
        np.random.shuffle(subgrp_ind)
        # pix = np.random.permutation(nl)
        # perm_ind = subgrp_ind[pix]
        ch = np.array_split(subgrp_ind, indices_or_sections=k)
        # shuffle, otherwise the last testing slots get no samples from smaller groups
        np.random.shuffle(ch)
        for j in range(k):
            chunks[j].extend(ch[j])

    # Xv
    train_idx = []
    test_idx = []
    train_labels = []
    test_labels = []

    n_err_test = pd.DataFrame(index=range(k), columns=deltas)
    n_err_train = pd.DataFrame(index=range(k), columns=deltas)
    clas = pd.DataFrame(index=range(k), columns=deltas)

    p = mp.Pool(n_core)
    jobs = []


    # for j in range(k):
    #     te = chunks[j]
    #     te_ind = X_train.columns[te]
    #     tr = reduce(operator.add, [chunks[z] for z in set(range(k)).difference([j])])
    #     tr_ind = X_train.columns[tr]
    #     test_idx.append(te_ind)
    #     train_idx.append(tr_ind)
    #
    #     test_data = X_train.loc[:, te_ind]
    #     train_data = X_train.loc[:, tr_ind]
    #     test_labels = meta_train.subgroup.loc[te_ind]
    #     train_labels = meta_train.subgroup.loc[tr_ind]
    #
    #     pargs = (deltas, train_data, train_labels, test_data, test_labels)
    #     pkwargs = dict(flat_prior=FLAT_PRIOR)
    #     jobs.append(
    #         p.apply_async(run_validation, args=pargs, kwds=pkwargs)
    #     )
    #
    # p.close()
    # for j in range(k):
    #     etrain, etest, cl = jobs[j].get(1e6)
    #     n_err_test.loc[j] = etest
    #     n_err_train.loc[j] = etrain
    #     clas.loc[j] = cl
    #     logger.info("Complete xv run %d", j)

    n_bootstrap = 500
    sample_names = np.array(X_train.columns.values, copy=True)
    te_cutoff = int(np.round(X_train.columns.size / float(k)))
    for j in range(n_bootstrap):
        np.random.shuffle(sample_names)

        te_ind = sample_names[:te_cutoff]
        tr_ind = sample_names[te_cutoff:]

        print "Run %d. Test samples: %s" % (
            j + 1,
            ', '.join(te_ind)
        )

        te = X_train.loc[:, te_ind]
        tr = X_train.loc[:, tr_ind]

        test_idx.append(te_ind.copy())
        train_idx.append(tr_ind.copy())

        te_labels = meta_train.loc[te_ind, 'subgroup'].copy()
        tr_labels = meta_train.loc[tr_ind, 'subgroup'].copy()

        print repr(te_labels)

        train_labels.append(tr_labels.copy())
        test_labels.append(te_labels.copy())

        pargs = (deltas, tr, tr_labels, te, te_labels)
        pkwargs = dict(flat_prior=FLAT_PRIOR)
        jobs.append(
            p.apply_async(run_validation, args=pargs, kwds=pkwargs)
        )

    p.close()
    for j in range(n_bootstrap):
        etrain, etest, cl = jobs[j].get(1e6)
        n_err_test.loc[j] = etest
        n_err_train.loc[j] = etrain
        clas.loc[j] = cl
        logger.info("Complete xv run %d / %d", j + 1, n_bootstrap)

    # final run with all training data for assessment on test
    etrain, etest, _ = run_validation(deltas, X_train, meta_train.subgroup, X_test, meta_test.subgroup, flat_prior=FLAT_PRIOR)
    obj = NearestCentroidClassifier(data=X_train, labels=meta_train.subgroup, flat_prior=FLAT_PRIOR)

    n_train = float(X_train.columns.size)
    bootstrap_error_rate = n_err_test.sum(axis=0) * k / (float(n_bootstrap * n_train))
    bootstrap_error_rate_std = n_err_test.std(axis=0) * np.sqrt(k / (float(n_bootstrap * n_train)))
    # xv_error_rate = n_err_test.sum(axis=0) / float(X_train.columns.size)
    train_error_rate_mean = n_err_train.sum(axis=0) / float(sum([t.size for t in train_idx]))
    train_error_rate = etrain / float(X_train.columns.size)
    test_error_rate = etest / float(X_test.columns.size)

    n_remain = pd.Series(index=deltas)
    for d in deltas:
        obj.shrink_centroids(d)
        n_remain.loc[d] = obj.num_nonzero_features

    if True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        train_error_rate.plot(ax=ax, label="Training", c='b')
        test_error_rate.plot(ax=ax, label="Test", c='k')
        bootstrap_error_rate.plot(ax=ax, label="Xv", lw=3, c='r')
        ax.fill_between(
            deltas,
            bootstrap_error_rate - bootstrap_error_rate_std,
            bootstrap_error_rate + bootstrap_error_rate_std,
            color='r',
            alpha=0.3)
        # xv_error_rate.plot(ax=ax, label="Xv", lw=3)
        ax.legend(loc='upper left')
        ax.set_xlabel("Shrinkage parameter, $\Delta$")
        ax.set_ylabel("Error rate")

    # find the minimum testing error, using the simplest model to resolve ties
    # min_delta = xv_error_rate[xv_error_rate == xv_error_rate.min()].sort_index(ascending=False).index[0]
    min_delta = bootstrap_error_rate[bootstrap_error_rate == bootstrap_error_rate.min()].sort_index(ascending=False).index[0]
    obj = NearestCentroidClassifier(data=X, labels=m.subgroup, flat_prior=FLAT_PRIOR)
    obj.shrink_centroids(min_delta)

    cm_train = obj.confusion_matrix_training()
    cm_robi = obj.confusion_matrix(robi, robi_meta.subgroup)
    cm_kool = obj.confusion_matrix(kool, kool_meta.subgroup)
    cm_ncott = obj.confusion_matrix(ncott, ncott_meta.subgroup)

    # classify RNA-Seq
    X_htseq = load_xz_rnaseq(kind='htseq', yugene=True, gene_symbols=X.index)
    X_cuff = load_xz_rnaseq(kind='cuff', yugene=True, gene_symbols=X.index)

    rna_cuff_class = [(c, obj.classify(X_cuff.loc[:, c])) for c in X_cuff.columns]
    rna_cuff_class_probs = [(c, obj.class_probabilities(X_cuff.loc[:, c])) for c in X_cuff.columns]
    rna_htseq_class = [(c, obj.classify(X_htseq.loc[:, c])) for c in X_htseq.columns]
    rna_htseq_class_probs = [(c, obj.class_probabilities(X_htseq.loc[:, c])) for c in X_htseq.columns]

    # load Xiao-Nan data and try to classify
    xnan_sample_names = ('Pt1299', 'ICb1299-I', 'ICb1299-III', 'ICb1299-IV')
    X_xnan, xnan_meta = load_xiaonan_microarray(yugene=True, gene_symbols=X.index, sample_names=xnan_sample_names)

    xnan_class = [(c, obj.classify(X_xnan.loc[:, c])) for c in X_xnan.columns]
    xnan_class_probs = [(c, obj.class_probabilities(X_xnan.loc[:, c])) for c in X_xnan.columns]

    # visualise classifier
    s = obj.remaining_class_centroid_diff.sort_values(by=['Group 4', 'Group 3', 'SHH', 'WNT'], axis=0, ascending=False).index
    a = obj.remaining_class_centroid_diff.loc[s, ['Group 4', 'Group 3', 'SHH', 'WNT']]
    ax, cax = heatmap.single_heatmap(
        a,
        vmin=-.4,
        vmax=.4,
        xticklabels=False,
        fig_kwargs=dict(figsize=(10, 3))
    )
    ax.set_xlabel('Gene')
    plt.tight_layout(rect=[0, 0, 1.1, 1])  # tweaked to fill horizontal space
    fig = ax.get_figure()
    fig.savefig('robi_classifier_all_heatmap.png', dpi=200)
    fig.savefig('robi_classifier_all_heatmap.pdf', dpi=200)

    # compare with RNA-Seq (all genes)
    b = X_cuff.subtract(obj.global_centroid, axis=0).loc[s]
    c = pd.concat((a, b), axis=1).loc[:, list(a.columns) + list(b.columns)]
    ax, cax = heatmap.single_heatmap(
        c,
        vmin=-.4,
        vmax=.4,
        xticklabels=False,
        fig_kwargs=dict(figsize=(10, 4))
    )

    # reduce to SHH, 3, and 4
    to_keep = (obj.remaining_class_centroid_diff.loc[:, ['Group 4', 'Group 3', 'SHH']].abs() > 0.001).any(axis=1)
    cl = obj.remaining_class_centroid_diff.loc[to_keep]
    s = cl.sort_values(by=['Group 4', 'Group 3', 'SHH', 'WNT'], axis=0, ascending=False).index

    ax, cax = heatmap.single_heatmap(
        cl.loc[s, ['Group 4', 'Group 3', 'SHH']],
        vmin=-.4,
        vmax=.4,
        cbar=False,
        fig_kwargs=dict(figsize=(13, 3))
    )
    fig = ax.get_figure()
    quadmesh = ax.collections[-1]
    cax = fig.colorbar(quadmesh, pad=0.02)
    utils.axis_border(cax.ax, c='0.3', lw=1)
    ax.xaxis.label.set_visible(False)
    plt.tight_layout(pad=0.0, rect=[.01, 0, 1.1, .98])

    fig = ax.get_figure()
    fig.savefig('robi_classifier_shh34_heatmap.png', dpi=200)
    fig.savefig('robi_classifier_shh34_heatmap.pdf', dpi=200)
