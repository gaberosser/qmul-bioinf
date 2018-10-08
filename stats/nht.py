from scipy import stats, asarray
import numpy as np
import warnings
from utils import rinterface

RFUNCTIONS_PRESENT = rinterface.RFUNCTIONS_PRESENT

if RFUNCTIONS_PRESENT:
    from rpy2 import robjects
    from rpy2.robjects import FloatVector, Formula
    qwilcox = lambda q, n1, n2: robjects.r('qwilcox')(q, n1, n2)[0]
    pwilcox = lambda u, n1, n2: robjects.r('pwilcox')(u, n1, n2)[0]

    def wilcoxsign_test_func(x, y, distribution='exact'):
        frm = Formula("y ~ x")
        frm.environment["x"] = FloatVector(x)
        frm.environment["y"] = FloatVector(y)
        res = robjects.r("wilcoxsign_test")(frm)
        return robjects.r('pvalue')(res)[0]

    wilcoxsign_test = rinterface.RFunctionDeferred(wilcoxsign_test_func, imports=['coin'])

else:
    r_absent_exc = AttributeError("Requires R functions to be accessible using rpy2")


def mannwhitneyu_statistic(x, y):
    """
    Compute the Mann Whitney U statistic for the two samples.
    Extracted from scipy.stats.mannwhitneyu
    :return:
    """
    x = np.asarray(x)
    y = np.asarray(y)

    n1 = len(x)
    n2 = len(y)

    ranked = stats.rankdata(np.concatenate((x, y)))
    rankx = ranked[0:n1]  # get the x-ranks

    r1 = rankx.sum(axis=0)

    u1 = r1 - n1 * (n1 + 1) / 2.0
    u2 = n1 * n2 - u1

    return min(u2, u1)


def mannwhitneyu_test(x, y, alternative='two-tailed'):
    """
    Exact Mann Whitney U test.
    TODO: when number of samples is >20, use normal approximation
    :param x:
    :param y:
    :param alternative:
    :return:
    """
    if not RFUNCTIONS_PRESENT:
        raise r_absent_exc
    if alternative not in ('two-tailed',):
        # TODO: support greater, less
        raise AttributeError("Unsupported alternative")
    u = mannwhitneyu_statistic(x, y)
    p = pwilcox(u, len(x), len(y))
    if alternative == 'two-tailed':
        p = min(2.0 * p, 1.0)
    return p


def wilcoxon_signed_rank_statistic(x, y=None, zero_method="wilcox"):
    """
    Compute the T statistic required in the Wilcoxon Signed Rank Test for paired samples.
    This is the same calculation used in scipy.stats.wilcoxon
    See docs for that function for further explanation
    """

    if zero_method not in ["wilcox", "pratt", "zsplit"]:
        raise ValueError("Zero method should be either 'wilcox' "
                         "or 'pratt' or 'zsplit'")

    if y is None:
        d = asarray(x)
    else:
        x, y = map(asarray, (x, y))
        if len(x) != len(y):
            raise ValueError('Unequal N in wilcoxon.  Aborting.')
        d = x - y

    if zero_method == "wilcox":
        # Keep all non-zero differences
        d = np.compress(np.not_equal(d, 0), d, axis=-1)

    count = len(d)
    if count < 10:
        warnings.warn("Warning: sample size too small for normal approximation.")

    r = stats.rankdata(abs(d))
    r_plus = np.sum((d > 0) * r, axis=0)
    r_minus = np.sum((d < 0) * r, axis=0)

    if zero_method == "zsplit":
        r_zero = np.sum((d == 0) * r, axis=0)
        r_plus += r_zero / 2.
        r_minus += r_zero / 2.

    T = min(r_plus, r_minus)
    return T


def wilcoxon_signed_rank_test(x, y):
    """
    Compute the Wilcoxon signed rank test for paired samples. We use an R package to make the computation exact.
    :param x:
    :param y:
    :return:
    """
    if not RFUNCTIONS_PRESENT:
        raise r_absent_exc
    return wilcoxsign_test(x, y)


def spearman_exact(a, b, nperm=1000, seed=None):
    """
    Compute Spearman's r (nonparametric correlation coefficient) and the pvalue (testing the null: random ordering)
    Notably, we use a permutation approach here to compute the pvalue, rather than the t statistic used in Scipy's
    imlementation.
    As a result, this test is *very* slow compared to the approximate method (and, indeed, uses it in the for loop).
    :param a:
    :param b:
    :param nperm: The maximum number of permutations. If the sample size is sufficiently small, we use an exhaustive
    set of permutations instead.
    :return:
    """
    from scipy import special, stats
    import itertools
    this_r, _ = stats.spearmanr(a, b)

    if seed is None:
        rst = np.random.RandomState()
    else:
        rst = np.random.RandomState(seed=seed)

    if np.isnan(this_r):
        # failed, possibly due to completely tied observations
        return np.nan, np.nan

    if len(a) != len(b):
        raise ValueError("Length of arrays must match")

    r_perms = []
    if special.factorial(len(a)) < nperm:
        nperm = special.factorial(len(a))
        # exhaustive testing
        for t in itertools.permutations(a):
            rp, _ = stats.spearmanr(t, b)
            r_perms.append(rp)
    else:
        # permutation testing
        a_p = np.array(a, copy=True)
        for i in range(nperm):
            rst.shuffle(a_p)
            rp, _ = stats.spearmanr(a_p, b)
            r_perms.append(rp)

    # get pvalue from permutations
    r_perms = np.array(sorted(r_perms))
    # P value calculation: how many perms have an absolute value this large or larger?
    this_p = (np.abs(r_perms) >= np.abs(this_r)).sum() / float(nperm)
    return this_r, this_p

