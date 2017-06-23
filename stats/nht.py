from scipy import stats, asarray
import numpy as np
import warnings


class RFunctionDeferred(object):
    """
    Used for R functions that also require library imports.
    Rather than importing all these libraries when the Python module is imported, we defer it until the function is
    actually required
    """
    def __init__(self, func, imports=None):
        self.imports = imports or []
        self.func = func
        self.ready = False

    def import_packages(self):
        for im in self.imports:
            if not rpackages.isinstalled(im):
                utils = rpackages.importr('utils')
                utils.chooseCRANmirror(ind=1)
                utils.install_packages(im)
            rpackages.importr(im)
        self.ready = True

    def __call__(self, *args, **kwargs):
        if not self.ready:
            self.import_packages()
        return self.func(*args, **kwargs)


RFUNCTIONS_PRESENT = True
try:
    from rpy2 import robjects
    from rpy2.robjects import FloatVector, Formula
    import rpy2.robjects.packages as rpackages
    qwilcox = lambda q, n1, n2: robjects.r('qwilcox')(q, n1, n2)[0]
    pwilcox = lambda u, n1, n2: robjects.r('pwilcox')(u, n1, n2)[0]

    def wilcoxsign_test_func(x, y, distribution='exact'):
        frm = Formula("y ~ x")
        frm.environment["x"] = FloatVector(x)
        frm.environment["y"] = FloatVector(y)
        res = robjects.r("wilcoxsign_test")(frm)
        return robjects.r('pvalue')(res)[0]

    wilcoxsign_test = RFunctionDeferred(wilcoxsign_test_func, imports=['coin'])

except Exception:
    RFUNCTIONS_PRESENT = False
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



