from scipy import stats
import numpy as np


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
    N = n1 + n2

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
