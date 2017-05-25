from scipy import stats
import numpy as np


RFUNCTIONS_PRESENT = True
try:
    from rpy2 import robjects
    qwilcox = lambda q, n1, n2: robjects.r('qwilcox')(q, n1, n2)[0]
    pwilcox = lambda u, n1, n2: robjects.r('pwilcox')(u, n1, n2)[0]
    global wilcoxonsign_test
    wilcoxonsign_test = None
    ## FIXME: this doesn't work?

    def define_wilcoxonsrt():
        global wilcoxonsign_test
        robjects.r("library('coin')")
        wilcoxonsign_test = lambda x, y: robjects.r('wilcoxonsign_test')(x, y, distribution='exact')

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
    if wilcoxonsign_test is None:
        define_wilcoxonsrt()
    return wilcoxonsign_test(x, y)
