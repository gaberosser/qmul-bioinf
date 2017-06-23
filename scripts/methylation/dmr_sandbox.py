import pickle
import pandas as pd
import numpy as np
from scipy import stats
from stats import nht
from methylation import dmr


def reshape_data(dat):
    if len(dat.shape) == 1:
        dat = dat.reshape((dat.size, 1))
    elif len(dat.shape) == 2:
        pass
    else:
        raise ValueError("Invalid data supplied with %d dimensions (must be either 1D or 2D)" % len(dat.shape))
    return dat


def wilcoxon_rank_sum_permutation(x, y, n_it=1999, return_stats=False):
    """
    :param x, y: k x n array, DataFrame or similar, where k is the number of probes and n is the number of patients
    """
    x = reshape_data(np.asarray(x))
    y = reshape_data(np.asarray(y))

    nx = x.shape[1]
    ny = y.shape[1]

    k = x.shape[0]
    # when we index into the flat data, the switchover point is the first of the y data points
    switch_idx = nx * k

    if k != y.shape[0]:
        raise ValueError("The number of probes in the two samples does not match.")

    flatdat = np.concatenate((x.flatten('F'), y.flatten('F')))

    stat = np.zeros(n_it)
    for i in range(n_it):
        # NB: many of these permutations are equivalent in the eyes of the MWU calculation: only permutations
        # *that switch data between the groups* have an effect. Still, it's as good a way as any to randomise?
        # In particular, I can't think of a better method when nx != ny
        this_dat = np.random.permutation(flatdat)
        u = nht.mannwhitneyu_statistic(this_dat[:switch_idx], this_dat[switch_idx:])
        stat[i] = u
    u = nht.mannwhitneyu_statistic(x.flatten(), y.flatten())
    p = (stat <= u).sum() / float(n_it)
    if return_stats:
        return (p, stat)
    return p


def wrs_exact_permutation(x, y, n_max=1999, return_stats=False):
    """
    :param n_max: The maximum number of iterations to run. If the exact test requires more than this, we revert to
    sampling. Setting this to None, 0 or negative forces exact sampling, but this might be very slow and expensive.

    This generates a strange result in some cases, e.g. the 1000th cluster of 018. The exact sampling approach generates
    a very different distribution of statistics from the approximate sampling approach.
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
        print "Sampling strategy (%d iterations)" % n_it
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
        print "Exact strategy (%d iterations)" % n_it
        stat = np.zeros(n_it)
        count = 0
        for i in range(n_it1):
            perm_x = x.copy()
            perm_y = y.copy()
            str_fmt = ("{0:0%db}" % k).format(i)
            idx1 = np.array(list(str_fmt)) == '1'
            for jx in range(nx):
                for jy in range(ny):
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

if __name__ == "__main__":
    with open('clusters.pkl', 'rb') as f:
        din = pickle.load(f)

    m = din['m']
    test_results = din['test_results']
    clusters = din['clusters']

    test_results_relevant = {}
    test_results_significant = {}

    for sid in test_results:
        test_results_relevant[sid] = dmr.filter_dictionary(
            test_results[sid],
            filt=lambda x: 'pval' in x,
            n_level=3,
            copy=False
        )
        test_results_significant[sid] = dmr.filter_dictionary(
            test_results[sid],
            filt=lambda x: x.get('rej_h0', False),
            n_level=3,
            copy=False
        )

    # get values associated with 018 repeats
    n_it = 1999
    samples = ('GBM018_P10', 'GBM018_P12', 'DURA018_NSC_N4_P4', 'DURA018_NSC_N2_P6')
    tmp = m.loc[:, samples]
    data018 = [tmp.loc[t[1]['probes']] for t in dmr.dict_iterator(clusters, n_level=3)]

    t = data018[0]
    xy = t.values.flatten(order='F')
    n_x = 2 * t.shape[0]
    perm_xy = [np.random.permutation(xy) for i in range(n_it)]
    u = nht.mannwhitneyu_statistic(xy[:n_x], xy[n_x:])
    perm_u = [nht.mannwhitneyu_statistic(p[:n_x], p[n_x:]) for p in perm_xy]
