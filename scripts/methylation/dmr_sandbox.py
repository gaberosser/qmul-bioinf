import pickle
import pandas as pd
import numpy as np
from scipy import stats
from stats import nht
from methylation import dmr
from matplotlib import pyplot as plt
import multiprocessing as mp


def reshape_data(dat):
    if len(dat.shape) == 1:
        dat = dat.reshape((dat.size, 1))
    elif len(dat.shape) == 2:
        pass
    else:
        raise ValueError("Invalid data supplied with %d dimensions (must be either 1D or 2D)" % len(dat.shape))
    return dat


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

    # get new results
    s = (('GBM018_P10', 'GBM018_P12'), ('DURA018_NSC_N4_P4', 'DURA018_NSC_N2_P6'))
    new_results = dmr.test_clusters(clusters, m, s, n_jobs=4)

    # get values associated with 018 repeats
    n_it = 1999
    samples = ('GBM018_P10', 'GBM018_P12', 'DURA018_NSC_N4_P4', 'DURA018_NSC_N2_P6')
    tmp = m.loc[:, samples]
    data018 = [tmp.loc[t[1]['probes']] for t in dmr.dict_iterator(clusters, n_level=3)]

    # permutation Mann-Whitney U
    pool = mp.Pool(4)
    jobs = []
    cl018_pval_mwu = []
    for t in data018:
        x = t.iloc[:, :2]
        y = t.iloc[:, 2:]
        jobs.append(pool.apply_async(dmr.wilcoxon_rank_sum_permutation, args=(x, y)))
        cl018_pval_mwu.append(nht.mannwhitneyu_test(x.values.flat, y.values.flat))
    pool.close()
    cl018_pval_mwu_perm = [j.get(1e6) for j in jobs]

    #
    # p, stat = wilcoxon_rank_sum_permutation(x, y, n_max=1999, return_stats=True)
    # p_exact, stat_exact = wilcoxon_rank_sum_permutation(x, y, n_max=None, return_stats=True)
    # u_obs = nht.mannwhitneyu_statistic(x.values.flat, y.values.flat)
    # res_obs = stats.mannwhitneyu(x.values.flat, y.values.flat)
    # u_obs = res_obs.statistic
    # p_obs = res_obs.pvalue
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.hist(stat, 40, normed=True, label="Sampling strategy")
    # ax.hist(stat_exact, 40, normed=True, alpha=0.5, label="Exhaustive strategy")
    # ax.legend(loc='upper left')

    # investigate discrepancy between distributions
    t = data018[1000]
    x = np.asarray(t.iloc[:, :2])
    y = np.asarray(t.iloc[:, 2:])
    k = x.shape[0]
    nx = x.shape[1]
    ny = y.shape[1]
    n_it1 = 2 ** k

    multipliers = stats.binom.rvs(1, 0.5, size=(k, n_it))
    n_moved = multipliers.sum(axis=0)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(n_moved, np.arange(-1, k + 1) + .5, normed=True, label="Sampling strategy")
    # the distribution is itself binomial(k, 0.5)
    theor_pmf = stats.binom.pmf(range(k + 1), k, 0.5)
    ax.plot(np.arange(k + 1), theor_pmf, 'k-', label="Expected (binomial)")

    n_moved_e = []
    for i in range(n_it1):
        str_fmt = ("{0:0%db}" % k).format(i)
        n_moved_e.append((np.array(list(str_fmt)) == '1').sum())
    # this has the SAME distribution
    ax.hist(n_moved_e, np.arange(-1, k + 1) + .5, normed=True, alpha=0.5, label="Complete strategy")
    ax.legend(loc='upper left')

    # so, why are the U distributions so very different??
    # think it's the jx and jy component
    # try an approximate complete sampling (all idx covered, but random choice of jx, jy)
    stat2 = np.zeros(n_it1)
    high_scorers2 = []
    all_params2 = []
    count = 0
    for i in range(n_it1):
        perm_x = x.copy()
        perm_y = y.copy()
        str_fmt = ("{0:0%db}" % k).format(i)
        idx1 = np.array(list(str_fmt)) == '1'
        jx = np.random.randint(nx)
        jy = np.random.randint(ny)
        all_params2.append((i, jx, jy))
        # perform the data swap
        perm_x[idx1, jx] = y[idx1, jy]
        perm_y[idx1, jy] = x[idx1, jx]
        this_u = nht.mannwhitneyu_statistic(perm_x.flat, perm_y.flat)
        if this_u > 280:
            high_scorers2.append((i, jx, jy))
        stat2[count] = this_u
        count += 1

    # repeat the exhaustive (and exhausting!!) process
    stat3 = np.zeros(n_it1 * nx * ny)
    high_scorers3 = []
    all_params3 = []
    count = 0
    for i in range(n_it1):
        str_fmt = ("{0:0%db}" % k).format(i)
        idx1 = np.array(list(str_fmt)) == '1'
        for jx in range(nx):
            for jy in range(ny):
                perm_x = x.copy()
                perm_y = y.copy()
                all_params3.append((i, jx, jy))
                # perform the data swap
                perm_x[idx1, jx] = y[idx1, jy]
                perm_y[idx1, jy] = x[idx1, jx]
                this_u = nht.mannwhitneyu_statistic(perm_x.flat, perm_y.flat)
                if this_u > 280:
                    high_scorers3.append((i, jx, jy))
                stat3[count] = this_u
                count += 1

    stat4 = np.zeros(n_it1)
    count = 0
    for i in range(n_it1):
        perm_x = x.copy()
        perm_y = y.copy()
        str_fmt = ("{0:0%db}" % k).format(i)
        idx1 = np.array(list(str_fmt)) == '1'
        # perform the data swap
        perm_x[idx1] = y[idx1]
        perm_y[idx1] = x[idx1]
        stat4[count] = nht.mannwhitneyu_statistic(perm_x.flat, perm_y.flat)
        count += 1

