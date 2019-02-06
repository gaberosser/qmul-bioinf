import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import multiprocessing as mp

from utils import setops, output
from plotting import venn, common


"""
Aim: compare two different Venn diagrams, each showing the overlap of DE genes in similar contexts (one is iAPC,
the other iOPC). We want some kind of metric / statistic to show whether they seem to arise from a similar process.
"""

def marginal_totals(dat, membership):
    n_tot = {}
    for pid in dat:
        n_tot[pid] = {}
        for k in membership:
            n_tot[pid][k] = 0
            for x in membership[k]:
                n_tot[pid][k] += dat[pid][x]
    return n_tot


if __name__ == "__main__":
    n_iter = 1000
    apc_data = {
        '019': {'100': 309, '010': 838, '001': 756, '110': 672, '011': 1034, '101': 602, '111': 155}
    }
    opc_data = {
        '019': {'100': 298, '010': 337, '001': 197, '110': 203, '011': 42, '101': 267, '111': 1}
    }

    # approach 1: APC is the baseline, fix the marginal number for each comparison
    membership = {
        'A': ['100', '110', '101', '111'],
        'B': ['010', '110', '011', '111'],
        'C': ['001', '101', '011', '111'],
    }

    n_tot_apc = marginal_totals(apc_data, membership)
    n_tot_opc = marginal_totals(opc_data, membership)

    pid = '019'
    apc_tot = sum(apc_data[pid].values())
    opc_tot = sum(opc_data[pid].values())
    rvs_apc = dict([
        (
            k,
            [np.random.choice(range(apc_tot), replace=False, size=n_tot_apc[pid][k]) for i in range(n_iter)]
        ) for k in membership
    ])

    pool = mp.Pool()
    jobs = {}

    simulated_set_sizes = {}
    for i in range(n_iter):
        jobs[i] = pool.apply_async(
            setops.venn_from_arrays,
            args=tuple(rvs_apc[k][i] for k in membership)
        )

    pool.close()
    pool.join()

    for i, j in jobs.items():
        _, this_vc = j.get()
        for k, v in this_vc.items():
            simulated_set_sizes.setdefault(k, []).append(v)

    df = pd.DataFrame(simulated_set_sizes).transpose()

    # we can see by eye that the intersection is too large!
    # move on to null 2: weighted gene selection
    weights = np.array(
        reduce(lambda x, y: x + y, [
            [float(sum([int(t) for t in k]))] * v for k, v in apc_data[pid].items()
        ])
    )
    weights = weights / float(weights.sum())

    rvs_apc = dict([
        (
            k,
            [np.random.choice(range(apc_tot), replace=False, size=n_tot_apc[pid][k], p=weights) for i in range(n_iter)]
        ) for k in membership
    ])

    simulated_set_sizes = {}
    for i in range(n_iter):
        _, this_vc = setops.venn_from_arrays(*[rvs_apc[k][i] for k in membership])
        for k, v in this_vc.items():
            simulated_set_sizes.setdefault(k, []).append(v)
    df = pd.DataFrame(simulated_set_sizes).transpose()
