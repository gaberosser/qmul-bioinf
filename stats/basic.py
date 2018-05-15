import numpy as np


def construct_contingency(x, y):
    return np.array([
        [((x < 0) & (y < 0)).sum(), ((x > 0) & (y < 0)).sum()],
        [((x < 0) & (y > 0)).sum(), ((x > 0) & (y > 0)).sum()],
    ])


def ecdf_func(x):
    x = np.sort(x)
    n = float(x.size)
    def result(v):
        return np.searchsorted(x, v, side='right') / n
    return result