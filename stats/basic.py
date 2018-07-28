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


def ecdf(dat, xi):
    """
    Calculate the ecdf of every column in dat (if it is matrix-like)
    :param dat:
    :param xi: The values to use as a lookup
    :return: M x N matrix, where M is len(xi) and N is the number of columns in dat
    """
    dat = np.array(dat)

    b_flatten = False
    if len(dat.shape) == 1:
        dat = np.array([dat]).transpose()
        b_flatten = True

    n = dat.shape[1]
    m = len(xi)

    res = np.zeros((m, n))

    for i in range(n):
        func = ecdf_func(dat[:, i])
        res[:, i] = func(xi)

    if b_flatten:
        res = res[:, 0]

    return res