import pandas as pd
import collections


def binary_combinations(n, include_zero=False):
    """
    Generator of binary strings n digits long.
    :param n: Number of bits.
    :param include_zero: If True, return zero string ('0000...') otherwise (default) start at 1.
    :return:
    """
    str_fmt = "{0:0%db}" % n
    start = 0 if include_zero else 1
    for j in range(start, 2 ** n):
        yield str_fmt.format(j)


def binary_combinations_sum_gte(n, nmin, **kwargs):
    """
    Generator of binary strings n digits long.
    The strings generated have at least nmin ones.
    :param n:
    :param nmin:
    :param kwargs: Passed to binary_combinations
    :return:
    """
    for bn in binary_combinations(n, **kwargs):
        if len(bn.replace('0', '')) >= nmin:
            yield bn


def venn_from_arrays(*args, **kwargs):
    """
    The input arguments *args contain a number of arrays. Each one is a list or similar, e.g. a list of strings
     representing genes. The entries are analysed for intersections and used to create the relevant venn sets.
     Counts are also computed for convenience
    :param args:
    :param kwargs:
    :return:
    """
    n = len(args)

    all_items = reduce(lambda x, y: set(x).union(y), args, set())
    venn_counts = {}
    venn_sets = {}
    for bn in binary_combinations(n):
        # generate binary representation
        this_add = set(all_items)
        this_subtract = set()
        for k in range(n):
            if bn[k] == '0':
                this_subtract = this_subtract.union(args[k])
            else:
                this_add = this_add.intersection(args[k])
        this_intersection = this_add.difference(this_subtract)
        venn_counts[bn] = len(this_intersection)
        venn_sets[bn] = list(this_intersection)

    return venn_sets, venn_counts


def intersection_with_threshold(min_n=None, *args):
    """
    Each element in args is an iterable. We compute the intersection amongst all of the items in those iterables.
    If min_n is supplied, we relax the requirement that an item must be in every iterable. Instead, it must be in
    gte min_n of them.
    NB if min_n == 1, we effectively compute a union.
    :param args:
    """
    if min_n is None or min_n == len(args):
        return reduce(lambda x, y: set(x).intersection(y), args)

    # count each item
    ct = collections.Counter()
    for x in args:
        for k in x:
            ct[k] += 1

    # keep with a threshold
    return set([k for k in ct if ct[k] >= min_n])