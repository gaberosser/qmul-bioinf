import pandas as pd

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
     representing genes. The entries are analysed for intersections and used to create a Venn plot.
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


# def venn_from_array_with_constraints(constraint_fns=None, attributes=None, combine_op='and', *args):
#     """
#     In addition to finding matching elements, also check one or more arbitrary constraints based on attributes.
#     :param constraint_fns: Iterable of functions, each one to be applied to the corresponding element in attributes
#     :param attributes: Iterable of array of attributes, supplied as a pd.Series indexed by the IDs in args.
#     :param combine_op: Either 'and' or 'or'
#     :param args: Input arguments. Each element contains an iterable of IDs. As passed to `venn_from_arrays`.
#     :return:
#     """
#     n = len(args)
#
#     if combine_op not in {'and', 'or'}:
#         raise AttributeError("combine_op unrecognised. Options are 'and', 'or'")
#     if attributes is None:
#         raise AttributeError("Must supply attributes")
#     if constraint_fns is None:
#         raise AttributeError("Must supply constraint functions")
#     if len(constraint_fns) != len(attributes):
#         raise AttributeError("Must supply matching numbers of constraint functions and attributes")
#
#     # check that all attributes are present.
#     all_items = reduce(lambda x, y: set(x).union(y), args, set())
#     for attr in attributes:
#         if any(~pd.Index(all_items).isin(attr.index)):
#             raise KeyError(
#                 "One or more of the IDs is not matched in the attributes: %s" %
#                 ', '.join([str(t) for t in all_items.difference(attr.index)])
#             )
#
#     # run the basic venn algorithm
#     venn_sets, venn_counts = venn_from_arrays(*args)
#
#     # now run back through and check matching in each case
#     for attr, fun in zip(attributes, constraint_fns):
#
#         # dat = attr.loc[ID]
#         pass