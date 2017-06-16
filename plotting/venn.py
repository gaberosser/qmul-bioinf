from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3


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
    str_fmt = "{0:0%db}" % n
    for j in range(1, 2 ** n):
        # generate binary representation
        this_add = set(all_items)
        this_subtract = set()
        bn = str_fmt.format(j)
        for k in range(n):
            if bn[k] == '0':
                this_subtract = this_subtract.union(args[k])
            else:
                this_add = this_add.intersection(args[k])
        this_intersection = this_add.difference(this_subtract)
        venn_counts[bn] = len(this_intersection)
        venn_sets[bn] = list(this_intersection)

    return venn_sets, venn_counts


def venn_diagram(*args, **kwargs):
    ax = kwargs.pop('ax', plt.gca())
    n = len(args)
    venn = None
    if n not in {2, 3}:
        raise NotImplementedError("At present, we only support 2 and 3 way Venn diagrams")
    venn_sets, venn_counts = venn_from_arrays(*args, **kwargs)
    if n == 2:
        venn = venn2(subsets=venn_counts, ax=ax, **kwargs)
    elif n == 3:
        venn = venn3(subsets=venn_counts, ax=ax, **kwargs)
    return venn, venn_sets, venn_counts
