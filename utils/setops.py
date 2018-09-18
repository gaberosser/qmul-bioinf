import pandas as pd
import numpy as np
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


def binary_combinations_sum_eq(n, neq, **kwargs):
    """
    Generator of binary strings n digits long.
    The strings generated have exactly neq ones.
    :param n:
    :param neq:
    :param kwargs: Passed to binary_combinations
    :return:
    """
    for bn in binary_combinations(n, **kwargs):
        if len(bn.replace('0', '')) == neq:
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


def intersection_with_threshold(*args, **kwargs):
    """
    Each element in args is an iterable. We compute the intersection amongst all of the items in those iterables.
    If min_n is supplied, we relax the requirement that an item must be in every iterable. Instead, it must be in
    gte min_n of them.
    NB if min_n == 1, we effectively compute a union.
    :param args:
    :param kwargs: If min_n is specified, use this
    """
    min_n = kwargs.get('min_n')
    if min_n is None or min_n == len(args):
        return reduce(lambda x, y: set(x).intersection(y), args)

    # count each item
    ct = collections.Counter()
    for x in args:
        for k in x:
            ct[k] += 1

    # keep with a threshold
    return set([k for k in ct if ct[k] >= min_n])


def pairwise_similarity(arr, method='intersection'):
    """
    Compute the pairwise similarity between every possible pairing in the supplied data
    The similarity between A and B is defined as the fraction of elements in A that are also in B. This is
    directional.
    :param arr: An iterable of iterables. Each element is a list representing items, e.g. genes.
    :param method: Either 'intersection' (use the intersection to compute the similarity) or 'union' (use the union)
    :return: An array of similarity scores, each of which is between 0 and 1.
    """
    n = len(arr)
    res = []
    for i in range(n):
        for j in range(i + 1, n):
            gl0 = arr[i]
            gl1 = arr[j]
            if method == 'intersection':
                a = len(set(gl0).intersection(gl1))
                res.append(a / float(len(gl0)))
                res.append(a / float(len(gl1)))
            elif method == 'union':
                a = float(len(set(gl0).union(gl1)))
                res.append(len(gl0) / a)
                res.append(len(gl1) / a)
            else:
                raise ValueError("Unrecognised method %s" % method)

    return res


def reduce_intersection(*args):
    intersecter = lambda x, y: set(x).intersection(y)
    return reduce(intersecter, args)


def reduce_union(*args):
    unioner = lambda x, y: set(x).union(y)
    return reduce(unioner, args)


def venn_set_to_wide_dataframe(
    data,
    venn_set,
    set_labels,
    include_sets=None,
    full_data=None,
    cols_to_include=None,
    static_cols_to_include=None,
    consistency_check_col=None,
    consistency_check_method=None,
    run_sanity_check=False
):
    """
    Given the input DMR data and Venn sets, generate a wide format dataframe containing all the data, one column
    per patient and one row per gene.
    Optionally filter the sets to include only a subset.
    Optionally include non-significant results too.
    :param data: Dict containing DMR results, keyed by the entries of set_labels.
    These are produced using DmrResults.to_table(), i.e. they are pd.DataFrames
    :param venn_set:
    :param set_labels:
    :param include_sets:
    :param full_data: If supplied, this has the same format as `data`, but the lists are complete so that even non-
    significant results can be accessed.
    :param cols_to_include: Iterable of columns to include in the output. These are the columns that differ between
    the different members.
    :param static_cols_to_include: Iterable of columns to incude in the output. These are the columns that are
    identical across all members. NB we will not check this, so if variable columns are included here then the output
    will probably be undesirable.
    :param consistency_check_col: The name of the column in the input data to use for determing consistency among
    members.
    :param consistency_check_method: Supported options ('sign', 'equal'). This is the method used to assess consistency.
    If None (default), use the data type to guess the best comparison.
    :param run_sanity_check: If True (default), apply sanity checks to the results before returning.
    :return:
    """
    if cols_to_include is None:
        cols_to_include = []

    if static_cols_to_include is None:
        static_cols_to_include = []

    if include_sets is not None:
        venn_set = dict([
            (k, v) for k, v in venn_set.items() if k in include_sets
        ])

    res = []

    for k in venn_set:
        ids = venn_set[k]

        # we only need to add static columns once per block, so keep track with this indicator
        add_static_from = None

        # populate with individual patient results
        blocks = []
        consistency_check = []
        for i, t in enumerate(k):
            pid = set_labels[i]
            cols = [pid] + ["%s_%s" % (pid, lbl) for lbl in cols_to_include]
            this_datum = pd.DataFrame(
                index=ids,
                columns=cols
            )
            if t == '1':
                # this member is included here
                this_datum.loc[ids, pid] = 'Y'
                for c in cols_to_include:
                    this_datum.loc[ids, "%s_%s" % (pid, c)] = data[pid].loc[ids, c]

                # consistency check
                if consistency_check_col is not None:
                    cc = data[pid].loc[ids, consistency_check_col]
                    cc.name = pid
                    consistency_check.append(cc)

                # add static columns (if required)
                add_static_from = pid
            else:
                this_datum.loc[ids, pid] = 'N'
                if full_data is not None:
                    for c in cols_to_include:
                        this_datum.loc[ids, "%s_%s" % (pid, c)] = full_data[pid].loc[ids, c]

            blocks.append(this_datum)

        if add_static_from is not None and len(static_cols_to_include) > 0:
            static_block = data[add_static_from].loc[ids, static_cols_to_include]
            blocks = [static_block] + blocks
        elif len(static_cols_to_include) > 0 and full_data is not None:
            # this must have been the null set, no pid has been found to add the static variables
            # in this case, we use the full data to fill in variables
            static_block = full_data.values()[0].loc[ids, static_cols_to_include]
            blocks = [static_block] + blocks

        core_block = pd.concat(blocks, axis=1)

        if consistency_check_col is not None:
            # assess consistency
            consist = pd.Series(index=ids)

            if len(consistency_check) > 0:
                consistency_check = pd.concat(consistency_check, axis=1)

                if consistency_check_method is None:
                    # figure out what kind of consistency check is required based on data type
                    if isinstance(consistency_check.values[0, 0], str):
                        consistency_check_method = 'equal'
                    else:
                        consistency_check_method = 'sign'

                if consistency_check_method == 'sign':
                    idx = consistency_check.apply(
                        lambda col: np.sign(col) == np.sign(consistency_check.iloc[:, 0])
                    ).all(axis=1)
                elif consistency_check_method == 'equal':
                    idx = consistency_check.apply(lambda col: col == consistency_check.iloc[:, 0]).all(axis=1)
                else:
                    raise NotImplementedError("Unsupported consistency check method %s." % consistency_check_method)

                consist.loc[idx] = 'Y'
                consist.loc[~idx] = 'N'

            core_block.insert(core_block.shape[1], 'consistent', consist)
        res.append(core_block)

    if run_sanity_check:
        # check: no features should be in more than one data entry
        for i, k in enumerate(venn_set):
            for j, k2 in enumerate(venn_set):
                if k == k2: continue
                bb = len(res[i].index.intersection(res[j].index))
                if bb > 0:
                    raise AttributeError("Identified %d features that are in BOTH %s and %s" % (bb, k, k2))

    res = pd.concat(res, axis=0)

    return res


## TODO: merge this with plotting.venn.upset_plot_with_groups (which looks more efficient?)
def full_partial_unique_other_sets_from_groups(
        set_labels,
        subgroup_dict
):
    subgroup_ind = dict([
        (k, pd.Index(set_labels).isin(v)) for k, v in subgroup_dict.items()
    ])
    subgroup_size = dict([
        (k, len(v)) for k, v in subgroup_dict.items()
    ])

    groups = subgroup_dict.keys()

    # colours and sets for UpsetR plot
    sets_full = {}
    sets_partial = {}
    sets_unique = []
    sets_other = []

    for k in binary_combinations(len(set_labels)):
        this_k = np.array([t for t in k]).astype(bool)
        # case 1: single member of the set
        if this_k.sum() == 1:
            sets_unique.append(k)
        # case 2: there is more than one member, so need to distinguish between full, partial and other
        else:
            n_in_subgroup = []
            full_match = False
            for grp in groups:
                grp_idx = subgroup_ind[grp]
                # case 2a: number of members matches the number in the subgroup
                # now we need to determine whether it's an actual match
                n_in_subgroup.append(this_k[grp_idx].sum())
                if this_k[grp_idx].sum() == this_k.sum():
                    # right number of members to be a full match - but is it?
                    if this_k[grp_idx].sum() == subgroup_size[grp]:
                        # yes - full match
                        sets_full.setdefault(grp, []).append(k)
                        full_match = True
                        # no need to continue going through the subgroups for this k
                        break
            if not full_match:
                n_in_subgroup = np.array(n_in_subgroup)
                ii = np.where(n_in_subgroup > 0)[0]
                # 2b: if only one element is non-zero, we have a partial match
                if ii.size == 1:
                    sets_partial.setdefault(groups[ii[0]], []).append(k)
                # 2c: must be in the other category
                else:
                    sets_other.append(k)

    return {
        'full': sets_full,
        'partial': sets_partial,
        'specific': sets_unique,
        'mixed': sets_other
    }


def specific_sets(set_labels):
    n = len(set_labels)
    res = {}
    for i, lbl in enumerate(set_labels):
        arr = np.array(['0'] * n)
        arr[i] = '1'
        res[lbl] = ''.join(arr)
    return res
