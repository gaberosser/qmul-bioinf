import copy


def dict_iterator(d, parents=None, level=1, n_level=None):
    if parents is None:
        parents = []
    if (n_level is not None) and (level > n_level):
        yield (parents, d)
    else:
        for k in d.iterkeys():
            if isinstance(d[k], dict):
               for x in dict_iterator(d[k], parents + [k], level=level+1, n_level=n_level):
                   yield x
            else:
                yield (parents + [k], d[k])


def dict_by_sublevel(d, level, key, n_level=None):
    """
    Given the sublevel with specified key, rebuild the dictionary
    :param d:
    :param level:
    :param key:
    :param n_level:
    :return:
    """
    # translate level for zero indexing
    j = level - 1
    g = dict_iterator(d, n_level=n_level)
    res = {}
    for k, val in g:
        # if len(k) <= level:
        if len(k) < level:
            raise ValueError("Requested level is too low for this dictionary")
        if k[j] == key:
            remaining_keys = k[:j] + k[(j+1):]
            par = res
            for rk in remaining_keys[:-1]:
                par = par.setdefault(rk, {})
            par[remaining_keys[-1]] = val
    return res


def filter_dictionary(d, filt, n_level=None, copy=False):
    g = dict_iterator(d, n_level=n_level)
    res = {}
    for k, val in g:
        if not filt(val):
            continue
        par = res
        for rk in k[:-1]:
            par = par.setdefault(rk, {})
        par[k[-1]] = val
    if copy:
        res = copy.deepcopy(res)
    return res


def nested_dict_to_flat(x):
    """
    From the supplied nested dictionary, generate a flat dictionary in which the keys are tuples showing the
    original structure. Any object that is *not* a dictionary is considered an end node.
    """
    res = {}
    stack = [(None, None)] + x.items()
    key = []
    while True:
        k, x = stack.pop()
        if k is None and x is None:
            # signal to move up one in the nesting hierarchy
            if len(stack) == 0:
                # all done
                break
            else:
                # remove from the key
                key.pop()
        elif isinstance(x, dict):
            # leave a marker telling us where this nesting happened
            stack += [(None, None)]
            # append items to the stack
            stack.extend(x.items())
            # add to the key
            key.append(k)
        else:
            res[tuple(key + [k])] = x
    return res


def flat_dict_to_nested(x):
    res = {}
    for k, v in x.items():
        # create nesting structure
        parent = res
        for (i, t) in enumerate(k):
            if i == (len(k) - 1):
                parent[t] = v
            else:
                parent.setdefault(t, {})
                parent = parent[t]
    return res
