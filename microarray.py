import pandas as pd


def aggregate_by_probe_set(marray_data, method='median', groupby='gene_symbol'):
    """
    Aggregate the microarray DataFrame using a pre-existing column (or index).
    For example, we commonly want to go from probe set -> gene activity.
    The pre-existing column becomes the new index.
    :param marray_data:
    :param lookup: Optionally supply a lookup list that is used to
    :param method: Either a string specifying a common method (max, mean, sum, median) or a vectorised function
    :param groupby:
    :return:
    """
    grp_data = marray_data.groupby(groupby)
    if method == 'max':
        data = grp_data.max()
    elif method == 'mean':
        data = grp_data.mean()
    elif method == 'sum':
        data = grp_data.sum()
    elif method == 'median':
        data = grp_data.median()
    else:
        # try calling the supplied method directly
        data = grp_data.agg(method)

    return data