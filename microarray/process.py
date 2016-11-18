import pandas as pd


def aggregate_by_probe_set(marray_data, method='median', groupby='gene_symbol', axis=0):
    """
    Aggregate the microarray DataFrame using a pre-existing column (or index).
    For example, we commonly want to go from probe set -> gene activity.
    The pre-existing column becomes the new index.
    :param marray_data:
    :param lookup: Optionally supply a lookup list that is used to
    :param method: Either a string specifying a common method (max, mean, sum, median) or a vectorised function
    :param groupby:
    :param axis: Axis to use when grouping. 0 denotes grouping by the values in a column.
    :return:
    """
    grp_data = marray_data.groupby(groupby, axis=axis)
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


def yugene_transform(marray_data):
    """
    Apply the YuGene transform to the supplied data.
    Le Cao, Kim-Anh, Florian Rohart, Leo McHugh, Othmar Korn, and Christine A. Wells.
    "YuGene: A Simple Approach to Scale Gene Expression Data Derived from Different Platforms for Integrated Analyses."
    Genomics 103, no. 4 (April 2014): 239-51. doi:10.1016/j.ygeno.2014.03.001.
    Assume the data are supplied with samples in columns and genes in rows
    """
    res = marray_data.copy()
    for t in marray_data.columns:
        col = res.loc[:, t].sort_values(ascending=False)
        a = 1 - col.cumsum() / col.sum()
        res.loc[a.index, t] = a
    return res