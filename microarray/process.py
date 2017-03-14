import uuid

import numpy as np
import pandas as pd

from utils.log import get_console_logger


def median_absolute_deviation(data, axis=1):
    """
    Compute the MAD on the supplied data
    :param data: A pandas DataFrame or an object that can be used to create one.
    :param axis: The axis along which we will compute the MAD
    :return: pandas Series with the same index as data containing the MAD
    """
    if axis == 0:
        axi = 0
        axj = 1
    elif axis == 1:
        axi = 1
        axj = 0
    else:
        ValueError("Axis must be 0 or 1")
    data = pd.DataFrame(data)
    return data.subtract(data.median(axis=axi), axis=axj).abs().median(axis=axi)


def aggregate_by_probe_set(marray_data, method='median', groupby='gene_symbol'):
    """
    Aggregate the microarray DataFrame using a pre-existing column.
    For example, we commonly want to go from probe set -> gene activity.
    The pre-existing column becomes the new index.
    :param marray_data:
    :param lookup: Optionally supply a lookup list that is used to
    :param method: Either a string specifying a common method (max, mean, sum, median) or a vectorised function
    :param groupby:
    :return:
    """
    # we need a temporary column name - pick a random one
    uniq_key = str(uuid.uuid4().get_hex())[:12]
    while uniq_key in marray_data.columns:
        uniq_key = str(uuid.uuid4().get_hex())[:12]

    grp_data = marray_data.groupby(groupby, axis=0)
    if method == 'max':
        data = grp_data.max()
    elif method == 'mean':
        data = grp_data.mean()
    elif method == 'sum':
        data = grp_data.sum()
    elif method == 'median':
        data = grp_data.median()
    elif method == 'max_std':
        # mask the groupby column
        dat = marray_data.drop(groupby, axis=1)
        s = dat.std(axis=1)
        s.name = uniq_key
        s = pd.concat((s, marray_data.loc[:, groupby]), axis=1)
        probes = s.groupby(groupby, axis=0).apply(lambda x: x.loc[:, uniq_key].idxmax())
        data = marray_data.loc[probes].set_index(groupby)
    elif method == 'max_mad':
        # maximum by median absolute deviation
        # mask the groupby column
        dat = marray_data.drop(groupby, axis=1)
        mad = median_absolute_deviation(dat, axis=1)
        mad.name = uniq_key
        mad = pd.concat((mad, marray_data.loc[:, groupby]), axis=1)
        probes = mad.groupby(groupby, axis=0).apply(lambda x: x.loc[:, uniq_key].idxmax())
        data = marray_data.loc[probes].set_index(groupby)
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
    Assume data are all positive
    """
    res = marray_data.copy()
    # add columnwise offset to ensure all positive values
    colmin = res.min(axis=0)
    neg_warn = False
    for i in np.where(colmin < 0)[0]:
        res.iloc[:, i] -= colmin[i]
        neg_warn = True
    if neg_warn:
        logger = get_console_logger(__name__)
        logger.warning("Data contained negative values. Columnwise shift applied to correct this.")

    for t in marray_data.columns:
        col = res.loc[:, t].sort_values(ascending=False)
        cs = col.cumsum()
        s = col.sum()
        # numerical error: the final value in cumsum() may not equal the sum
        if cs[-1] != s:
            cs[cs == cs[-1]] = s
        a = 1 - cs / s
        res.loc[a.index, t] = a

    # a numerical error in cumsum() may result in some small negative values. Zero these.
    res[res < 0] = 0.

    # colmin = res.min(axis=0)
    # colmin[colmin >= 0] = 0.
    # res = res.subtract(colmin, axis=1)

    return res


def variance_stabilizing_transform(marray_data):
    """
    Requires rpy2 and the `vsn` package.
    Use the vsn package in R to compute the variance stabilised transform of the supplied raw data.
    :param marray_data: Must contain only numeric data - no gene symbol columns or similar
    """
    from rpy2 import robjects
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    robjects.r("library('vsn')")
    rmat = pandas2ri.py2ri(marray_data)
    rmat = robjects.r['data.matrix'](rmat)
    v = robjects.r['vsn2'](rmat)
    v = robjects.r['predict'](v, newdata=rmat)
    dat = np.asarray(v)
    return pd.DataFrame(dat, index=marray_data.index, columns=marray_data.columns)
