import pandas as pd
from utils.log import get_console_logger
import numpy as np
from utils import rinterface
# from rpy2 import robjects
# from rpy2.robjects import pandas2ri, r
logger = get_console_logger(__name__)


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



def yugene_transform(marray_data, resolve_ties=True):
    """
    Apply the YuGene transform to the supplied data.
    Le Cao, Kim-Anh, Florian Rohart, Leo McHugh, Othmar Korn, and Christine A. Wells.
    "YuGene: A Simple Approach to Scale Gene Expression Data Derived from Different Platforms for Integrated Analyses."
    Genomics 103, no. 4 (April 2014): 239-51. doi:10.1016/j.ygeno.2014.03.001.
    Assume the data are supplied with samples in columns and genes in rows
    :param resolve_ties: If True (default), replace all tied values with the mean. This is especially significant at
    low count values, which are often highly degenerate.
    """

    res = marray_data.copy()
    # add columnwise offset to ensure all positive values
    colmin = res.min(axis=0)
    neg_warn = False
    for i in np.where(colmin < 0)[0]:
        res.iloc[:, i] -= colmin[i]
        neg_warn = True
    if neg_warn:
        logger.warning("Data contained negative values. Columnwise shift applied to correct this.")

    for t in marray_data.columns:
        col = res.loc[:, t].sort_values(ascending=False)
        cs = col.cumsum()
        s = col.sum()
        # numerical error: the final value in cumsum() may not equal the sum
        if cs[-1] != s:
            cs[cs == cs[-1]] = s
        a = 1 - cs / s

        if resolve_ties:
            # FIXME: this is tediously slow; can definitely improve it!
            # find tied values in the input data
            tied = np.unique(col.loc[col.duplicated()].values)
            if tied.size > 1:
                logger.info("Resolving %d ties in column %s.", tied.size - 1, t)
                for i in tied[tied > 0]:
                    a[col == i] = a[col == i].mean()
            else:
                logger.info("No ties to resolve in column %s.", t)

        res.loc[a.index, t] = a

    # a numerical error in cumsum() may result in some small negative values. Zero these.
    res[res < 0] = 0.

    # colmin = res.min(axis=0)
    # colmin[colmin >= 0] = 0.
    # res = res.subtract(colmin, axis=1)

    return res


def variance_stabilizing_transform(data):
    """
    Requires rpy2 and the `vsn` package.
    Use the vsn package in R to compute the variance stabilised transform of the supplied raw data.
    :param data: Must contain only numeric data - no gene symbol columns or similar
    """
    from rpy2 import robjects, rinterface
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    try:
        robjects.r("library('vsn')")
    except rinterface.RRuntimeError:
        logger.warning("Failed to import library vsn, attempting install")
        # attempt install
        rcmd = 'source("http://www.bioconductor.org/biocLite.R")'
        robjects.r(rcmd)
        robjects.r('biocLite("vsn")')

    rmat = pandas2ri.py2ri(data)
    rmat = robjects.r['data.matrix'](rmat)
    v = robjects.r['vsn2'](rmat, blind=True)
    v = robjects.r['predict'](v, newdata=rmat)
    dat = np.asarray(v)
    return pd.DataFrame(dat, index=data.index, columns=data.columns)


def _rlog_blind(the_data, remove_zero_counts=True):
    robjects = rinterface.robjects
    pandas2ri = rinterface.pandas2ri
    env = robjects.globalenv

    rdata = pandas2ri.py2ri(the_data)
    env['dat'] = rdata
    robjects.r("meta <- data.frame(rep('foo', ncol(dat)), row.names = colnames(dat))")
    robjects.r("dds <- DESeqDataSetFromMatrix(dat, meta, ~1)")
    if remove_zero_counts:
        robjects.r("dds <- dds[rowSums(counts(dds)) > 0,]")
    robjects.r("dds <- estimateSizeFactors(dds)")
    robjects.r("rld <- rlogTransformation(dds, blind = T)")

    rld = robjects.r("assay(rld)")

    return pd.DataFrame(rld, index=robjects.r("rownames(dds)"), columns=the_data.columns)


rlog_blind = rinterface.RFunctionDeferred(_rlog_blind, imports=['DESeq2'])


def _vst_blind(the_data, remove_zero_counts=True):
    robjects = rinterface.robjects
    pandas2ri = rinterface.robjects.pandas2ri
    env = robjects.globalenv

    rdata = pandas2ri.py2ri(the_data)
    env['dat'] = rdata
    robjects.r("meta <- data.frame(rep('foo', ncol(dat)), row.names = colnames(dat))")
    robjects.r("dds <- DESeqDataSetFromMatrix(dat, meta, ~1)")
    if remove_zero_counts:
        robjects.r("dds <- dds[rowSums(counts(dds)) > 0,]")
    robjects.r("dds <- estimateSizeFactors(dds)")
    robjects.r("rld <- vst(dds, blind = T)")

    rld = robjects.r("assay(rld)")

    return pd.DataFrame(rld, index=robjects.r("rownames(dds)"), columns=the_data.columns)


vst_blind = rinterface.RFunctionDeferred(_vst_blind, imports=['DESeq2'])
