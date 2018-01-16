import pandas as pd
from utils import rinterface
from rpy2 import robjects
from rpy2.robjects import pandas2ri, r
pandas2ri.activate()

toptags_cols = [
    "logFC", "unshrunk.logFC", "logCPM", "PValue", "FDR"
]

def _edger_func_glmqlfit(the_data, the_groups, the_contrast, fdr=0.01, lfc=1, return_full=False):
    rdata = pandas2ri.py2ri(the_data)
    rgroups = robjects.FactorVector(the_groups)
    y = r("DGEList")(rdata)
    y = r("calcNormFactors")(y)
    formula = robjects.Formula("~0 + groups")
    formula.environment['groups'] = rgroups
    design = r("model.matrix")(formula)
    design.colnames = r('levels')(rgroups)
    y = r("estimateDisp")(y, design)
    rcontrast = r('makeContrasts')(robjects.StrVector([the_contrast]), levels=design)
    fit = r('glmQLFit')(y, design)
    lrt = r('glmTreat')(fit, contrast=rcontrast, lfc=lfc)
    if return_full:
        toptags = r('topTags')(lrt, n=r('Inf'), **{'p.value': 1.})
    else:
        toptags = r('topTags')(lrt, n=r('Inf'), **{'p.value': fdr})
    if len(toptags) == 0:
        return pd.DataFrame(columns=toptags_cols)
    else:
        return pandas2ri.ri2py_dataframe(toptags[toptags.names.index('table')])


edger_glmqlfit = rinterface.RFunctionDeferred(_edger_func_glmqlfit, imports=['edgeR'])


def _edger_func_glmfit(the_data, the_groups, the_contrast, fdr=0.01, lfc=1, return_full=False):
    rdata = pandas2ri.py2ri(the_data)
    rgroups = robjects.FactorVector(the_groups)
    y = r("DGEList")(rdata)
    y = r("calcNormFactors")(y)
    formula = robjects.Formula("~0 + groups")
    formula.environment['groups'] = rgroups
    design = r("model.matrix")(formula)
    design.colnames = r('levels')(rgroups)
    y = r("estimateDisp")(y, design)
    rcontrast = r('makeContrasts')(robjects.StrVector([the_contrast]), levels=design)
    fit = r('glmFit')(y, design)
    lrt = r('glmTreat')(fit, contrast=rcontrast, lfc=lfc)
    if return_full:
        toptags = r('topTags')(lrt, n=r('Inf'), **{'p.value': 1.})
    else:
        toptags = r('topTags')(lrt, n=r('Inf'), **{'p.value': fdr})
    if len(toptags) == 0:
        return pd.DataFrame(columns=toptags_cols)
    else:
        return pandas2ri.ri2py_dataframe(toptags[toptags.names.index('table')])


edger_glmfit = rinterface.RFunctionDeferred(_edger_func_glmfit, imports=['edgeR'])


def _edger_func_exacttest(the_data, the_groups, fdr=0.01, lfc=1, pair=None, return_full=False):
    """
    Run edgeR DE analysis without fitting a GLM. Instead, we just compare two groups. Only a single factor is supported.
    :param the_data:
    :param the_groups:
    :param fdr:
    :param lfc:
    :param pair: An iterable of two group names. If None, compare the first two groups.
    :return:
    """
    if pair is None:
        lvl, fct = pd.factorize(the_groups)
        pair = fct[:2]
    rpair = robjects.StrVector(pair)
    rdata = pandas2ri.py2ri(the_data)
    rgroups = robjects.FactorVector(the_groups)
    y = r("DGEList")(rdata, group=rgroups)
    y = r("calcNormFactors")(y)
    y = r("estimateDisp")(y)
    et = r('exactTest')(y, rpair)
    if return_full:
        toptags = r('topTags')(et, n=r('Inf'), **{'p.value': 1.})
    else:
        toptags = r('topTags')(et, n=r('Inf'), **{'p.value': fdr})
    if len(toptags) == 0:
        return pd.DataFrame(columns=toptags_cols)
    else:
        tt = pandas2ri.ri2py_dataframe(toptags[toptags.names.index('table')])
        if lfc is not None:
            tt = tt.loc[tt.loc[:, 'logFC'].abs() >= lfc]
        return tt


edger_exacttest = rinterface.RFunctionDeferred(_edger_func_exacttest, imports=['edgeR'])