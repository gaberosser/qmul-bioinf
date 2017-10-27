import pandas as pd
from utils import rinterface
from rpy2 import robjects
from rpy2.robjects import pandas2ri, r
pandas2ri.activate()

toptags_cols = [
    "logFC", "unshrunk.logFC", "logCPM", "PValue", "FDR"
]

def _edger_func(the_data, the_groups, the_contrast, fdr=0.01, lfc=1):
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
    toptags = r('topTags')(lrt, n=r('Inf'), **{'p.value': fdr})
    if len(toptags) == 0:
        return pd.DataFrame(columns=toptags_cols)
    else:
        return pandas2ri.ri2py_dataframe(toptags[toptags.names.index('table')])

edger = rinterface.RFunctionDeferred(_edger_func, imports=['edgeR'])

