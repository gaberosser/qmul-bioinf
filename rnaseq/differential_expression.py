from utils import rinterface, output
from rnaseq.filter import filter_by_cpm
from rpy2 import robjects
from rpy2.robjects import pandas2ri, r
pandas2ri.activate()


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
    return pandas2ri.ri2py_dataframe(toptags[toptags.names.index('table')])

edger = rinterface.RFunctionDeferred(_edger_func, imports=['edgeR'])


if __name__ == '__main__':
    ## TODO: move this to a separate script
    from load_data import rnaseq_data
    pids = ['017', '050', '054', '061']
    obj = rnaseq_data.load_by_patient(pids, annotate_by='Ensembl Gene ID')

    de = {}
    de_gibco = {}

    for pid in pids:
        the_idx = obj.meta.index.str.contains(pid)
        the_data = obj.data.loc[:, the_idx]

        the_data = filter_by_cpm(the_data.loc[obj.data.index.str.contains('ENSG')], min_n_samples=1)
        the_genes = the_data.index

        the_groups = obj.meta.loc[the_idx, 'type'].values
        the_contrast = "GBM - iNSC"

        de[pid] = edger(the_data, the_groups, the_contrast)

        # repeat with gibco reference
        # use the same genes, rather than filtering again
        the_idx = (obj.meta.index.str.contains(pid) & (obj.meta.type == 'GBM')) | obj.meta.index.str.contains('GIBCO')
        the_data = obj.data.loc[the_genes, the_idx]
        the_groups = obj.meta.loc[the_idx, 'type'].values
        the_contrast = "GBM - NSC"

        de_gibco[pid] = edger(the_data, the_groups, the_contrast)


