import pandas as pd
import general
import multiprocessing as mp
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


def run_one_de(the_data, the_groups, the_comparison, lfc=1, fdr=0.01, method='QLGLM', return_full=False, tax_id=9606):
    the_contrast = "%s - %s" % (the_comparison[0], the_comparison[1])
    if method == 'QLGLM':
        res = edger_glmqlfit(
            the_data,
            the_groups,
            the_contrast,
            lfc=lfc,
            fdr=fdr,
            return_full=return_full
        )
    elif method == 'GLM':
        res = edger_glmfit(
            the_data,
            the_groups,
            the_contrast,
            lfc=lfc,
            fdr=fdr,
            return_full=return_full
        )
    elif method == 'exact':
        res = edger_exacttest(
            the_data,
            the_groups,
            pair=the_comparison[::-1],
            lfc=lfc,
            fdr=fdr,
            return_full=return_full
        )
    else:
        raise AttributeError("Unrecognised method %s." % method)

    general.add_gene_symbols_to_ensembl_data(res, tax_id=tax_id)
    general.add_fc_direction(res)

    return res


def compute_cross_de(
        rnaseq_obj,
        pids,
        external_references=(('GIBCO', 'NSC'),),
        lfc=1,
        fdr=0.01,
        method='QLGLM',
        njob=None,
        tax_id=9606,
        type_field='type',
        default_types=('GBM', 'iNSC'),
):
    """
    Compute DE between every patient GBM sample and every _other_ healthy patient sample, in addition to paired DE.
    We can also include one or more external references (e.g. Gibco, the default).
    :param rnaseq_obj:
    :param pids:
    :param external_references: List of 2-tuples ('sample name contains', 'cell type')
    :param lfc:
    :param fdr:
    :param method:
    :param njob: Number of multiprocessing workers to use (default: None - automatically selected)
    :param tax_id: Taxonomy ID, used for labelling genes (default: human)
    :param type_field: The name of the column used to lookup the cell type
    :param default_types: The types to be used for the basic GBM vs iNSC comparison. The first entry is then also used
    in the external reference comparison.
    :return:
    """
    if method not in {'QLGLM', 'GLM', 'exact'}:
        raise NotImplementedError("Unsupported method.")
    de = {}

    pool = mp.Pool()
    jobs = {}

    for pid in pids:

        # cross comparison
        for pid2 in pids:
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, type_field] == default_types[0])) | \
                      (rnaseq_obj.meta.index.str.contains(pid2) & (rnaseq_obj.meta.loc[:, type_field] == default_types[1]))
            the_data = rnaseq_obj.data.loc[:, the_idx]
            the_groups = rnaseq_obj.meta.loc[the_idx, type_field].values
            the_comparison = list(default_types)
            if njob == 1:
                de[(pid, pid2)] = run_one_de(
                    the_data,
                    the_groups,
                    the_comparison,
                    lfc=lfc,
                    fdr=fdr,
                    method=method,
                    tax_id=tax_id
                )
            else:
                jobs[(pid, pid2)] = pool.apply_async(
                    run_one_de,
                    args=(the_data, the_groups, the_comparison),
                    kwds=dict(lfc=lfc, fdr=fdr, method=method, tax_id=tax_id)
                )

        # external reference comparison
        for er, er_type in external_references:
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, type_field] == default_types[0])) | \
                      (rnaseq_obj.meta.index.str.contains(er) & (rnaseq_obj.meta.loc[:, type_field] == er_type))
            the_data = rnaseq_obj.data.loc[:, the_idx]
            the_groups = rnaseq_obj.meta.loc[the_idx, type_field].values
            the_comparison = [default_types[0], er_type]
            if njob == 1:
                de[(pid, er)] = run_one_de(
                    the_data,
                    the_groups,
                    the_comparison,
                    lfc=lfc,
                    fdr=fdr,
                    method=method,
                    tax_id=tax_id
                )
            else:
                jobs[(pid, er)] = pool.apply_async(
                    run_one_de,
                    args=(the_data, the_groups, the_comparison),
                    kwds=dict(lfc=lfc, fdr=fdr, method=method, tax_id=tax_id)
                )

    if njob != 1:
        pool.close()
        pool.join()
        for k, j in jobs.items():
            de[k] = j.get(1e6)

    return de