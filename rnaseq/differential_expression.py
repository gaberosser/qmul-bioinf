import pandas as pd
import general
import multiprocessing as mp
from utils import setops
from utils import rinterface
from rpy2 import robjects
from rpy2.robjects import pandas2ri, r
pandas2ri.activate()

toptags_cols = [
    "logFC", "unshrunk.logFC", "logCPM", "PValue", "FDR"
]


def model_matrix_from_formula(formula, **vars):
    formula = robjects.Formula(formula)
    for k, v in vars.items():
        ## TODO: this only works with discrete groups, but edgeR also supports continuous variables
        formula.environment[k] = robjects.FactorVector(v)
    return r("model.matrix")(formula)


def _edger_func_compute_DGEList_with_dispersion(the_data, the_formula=None, common_disp=False, **vars):
    if the_formula is None and not common_disp:
        raise AttributeError("Either supply a formula OR specify common_disp = True.")
    rdata = pandas2ri.py2ri(the_data)

    y = r("DGEList")(rdata)
    y = r("calcNormFactors")(y)

    if the_formula is not None:
        design = model_matrix_from_formula(the_formula, **vars)
        y = r("estimateDisp")(y, design)
        return y, design

    else:
        # use a common estimate of the dispersion rather than using experimental structure
        # this is helpful where we have no replicates
        y = r("estimateGLMCommonDisp")(y, method='deviance', robust=True, subset=robjects.NULL)
        return y


edger_dgelist = rinterface.RFunctionDeferred(_edger_func_compute_DGEList_with_dispersion, imports=['edgeR'])


def _edger_func_fit_to_dgelist(dgelist, the_method, design=None, the_formula=None, **vars):
    if design is None and the_formula is None:
        raise AttributeError("Must either supply a design matrix OR a formula")
    if the_method not in {'GLM', 'QLGLM'}:
        raise NotImplementedError("Only GLM and QLGLM methods are supported at present")

    if design is None:
        design = model_matrix_from_formula(the_formula, **vars)

    fit = None
    if the_method == 'GLM':
        fit = r('glmFit')(dgelist, design)
    elif the_method == 'QLGLM':
        fit = r('glmQLFit')(dgelist, design)

    return fit, design


edger_fit_dgelist = rinterface.RFunctionDeferred(_edger_func_fit_to_dgelist, imports=['edgeR'])


def _edger_func_fit_glm(the_data, the_method, the_formula, common_disp=False, **vars):
    if the_method not in {'GLM', 'QLGLM'}:
        raise NotImplementedError("Only GLM and QLGLM methods are supported at present")
    fit = None
    rdata = pandas2ri.py2ri(the_data)

    formula = robjects.Formula(the_formula)
    for k, v in vars.items():
        formula.environment[k] = robjects.FactorVector(v)

    y = r("DGEList")(rdata)
    y = r("calcNormFactors")(y)
    design = r("model.matrix")(formula)

    if common_disp:
        # use a common estimate of the dispersion rather than using experimental structure
        # this is helpful where we have no replicates
        y = r("estimateGLMCommonDisp")(y, method='deviance', robust=True, subset=robjects.NULL)
    else:
        y = r("estimateDisp")(y, design)
    if the_method == 'GLM':
        fit = r('glmFit')(y, design)
    elif the_method == 'QLGLM':
        fit = r('glmQLFit')(y, design)
    return fit, design


edger_fit_glm = rinterface.RFunctionDeferred(_edger_func_fit_glm, imports=['edgeR'])


def _edger_func_fit(the_data, the_groups, the_method):
    if the_method not in {'GLM', 'QLGLM'}:
        raise NotImplementedError("Only GLM and QLGLM methods are supported at present")
    fit = None
    rdata = pandas2ri.py2ri(the_data)
    rgroups = robjects.FactorVector(the_groups)
    y = r("DGEList")(rdata)
    y = r("calcNormFactors")(y)
    formula = robjects.Formula("~0 + groups")
    formula.environment['groups'] = rgroups
    design = r("model.matrix")(formula)
    design.colnames = r('levels')(rgroups)
    y = r("estimateDisp")(y, design)
    if the_method == 'GLM':
        fit = r('glmFit')(y, design)
    elif the_method == 'QLGLM':
        fit = r('glmQLFit')(y, design)
    return fit, design


edger_fit = rinterface.RFunctionDeferred(_edger_func_fit, imports=['edgeR'])


def _edger_func_test(fit, design, contrast_str, fdr=0.01, lfc=1, return_full=False):
    rcontrast = r('makeContrasts')(contrast_str, levels=design)
    lrt = r('glmTreat')(fit, contrast=rcontrast, lfc=lfc)
    if return_full:
        toptags = r('topTags')(lrt, n=r('Inf'), **{'p.value': 1.})
    else:
        toptags = r('topTags')(lrt, n=r('Inf'), **{'p.value': fdr})
    if len(toptags) == 0:
        return pd.DataFrame(columns=toptags_cols)
    else:
        return pandas2ri.ri2py_dataframe(toptags[toptags.names.index('table')])


edger_test = rinterface.RFunctionDeferred(_edger_func_test, imports=['edgeR'])


def _edger_func_glmqlfit(the_data, the_groups, the_contrast, fdr=0.01, lfc=1, return_full=False):
    fit, design = _edger_func_fit(the_data, the_groups, 'QLGLM')
    return _edger_func_test(fit, design, the_contrast, fdr=fdr, lfc=lfc, return_full=return_full)


edger_glmqlfit = rinterface.RFunctionDeferred(_edger_func_glmqlfit, imports=['edgeR'])


def _edger_func_glmfit(the_data, the_groups, the_contrast, fdr=0.01, lfc=1, return_full=False):
    fit, design = _edger_func_fit(the_data, the_groups, 'GLM')
    return _edger_func_test(fit, design, the_contrast, fdr=fdr, lfc=lfc, return_full=return_full)


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


def run_one_de(
        the_data,
        the_groups,
        the_comparison,
        lfc=1,
        fdr=0.01,
        method='QLGLM',
        return_full=False,
        tax_id=9606,
        add_gene_symbols=True,
):
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

    if add_gene_symbols:
        general.add_gene_symbols_to_ensembl_data(res, tax_id=tax_id)
    general.add_fc_direction(res)

    return res


def run_multiple_de(
        the_data,
        the_groups,
        comparisons,
        lfc=1,
        fdr=0.01,
        method='QLGLM',
        return_full=False,
        tax_id=9606,
):
    """
    :param comparisons: Dictionary of comparisons, each of which is a string formula, e.g. "groupB - groupA".
    The keys are used to organise the resulting dictionary.
    :return: Dict with same keys as comparisons containing DE result tables.
    """
    if method not in {'GLM', 'QLGLM'}:
        raise NotImplementedError("Only GLM and QLGLM methods are supported at present")

    fit, design = edger_fit(the_data, the_groups, method)
    res = {}

    for k, contrast_str in comparisons.items():
        res[k] = edger_test(fit, design, contrast_str, fdr=fdr, lfc=lfc, return_full=return_full)

    for k in res:
        general.add_gene_symbols_to_ensembl_data(res[k], tax_id=tax_id)
        general.add_fc_direction(res[k])

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
        return_full=False,
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

    if njob is None:
        njob = mp.cpu_count()

    if njob != 1:
        pool = mp.Pool(processes=njob)
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
                    tax_id=tax_id,
                    return_full=return_full
                )
            else:
                jobs[(pid, pid2)] = pool.apply_async(
                    run_one_de,
                    args=(the_data, the_groups, the_comparison),
                    kwds=dict(lfc=lfc, fdr=fdr, method=method, tax_id=tax_id, return_full=return_full)
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
                    tax_id=tax_id,
                    return_full=return_full,
                )
            else:
                jobs[(pid, er)] = pool.apply_async(
                    run_one_de,
                    args=(the_data, the_groups, the_comparison),
                    kwds=dict(lfc=lfc, fdr=fdr, method=method, tax_id=tax_id, return_full=return_full)
                )

    if njob != 1:
        pool.close()
        pool.join()
        for k, j in jobs.items():
            de[k] = j.get(1e6)

    return de


def venn_set_to_dataframe(
        data,
        venn_set,
        set_labels,
        include_sets=None,
        full_data=None,
        logfc_col='logFC',
        fdr_col='FDR',
        run_sanity_check=False,
        add_null_set=False,
):
    """
    Given the input DE data and Venn sets, generate a wide format dataframe containing all the data, one column
    per patient and one row per gene.
    Optionally filter the sets to include only a subset.
    Optionally include non-significant results too.
    :param data: Dict containing DE results, keyed by the entries of set_labels
    :param venn_set:
    :param set_labels:
    :param include_sets:
    :param full_data: If supplied, this has the same format as `data`, but the lists are complete so that even non-
    significant results can be accessed.
    :param logfc_col: The name of the log fold change column in the input data. Also used to name columns in the df.
    :param fdr_col: The name of the FDR column in the input data. Also used to name columns in the df.
    :param run_sanity_check: (default: False) If True, run an additional sanity check at the end. This *should* be
    unnecessary. It's slow for larger numbers of members.
    :return:
    """
    if add_null_set and full_data is None:
        raise ValueError("Can only add_null_set if full_data is supplied.")
    if include_sets is not None:
        venn_set = dict([
            (k, v) for k, v in venn_set.items() if k in include_sets
        ])


    # precompute columns
    cols = reduce(
        lambda x, y: x + y,
        [[t, "%s_%s" % (t, logfc_col), "%s_%s" % (t, fdr_col)] for t in set_labels]
    ) + ['consistency']

    res = []
    genes_seen = set()
    for k in venn_set:
        the_genes = venn_set[k]
        genes_seen.update(the_genes)

        # populate with individual patient results
        this_block = pd.DataFrame(index=the_genes, columns=cols)
        # blocks = []
        consistency_check = []
        for i, t in enumerate(k):
            pid = set_labels[i]

            if t == '1':
                this_block.loc[:, pid] = 'Y'
                this_block.loc[the_genes, "%s_%s" % (pid, logfc_col)] = data[pid].loc[the_genes, logfc_col]
                this_block.loc[the_genes, "%s_%s" % (pid, fdr_col)] = data[pid].loc[the_genes, fdr_col]

                cc = data[pid].loc[the_genes, 'Direction']
                cc.name = pid
                consistency_check.append(cc)
            else:
                this_block.loc[:, pid] = 'N'

                # this_datum.loc[the_genes, pid] = 'N'
                if full_data is not None:
                    # we can't guarantee there will be entries for all genes, as filtering removes some
                    # therefore find matches in advance and only fill in those rows
                    the_genes_present = pd.Index(the_genes).intersection(full_data[pid].index)
                    this_block.loc[the_genes_present, "%s_%s" % (pid, logfc_col)] = full_data[pid].loc[the_genes_present, logfc_col]
                    this_block.loc[the_genes_present, "%s_%s" % (pid, fdr_col)] = full_data[pid].loc[the_genes_present, fdr_col]

        # assess consistency of DE direction
        consist = pd.Series(index=the_genes)

        if len(consistency_check) > 0:
            consistency_check = pd.concat(consistency_check, axis=1)
            idx = consistency_check.apply(lambda col: col == consistency_check.iloc[:, 0]).all(axis=1)
            consist.loc[idx] = 'Y'
            consist.loc[~idx] = 'N'

        this_block.loc[:, 'consistency'] = consist

        res.append(this_block)

    # check: no genes should be in more than one data entry
    if run_sanity_check:
        for i, k in enumerate(venn_set):
            for j, k2 in enumerate(venn_set):
                if k == k2: continue
                bb = len(res[i].index.intersection(res[j].index))
                if bb > 0:
                    raise AttributeError("Identified %d genes that are in BOTH %s and %s" % (bb, k, k2))

    if add_null_set:
        all_genes = setops.reduce_union(*[t.index for t in full_data.values()])
        add_genes = all_genes.difference(genes_seen)
        this_block = pd.DataFrame(index=add_genes, columns=cols)
        for pid in set_labels:
            # by definition, no samples are DE positive in the null set
            this_block.loc[:, pid] = 'N'
            the_genes_present = add_genes.intersection(full_data[pid].index)
            this_block.loc[the_genes_present, "%s_%s" % (pid, logfc_col)] = full_data[pid].loc[the_genes_present, logfc_col]
            this_block.loc[the_genes_present, "%s_%s" % (pid, fdr_col)] = full_data[pid].loc[the_genes_present, fdr_col]
        res.append(this_block)

    res = pd.concat(res, axis=0)

    # add gene symbols
    general.add_gene_symbols_to_ensembl_data(res)

    return res