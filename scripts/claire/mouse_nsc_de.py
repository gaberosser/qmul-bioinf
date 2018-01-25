from load_data import rnaseq_data
from rnaseq import differential_expression
from utils.output import unique_output_dir
import references
import os
import pandas as pd
import numpy as np


if __name__ == '__main__':
    outdir = unique_output_dir("mouse_NSC_DE", reuse_empty=True)

    lfc = 0
    fdr = 0.01

    # load our data

    obj = rnaseq_data.mouse_nsc_validation_samples(annotate_by='Ensembl Gene ID')
    dat = obj.data
    dat = dat.loc[dat.index.str.contains('ENS')]
    idx = dat.columns.str.contains(r'eNSC[0-9]med') | dat.columns.str.contains(r'mDura[0-9AN]*human')
    dat = dat.loc[:, idx]
    the_groups = pd.Series('eNSC', index=dat.columns)
    the_groups[dat.columns.str.contains('mDura')] = 'iNSC'
    the_contrast = 'iNSC - eNSC'

    res = differential_expression.edger_glmqlfit(dat, the_groups, the_contrast, lfc=lfc, fdr=fdr, return_full=True)

    gs = references.ensembl_to_gene_symbol(res.index, tax_id=10090)
    gs = gs.loc[~gs.index.duplicated()]
    res.insert(0, "Gene symbol", gs.values)

    out_fn = os.path.join(outdir, "de_insc-ensc.xlsx")

    xl_writer = pd.ExcelWriter(out_fn)
    res.to_excel(xl_writer, "DE results", index=True)
    xl_writer.save()

    # load the GBM mouse model data and combine
    obj_gbm = rnaseq_data.mouse_gbm_pten_p53(annotate_by='Ensembl Gene ID')
    obj = rnaseq_data.MultipleBatchLoader([obj, obj_gbm])
    dat = obj.data
    dat = dat.loc[dat.index.str.contains('ENS')]

    idx = (
        dat.columns.str.contains(r'eNSC[0-9]med')
        | dat.columns.str.contains(r'mDura[0-9AN]*human')
        | dat.columns.str.contains(r'GBM')
        | dat.columns.str.contains(r'WT')
    )
    dat = dat.loc[:, idx]
    groups = pd.Series('eNSC', index=dat.columns)
    groups[dat.columns.str.contains('mDura')] = 'iNSC'
    groups[dat.columns.str.contains('GBM')] = 'GBM'
    groups[dat.columns.str.contains('WT')] = 'eNSC2'

    # expt 1: eNSC vs iNSC
    res = {}
    the_dat = dat.loc[:, groups.index[(groups == 'iNSC') | (groups == 'eNSC')]]
    the_groups = pd.Series('eNSC', index=the_dat.columns)
    the_groups[the_dat.columns.str.contains('mDura')] = 'iNSC'
    the_contrast = 'iNSC - eNSC'

    res['iNSC_vs_eNSC'] = differential_expression.edger_glmqlfit(the_dat, the_groups, the_contrast, lfc=lfc, fdr=fdr, return_full=True)