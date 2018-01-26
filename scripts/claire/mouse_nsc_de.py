from load_data import rnaseq_data
from rnaseq import differential_expression, general
from utils.output import unique_output_dir
from utils import excel
import references
import os
import pandas as pd
import numpy as np


def run_one_de(dat, groups, contrasts, **kwargs):
    the_idx = groups.isin(contrasts)
    the_dat = dat.loc[:, groups.index[the_idx]]
    the_groups = groups.loc[the_idx]
    the_contrast = ' - '.join(contrasts)
    print "DE"
    print "(%s)\n vs\n(%s)" % (
        ", ".join(the_dat.columns[the_groups == contrasts[0]]),
        ", ".join(the_dat.columns[the_groups == contrasts[1]]),
    )
    return differential_expression.edger_glmqlfit(the_dat, the_groups, the_contrast, **kwargs)


if __name__ == '__main__':
    outdir = unique_output_dir("mouse_NSC_DE", reuse_empty=True)

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'return_full': True
    }

    # load our data
    obj = rnaseq_data.mouse_nsc_validation_samples(annotate_by='Ensembl Gene ID')

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

    res = {}

    # expt 1: eNSC vs iNSC
    res['iNSC_vs_eNSC'] = run_one_de(dat, groups, ('iNSC', 'eNSC'), **de_params)
    print "%d DE genes\n" % (res['iNSC_vs_eNSC'].FDR <= de_params['fdr']).sum()

    # expt 2: GBM (Pten/P53) vs WT NSC
    res['GBM_vs_eNSC_WT'] = run_one_de(dat, groups, ('GBM', 'eNSC2'), **de_params)
    print "%d DE genes\n" % (res['GBM_vs_eNSC_WT'].FDR <= de_params['fdr']).sum()

    # expt 3: GBM (Pten/P53) vs our eNSC
    res['GBM_vs_eNSC'] = run_one_de(dat, groups, ('GBM', 'eNSC'), **de_params)
    print "%d DE genes\n" % (res['GBM_vs_eNSC'].FDR <= de_params['fdr']).sum()

    # expt 4: GBM (Pten/P53) vs our iNSC
    res['GBM_vs_iNSC'] = run_one_de(dat, groups, ('GBM', 'iNSC'), **de_params)
    print "%d DE genes\n" % (res['GBM_vs_iNSC'].FDR <= de_params['fdr']).sum()

    # expt 5: Our eNSC vs WT eNSC
    res['eNSC_vs_eNSC_WT'] = run_one_de(dat, groups, ('eNSC', 'eNSC2'), **de_params)
    print "%d DE genes\n" % (res['eNSC_vs_eNSC_WT'].FDR <= de_params['fdr']).sum()

    res_sign = {}
    for k, v in res.items():
        general.add_gene_symbols_to_ensembl_data(v, tax_id=10090)
        res_sign[k] = v.loc[v.FDR <= de_params['lfc']]

    excel.pandas_to_excel(res, os.path.join(outdir, "mouse_GBM_NSC_DE_all.xlsx"))
    excel.pandas_to_excel(res_sign, os.path.join(outdir, "mouse_GBM_NSC_DE_significant.xlsx"))

    for_export = dat.copy()
    general.add_gene_symbols_to_ensembl_data(for_export, tax_id=10090)
    for_export.to_excel(os.path.join(outdir, 'gene_counts.xlsx'))