import itertools
import multiprocessing as mp
import os

import pandas as pd

from rnaseq import differential_expression, general, loader
from utils import excel
from utils.output import unique_output_dir


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

    # all samples
    # samples = ['mDura%shuman' % i for i in ('3N1', '5N24A', '6N6')]

    # drop mDura3N1 as it is significantly different
    samples = ['mDura%shuman' % i for i in ('5N24A', '6N6')]

    obj1 = loader.load_references('wtchg_p170390', source='star', tax_id=10090, samples=samples, strandedness='r')
    obj1_s = loader.load_references('wtchg_p170390', source='salmon', tax_id=10090, samples=samples)

    samples = ['eNSC%dmed' % i for i in (3, 5, 6)]
    obj2 = loader.load_references('wtchg_p170506', source='star', tax_id=10090, samples=samples, strandedness='r')
    obj2_s = loader.load_references('wtchg_p170506', source='salmon', tax_id=10090, samples=samples)

    # load the GBM mouse model data and combine
    samples = ['GBM209', 'GBM328', 'WT383', 'WT384']
    obj_gbm = loader.load_references(
        'GBM_Pten_P53', source='star', tax_id=10090, samples=samples, batch_names=['BR'], strandedness='u'
    )
    obj_gbm_s = loader.load_references('GBM_Pten_P53', source='salmon', tax_id=10090, samples=samples, batch_names=['BR'])

    # Ying / Brandner data
    # Drop the 3rd replicate in both GIC runs as it appears different from the other two
    samples = ['NSC repl %d' % i for i in range(1, 4)] + \
        ['P53-Pten GIC repl %d' % i for i in range(1, 3)] + \
        ['Idh1-P53-Pten GIC repl %d' % i for i in range(1, 3)]

    # All samples included
    # samples = ['NSC repl %d' % i for i in range(1, 4)] + \
    #     ['P53-Pten GIC repl %d' % i for i in range(1, 4)] + \
    #     ['Idh1-P53-Pten GIC repl %d' % i for i in range(1, 4)]

    obj_ying = loader.load_references(
        'brandner_mouse_gic', source='star', tax_id=10090, batch_names=['Ying'], strandedness='u', samples=samples
    )
    obj_ying_s = loader.load_references(
        'brandner_mouse_gic', source='salmon', tax_id=10090, batch_names=['Ying'], samples=samples
    )

    samples = ['GBM repl %d' % i for i in range(1, 6)]
    obj_gse73127 = loader.load_references(
        'GSE73127',
        source='star',
        tax_id=10090,
        batch_names=['GSE73127'],
        strandedness='r',
        samples=samples
    )
    obj_gse73127_s = loader.load_references(
        'GSE73127',
        source='salmon',
        tax_id=10090,
        batch_names=['GSE73127'],
        samples=samples
    )

    samples = ['Embryonic NSC repl 1', 'Embryonic NSC repl 2', 'Tumour']
    obj_gse64411 = loader.load_references(
        'GSE64411/trimgalore',
        source='star',
        tax_id=10090,
        batch_names=['GSE64411'],
        strandedness='u',
        samples=samples
    )
    obj_gse64411_s = loader.load_references(
        'GSE64411/trimgalore',
        source='salmon',
        tax_id=10090,
        batch_names=['GSE64411'],
        samples=samples
    )

    samples = ['GXP GBM%d' % i for i in range(1, 7)]
    obj_gse75300 = loader.load_references(
        'GSE75300',
        source='star',
        tax_id=10090,
        batch_names=['GSE75300'],
        strandedness='r',
        samples=samples
    )
    # TODO
    # obj_gse75300_s = loader.load_references(
    #     'GSE75300/trim_galore',
    #     source='salmon',
    #     tax_id=10090,
    #     batch_names=['GSE64411'],
    #     samples=samples
    # )

    obj = loader.MultipleBatchLoader([obj1, obj2, obj_gbm, obj_ying, obj_gse64411, obj_gse73127, obj_gse75300])
    obj_s = loader.MultipleBatchLoader([obj1_s, obj2_s, obj_gbm_s, obj_ying_s, obj_gse64411_s, obj_gse73127_s])

    dat = obj.data

    groups = pd.Series('eNSC', index=dat.columns)
    groups[dat.columns.str.contains('mDura')] = 'iNSC'
    groups[dat.columns.str.contains('GBM') & (obj.batch_id == 'BR')] = 'GIC_BR'
    groups[dat.columns.str.contains('WT') & (obj.batch_id == 'BR')] = 'eNSC_BR'
    groups[dat.columns.str.contains(r'NSC.*\(Ying\)')] = 'eNSC_Ying'
    groups[dat.columns.str.contains(r'^P53-Pten GIC')] = 'GIC_Ying'
    groups[dat.columns.str.contains(r'Idh1-P53-Pten GIC')] = 'IDH1_GIC_Ying'
    groups[dat.columns.str.contains('GBM') & (obj.batch_id == 'GSE73127')] = 'GBM_GSE73127'
    groups[dat.columns.str.contains('Tumour') & (obj.batch_id == 'GSE64411')] = 'GBM_GSE64411'
    groups[dat.columns.str.contains('NSC') & (obj.batch_id == 'GSE64411')] = 'ENSC_GSE64411'
    groups[dat.columns.str.contains('GBM') & (obj.batch_id == 'GSE75300')] = 'GBM_GSE75300'

    gic_lines = ['GIC_BR', 'GIC_Ying', 'IDH1_GIC_Ying']
    nsc_lines = ['eNSC', 'iNSC', 'eNSC_BR', 'eNSC_Ying']

    comparisons = list(itertools.product(gic_lines, nsc_lines)) + [
        ('iNSC', 'eNSC'),
        ('iNSC', 'eNSC_BR'),
        ('iNSC', 'eNSC_Ying'),
        ('eNSC', 'eNSC_BR'),
        ('eNSC', 'eNSC_Ying'),
        ('GBM_GSE73127', 'eNSC'),
        ('GBM_GSE73127', 'iNSC'),
        ('GBM_GSE64411', 'eNSC'),
        ('GBM_GSE64411', 'iNSC'),
        ('GBM_GSE64411', 'ENSC_GSE64411'),
        ('GBM_GSE75300', 'eNSC'),
        ('GBM_GSE75300', 'iNSC'),
    ]

    res = {}
    res_sign = {}
    res_lfc0 = {}

    jobs = {}
    pool = mp.Pool()

    for cmp in comparisons:
        lbl = "%s_vs_%s" % cmp
        jobs[lbl] = pool.apply_async(run_one_de, args=(dat, groups, cmp), kwds=de_params)
        # res[lbl] = run_one_de(dat, groups, cmp, **de_params)
        # print "%d DE genes\n" % (res[lbl].FDR <= de_params['fdr']).sum()

    for lbl in jobs:
        res[lbl] = jobs[lbl].get(1e6)
        print lbl
        print "%d DE genes\n" % (res[lbl].FDR <= de_params['fdr']).sum()

    for k, v in res.items():
        general.add_gene_symbols_to_ensembl_data(v, tax_id=10090)
        res_sign[k] = v.loc[v.FDR <= de_params['fdr']]

    excel.pandas_to_excel(res, os.path.join(outdir, "mouse_GBM_NSC_DE_all.xlsx"))
    excel.pandas_to_excel(res_sign, os.path.join(outdir, "mouse_GBM_NSC_DE_significant.xlsx"))

    # finally, re-run with a lfc of zero
    # disabled for now to speed things up
    if False:
        de_params['lfc'] = 0

        jobs2 = {}
        print "No logFC requirement"

        for cmp in comparisons:
            lbl = "%s_vs_%s" % cmp
            jobs2[lbl] = pool.apply_async(run_one_de, args=(dat, groups, cmp), kwds=de_params)

        for lbl in jobs2:
            res_lfc0[lbl] = jobs2[lbl].get(1e6)
            print lbl
            print "%d DE genes\n" % (res_lfc0[lbl].FDR <= de_params['fdr']).sum()
            general.add_gene_symbols_to_ensembl_data(res_lfc0[lbl], tax_id=10090)

        excel.pandas_to_excel(res_lfc0, os.path.join(outdir, "mouse_GBM_NSC_DE_all_nologfc.xlsx"))

    for_export = dat.copy()
    general.add_gene_symbols_to_ensembl_data(for_export, tax_id=10090)
    for_export.to_excel(os.path.join(outdir, 'gene_counts.xlsx'))

    for_export_s = obj_s.data.copy()
    general.add_gene_symbols_to_ensembl_data(for_export_s, tax_id=10090)
    for_export_s.to_excel(os.path.join(outdir, 'tpm_values.xlsx'))

