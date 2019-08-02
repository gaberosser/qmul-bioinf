from utils import output
from rnaseq import gsea
from settings import DATA_DIR, HGIC_LOCAL_DIR
import references
from scripts.wang_gbm_classifier.classify_our_gbm_samples import simplicity_score, load_pvalue_results

import pandas as pd
import os


if __name__ == '__main__':
    ## TCGA
    # pval cutoff to declare a match
    alpha = 0.05

    # rnaseq_type = 'counts'
    rnaseq_type = 'fpkm'
    # rnaseq_type = 'gliovis'

    # remove_idh1 = False
    remove_idh1 = True

    outdir = output.unique_output_dir()

    basedir = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/tcga'
    )

    brennan_s7_fn = os.path.join(basedir, "brennan_s7.csv")
    brennan_s7 = pd.read_csv(brennan_s7_fn, header=0, index_col=0)

    if rnaseq_type == 'counts':
        rnaseq_dat_fn = os.path.join(basedir, 'rnaseq.xlsx')
        rnaseq_meta_fn = os.path.join(basedir, 'rnaseq.meta.xlsx')
        sheet_name = 'htseq'
    elif rnaseq_type == 'fpkm':
        rnaseq_dat_fn = os.path.join(basedir, 'rnaseq.xlsx')
        rnaseq_meta_fn = os.path.join(basedir, 'rnaseq.meta.xlsx')
        sheet_name = 'fpkm'
    elif rnaseq_type == 'gliovis':
        rnaseq_dat_fn = os.path.join(basedir, 'gliovis', 'gliovis_tcga_gbm_rnaseq.xlsx')
        rnaseq_meta_fn = os.path.join(basedir, 'gliovis', 'GlioVis_TCGA_GBMLGG.meta.xlsx')
        sheet_name = 0
    else:
        raise NotImplementedError("Unrecognised rnaseq data type")


    rnaseq_dat_raw = pd.read_excel(rnaseq_dat_fn, header=0, index_col=0, sheet_name=sheet_name)
    rnaseq_meta = pd.read_excel(rnaseq_meta_fn, header=0, index_col=0)

    if rnaseq_type == 'gliovis':
        # filter only GBM
        rnaseq_meta = rnaseq_meta.loc[rnaseq_meta.Histology == 'GBM']
        rnaseq_dat_raw = rnaseq_dat_raw.loc[:, rnaseq_meta.index]

        idh1_status = pd.Series(data='Mut', index=rnaseq_meta.index, name='idh1_status')
        idh1_status.loc[rnaseq_meta.loc[rnaseq_meta.loc[:, 'IDH.status'] == 'WT'].index] = 'WT'
        rnaseq_meta.insert(0, 'idh1_status', idh1_status)
    else:
        # simplify sample naming
        new_cols = rnaseq_dat_raw.columns.str.replace(r'(?P<x>TCGA-[0-9]{2}-[0-9]{4})-.*', r'\g<x>')
        rnaseq_meta = rnaseq_meta.loc[~new_cols.duplicated()]
        rnaseq_dat_raw = rnaseq_dat_raw.loc[:, rnaseq_meta.index]
        rnaseq_meta.index = new_cols[~new_cols.duplicated()]
        rnaseq_dat_raw.columns = rnaseq_meta.index

        # add IDH1 status from Brennan metadata
        idh1_status = brennan_s7.reindex(rnaseq_dat_raw.columns)['idh1_status']
        rnaseq_meta.insert(0, 'idh1_status', idh1_status)

    if remove_idh1:
        # filter IDH1 mutants
        idh1_wt = (~rnaseq_meta.idh1_status.isnull()) & (rnaseq_meta.idh1_status == 'WT')

        rnaseq_meta = rnaseq_meta.loc[idh1_wt]
        rnaseq_dat = rnaseq_dat_raw.loc[:, rnaseq_meta.index]
    else:
        rnaseq_dat = rnaseq_dat_raw.loc[:, rnaseq_dat_raw.columns.str.contains('TCGA')]

    if rnaseq_type != 'gliovis':
        # add gene symbols for gene signature scoring?
        gs = references.ensembl_to_gene_symbol(rnaseq_dat.index).dropna()
        rnaseq_dat = rnaseq_dat.loc[gs.index]
        rnaseq_dat.index = gs.values

    if rnaseq_type == 'counts':
        # convert to CPM
        rnaseq_dat = rnaseq_dat.divide(rnaseq_dat.sum(axis=0), axis=1) * 1e6

    fn = os.path.join(outdir, "tcga_%s.gct" % rnaseq_type)
    gsea.data_to_gct(rnaseq_dat, fn)

    gsea.wang_ssgsea_classification(fn)
    the_dir, the_stem = os.path.split(fn)
    outfn = os.path.join(the_dir, "p_result_%s.txt" % the_stem)

    pval = load_pvalue_results(outfn)
    ss = simplicity_score(pval)

    # create one more output file, containing the best match for each sample
    num_match = (pval < alpha).sum(axis=1)
    for_export = pd.DataFrame(index=pval.index, columns=['Wang subclass', 'Number of matches'])
    for_export.loc[:, 'Number of matches'] = num_match
    min_match = pval.idxmin(axis=1)
    for_export.loc[num_match == 1, 'Wang subclass'] = min_match.loc[num_match == 1]

    for row in num_match.index[num_match > 1]:
        this = pval.loc[row]
        for_export.loc[row, 'Wang subclass'] = ','.join(this.index[this < alpha])
    for_export.insert(2, 'Simplicity score', ss)

    for_export.to_csv(os.path.join(outdir, "tcga_%s_wang_classification.csv" % rnaseq_type))