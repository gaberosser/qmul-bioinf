from utils import output
from rnaseq import gsea
from settings import DATA_DIR_NON_GIT, HGIC_LOCAL_DIR
import references
from scripts.wang_gbm_classifier.classify_our_gbm_samples import simplicity_score, load_pvalue_results

import pandas as pd
import os


if __name__ == '__main__':
    ## TCGA
    # pval cutoff to declare a match
    alpha = 0.05

    # rnaseq_type = 'counts'
    rnaseq_type = 'gliovis'

    # remove_idh1 = False
    remove_idh1 = True

    outdir = output.unique_output_dir()

    if rnaseq_type == 'counts':
        rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count')
        rnaseq_dat_fn = os.path.join(rnaseq_dir, 'counts.csv')
        rnaseq_meta_fn = os.path.join(rnaseq_dir, 'sources.csv')
        reader = pd.read_csv
    elif rnaseq_type == 'fpkm':
        rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count_fpkm')
        rnaseq_dat_fn = os.path.join(rnaseq_dir, 'fpkm.csv')
        rnaseq_meta_fn = os.path.join(rnaseq_dir, 'sources.csv')
        reader = pd.read_csv
    elif rnaseq_type == 'gliovis':
        rnaseq_dir = os.path.join(
            HGIC_LOCAL_DIR,
            'current/input_data/tcga/gliovis'
        )
        rnaseq_dat_fn = os.path.join(rnaseq_dir, 'gliovis_tcga_gbm_rnaseq.xlsx')
        rnaseq_meta_fn = os.path.join(rnaseq_dir, 'GlioVis_TCGA_GBMLGG.meta.xlsx')
        reader = pd.read_excel

        # rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm')
        # rnaseq_dat_fn = os.path.join(rnaseq_dir, 'gliovis_tcga_gbmlgg_expression.csv')
        # rnaseq_meta_fn = os.path.join(rnaseq_dir, 'gliovis_tcga_gbmlgg_meta.csv')
    else:
        raise NotImplementedError("Unrecognised rnaseq data type")


    rnaseq_dat_raw = reader(rnaseq_dat_fn, header=0, index_col=0)
    rnaseq_meta = reader(rnaseq_meta_fn, header=0, index_col=0)

    if rnaseq_type == 'gliovis':
        # filter only GBM
        rnaseq_meta = rnaseq_meta.loc[rnaseq_meta.Histology == 'GBM']
        # rnaseq_dat_raw = rnaseq_dat_raw.transpose().loc[:, rnaseq_meta.index]
        rnaseq_dat_raw = rnaseq_dat_raw.loc[:, rnaseq_meta.index]
        # add meta columns for compatibility
        idh1_status = pd.Series(data='Mut', index=rnaseq_meta.index, name='idh1_status')
        # idh1_status.loc[rnaseq_meta.loc[rnaseq_meta.loc[:, 'IDH_codel.subtype'] == 'IDHwt'].index] = 'WT'
        idh1_status.loc[rnaseq_meta.loc[rnaseq_meta.loc[:, 'IDH.status'] == 'WT'].index] = 'WT'
        rnaseq_meta.loc[:, 'idh1_status'] = idh1_status
        rnaseq_meta.loc[:, 'expression_subclass'] = rnaseq_meta.loc[:, 'Subtype.original']

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