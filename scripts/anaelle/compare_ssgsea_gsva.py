from rnaseq.gsea import ssgsea, run_one_ssgsea
from rnaseq import loader, gsva
import pandas as pd
from settings import DATA_DIR, GIT_LFS_DATA_DIR, HGIC_LOCAL_DIR
import os
import csv
import references
import datetime
from matplotlib import pyplot as plt
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import seaborn as sns
import numpy as np
import collections
from scipy import stats

from utils.output import unique_output_dir
from scripts.hgic_final import consts


if __name__ == '__main__':
    # Step 1: compare RNA-Seq count data and some known gene signatures
    # We just want to demonstrate that ssGSEA and GSVA are similar
    outdir = unique_output_dir()

    ipa_pathway_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/ipa_pathways',
        'ipa_exported_pathways_ensembl_ids.csv'
    )
    ipa_pathways = {}
    with open(ipa_pathway_fn, 'rb') as f:
        c = csv.reader(f)
        for row in c:
            ipa_pathways[row[0]] = row[2:]

    # load TCGA count data (for GSVA)
    tcga_dir = os.path.join(DATA_DIR, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count')
    tcga_dat_fn = os.path.join(tcga_dir, 'counts.csv')
    tcga_meta_fn = os.path.join(tcga_dir, 'sources.csv')

    tcga_counts = pd.read_csv(tcga_dat_fn, header=0, index_col=0)
    tcga_gene_symbol = tcga_counts[['Approved Symbol']]
    tcga_counts.drop('Approved Symbol', axis=1, inplace=True)
    tcga_meta = pd.read_csv(tcga_meta_fn, header=0, index_col=0)

    # load TCGA FPKM data (for ssGSEA)
    tcga_fpkm_dir = os.path.join(DATA_DIR, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count_fpkm')
    tcga_fpkm_fn = os.path.join(tcga_fpkm_dir, 'fpkm.csv')
    tcga_fpkm = pd.read_csv(tcga_fpkm_fn, header=0, index_col=0).drop('Approved Symbol', axis=1)

    # for some reason, the samples are slightly different in the two TCGA datasets (more in the FPKM set)
    # harmonise now
    common_ix = tcga_counts.columns.intersection(tcga_fpkm.columns)
    tcga_counts = tcga_counts[common_ix]
    tcga_fpkm = tcga_fpkm[common_ix]
    tcga_meta = tcga_meta.loc[common_ix]

    # remove IDH1 mut
    ix = tcga_meta.idh1_status == 'WT'
    tcga_counts = tcga_counts.loc[:, ix]
    tcga_fpkm = tcga_fpkm.loc[:, ix]

    # run ssGSEA
    ssgsea_tcga = ssgsea(tcga_fpkm, ipa_pathways)
    ssgsea_tcga.to_csv(os.path.join(outdir, "tcga_fpkm_ssgsea.csv"))

    # run GSVA on count data
    gsva_obj = gsva.GSVAForCounts(tcga_counts, gene_sets=ipa_pathways)
    gsva_res = gsva_obj.run_enrichment()
    gsva_res['es1_jk'].to_csv(os.path.join(outdir, "tcga_counts_gsva_es1.csv"))
    gsva_res['es2_jk'].to_csv(os.path.join(outdir, "tcga_counts_gsva_es2.csv"))
    gsva_res['es3_jk'].to_csv(os.path.join(outdir, "tcga_counts_gsva_es3.csv"))

    # run GSVA on FPKM data
    gsva_obj2 = gsva.GSVAForNormedData(tcga_fpkm, gene_sets=ipa_pathways)
    gsva_res2 = gsva_obj2.run_enrichment()
    gsva_res2['es1_jk'].to_csv(os.path.join(outdir, "tcga_fpkm_gsva_es1.csv"))
    gsva_res2['es2_jk'].to_csv(os.path.join(outdir, "tcga_fpkm_gsva_es2.csv"))
    gsva_res2['es3_jk'].to_csv(os.path.join(outdir, "tcga_fpkm_gsva_es3.csv"))

    # load (our) data
    obj_star = loader.load_by_patient(consts.PIDS, source='star', include_control=False)
    obj_salmon = loader.load_by_patient(consts.PIDS, source='salmon', include_control=False)


    pass