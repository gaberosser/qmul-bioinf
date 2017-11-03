import pandas as pd
from rnaseq import gsea
from load_data import rnaseq_data
from utils.output import unique_output_dir
import os
import references


def prepare_gct_files(outdir=None):
    """
    Prepare the GCT files required to perform classification:
    - Our GBM FFPE and cell culture samples
    - TCGA RNA-Seq cohort
    - Both combined
    In all cases, use FPKM units and gene symbols, as these are used by Wang
    """
    if outdir is None:
        outdir = unique_output_dir("gct_files_for_wang")

    infiles = []

    # 1) Our data
    obj_ffpe = rnaseq_data.load_by_patient('all', type='ffpe')
    dat_ffpe = obj_ffpe.get_fpkm()
    dat_ffpe.columns = ['%s_FFPE' % t for t in obj_ffpe.meta.reference_id]
    obj_cc = rnaseq_data.load_by_patient(patient_ids='all')
    dat_cc = obj_cc.get_fpkm()
    dat_cc = dat_cc.loc[:, obj_cc.meta.type == 'GBM']
    dat_all = pd.concat((dat_cc, dat_ffpe), axis=1)
    idx = references.ensembl_to_gene_symbol(dat_all.index).dropna()
    dat_all = dat_all.loc[idx.index]
    dat_all.index = idx
    fn = os.path.join(outdir, "gbm_ffpe_cc_fpkm.gct")
    gsea.data_to_gct(dat_all, fn)
    infiles.append(fn)

    # 2) TCGA (IDH1 WT only)
    tcga_dat, tcga_meta = rnaseq_data.tcga_primary_gbm(units='fpkm')
    tcga_dat = tcga_dat.loc[:, tcga_meta.idh1_status == 'WT']
    idx = references.ensembl_to_gene_symbol(tcga_dat.index).dropna()
    idx = idx.loc[~idx.index.duplicated()]
    tcga_dat = tcga_dat.loc[idx.index]
    tcga_dat.index = idx
    fn = os.path.join(outdir, "tcga_idh1_wt_fpkm.gct")
    gsea.data_to_gct(tcga_dat, fn)
    infiles.append(fn)

    # 3) Combined
    dat = gsea.combine_gct_files(*infiles)
    fn = os.path.join(outdir, "tcga_idh1_wt_and_gbm_ffpe_cc_fpkm.gct")
    gsea.data_to_gct(tcga_dat, fn)


if __name__ == '__main__':
    pass