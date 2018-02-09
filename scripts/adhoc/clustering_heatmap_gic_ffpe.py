from rnaseq import loader
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import clustering
import pandas as pd
import numpy as np
from utils import output
import os
import references
from stats import transformations
import hgic_consts


if __name__ == '__main__':
    outdir = output.unique_output_dir("cluster_gic_ffpe")
    pids = ['018', '019', '031', '017', '050', '054']
    obj_ffpe = loader.load_by_patient(pids, type='ffpe', include_control=False)
    obj_gic = loader.load_by_patient(pids, type='cell_culture', include_control=False)
    obj = loader.loader.MultipleBatchLoader([obj_ffpe, obj_gic])

    # drop iNSC, iPSC
    obj.meta = obj.meta.loc[~obj.meta.index.str.contains('DURA')]
    obj.data = obj.data.loc[:, obj.meta.index]

    # relabel the FFPE samples
    idx = obj.meta.index.tolist()
    for k, v in hgic_consts.NH_TO_PATIENT_ID.items():
        for i, t in enumerate(idx):
            if k in t:
                idx[i] = "FFPE GBM%s" % v
    obj.meta.index = idx
    obj.data.columns = idx

    cpm = obj.data.divide(obj.data.sum(), axis=1)
    lcpm = np.log2((obj.data + 1).divide((obj.data + 1).sum(), axis=1))

    mad = transformations.median_absolute_deviation(lcpm).sort_values(ascending=False)

    cg = clustering.plot_clustermap(lcpm.loc[mad.index[:3000]], cmap='RdBu_r')
    cg.gs.update(bottom=0.15)
    cg.savefig(os.path.join(outdir, "cluster_by_top_3000_genes.png"), dpi=200)

    # load Verhaak signatures
    # manual amendments:
    # Classical. C14orf159 -> DGLUCY, KIAA0494 -> EFCAB14, LHFP -> LHFPL6
    # Proneural. HN1 -> JPT1, PAK7 -> PAK5, ZNF643 -> ZFP69B

    cl = ['PTPRA', 'ELOVL2', 'MLC1', 'SOX9', 'ARNTL', 'DENND2A', 'BBS1', 'ABLIM1', 'PAX6', 'ZHX3', 'USP8', 'PLCG1', 'CDH4',
     'RASGRP1', 'ACSBG1', 'CST3', 'BCKDHB', 'LHFPL6', 'VAV3', 'ACSL3', 'EYA2', 'SEPT11', 'SLC4A4', 'SLC20A2', 'DGLUCY',
     'CTNND1', 'ZFHX4', 'SPRY2', 'ZNF45', 'NCOA1', 'PLCE1', 'DTNA', 'POLRMT', 'SALL1', 'TYK2', 'TJP1', 'MEOX2', 'FGFR3',
     'STXBP3', 'GRIK1', 'GATM', 'UPF1', 'NPEPL1', 'EFCAB14', 'RBCK1', 'PHKB', 'SLC3A2', 'PPARGC1A', 'PNPLA6', 'MYO5C']

    mes = ['ARPC1B', 'S100A11', 'CTSC', 'GLIPR1', 'NNMT', 'VDR', 'RGS2', 'CTSB', 'TGFBI', 'PLAUR', 'LY96', 'BCL3', 'TNFAIP8',
     'IER3', 'PRSS23', 'IL7R', 'RAB27A', 'RUNX1', 'P4HA2', 'CYP1B1', 'BACE2', 'ACPP', 'FTL', 'SLPI', 'RAC2', 'RARRES1',
     'SYNGR2', 'THBS1', 'IL6', 'CAV1', 'PI3', 'CDCP1', 'ITGB1', 'LOX', 'CD72', 'COL1A2', 'ANPEP', 'MMP7', 'SPAG4',
     'BNC2', 'NDRG1', 'CNN2', 'LUM', 'PTGS2', 'COL3A1', 'COL5A1', 'SDC1', 'COL1A1', 'GPRC5A', 'COL15A1']

    pn = ['JPT1', 'RAB33A', 'HDAC2', 'MYT1', 'MTSS1', 'HOXD3', 'GPR17', 'PTTG1', 'KLRC3', 'HRASLS', 'TCP1', 'NPPA', 'PFDN2',
     'CA10', 'EPHB1', 'UGT8', 'PAK5', 'SLC1A1', 'NARF', 'DCTN3', 'SMPD3', 'ZNF804A', 'RASL11B', 'MYB', 'PDGFRA',
     'ERBB3', 'CLGN', 'SOX10', 'BCL11A', 'NMU', 'ZFP69B', 'CDKN1C', 'JPH3', 'PCDHA9', 'IL1RAPL1', 'MAST1', 'VIPR2',
     'SIM2', 'BAMBI', 'PKMYT1', 'PLCB4', 'SLC17A6', 'KLRK1', 'CENPJ', 'NHLH1', 'GABRB3', 'KLRC4', 'KCNK3', 'GRID2',
     'DACH1']

    lcpm_gene = references.translate_quantification_resolving_duplicates(
        lcpm,
        from_field='Ensembl Gene ID',
        to_field = 'Approved Symbol'
    )

    cg = clustering.plot_clustermap(lcpm_gene.loc[cl + mes + pn], cmap='RdBu_r')
    cg.gs.update(bottom=0.15)
    cg.savefig(os.path.join(outdir, "cluster_by_verhaak signature.png"), dpi=200)

    # try comparing every FFPE to all GBMs
    ng = 3000
    meth = 'pearson'

    cc = lcpm.loc[mad.index[:ng]].corr(method=meth)
    cc = cc.loc[cc.index.str.contains('FFPE'), ~cc.columns.str.contains('FFPE')]
    cc = cc.loc[~cc.index.str.contains('GBM050')]
    cc = cc.loc[:, ~cc.columns.str.contains('GBM050')]
    fig = plt.figure(figsize=(9, 4.5))
    ax = fig.add_subplot(111)
    ax = sns.heatmap(cc, cmap='RdBu_r', ax=ax)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    ax.set_aspect('equal')
    ax.figure.tight_layout()
    fig.savefig(os.path.join(outdir, "%s_correlation_top_%d_genes.png" % (meth, ng)), dpi=200)
    fig.savefig(os.path.join(outdir, "%s_correlation_top_%d_genes.tiff" % (meth, ng)), dpi=200)

    # same but using only the Verhaak genes
    cc = lcpm_gene.loc[cl + mes + pn].corr(method=meth)
    cc = cc.loc[cc.index.str.contains('FFPE'), ~cc.columns.str.contains('FFPE')]
    cc = cc.loc[~cc.index.str.contains('GBM050')]
    cc = cc.loc[:, ~cc.columns.str.contains('GBM050')]
    fig = plt.figure(figsize=(8, 4.5))
    ax = fig.add_subplot(111)
    ax = sns.heatmap(cc, cmap='RdBu_r', ax=ax)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    ax.set_aspect('equal')
    ax.figure.tight_layout()
    fig.savefig(os.path.join(outdir, "%s_correlation_verhaak_genes.png" % meth), dpi=200)
