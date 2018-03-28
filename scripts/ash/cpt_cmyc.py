import pandas as pd
import numpy as np
from scipy import stats
import os
from settings import DATA_DIR_NON_GIT
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import heatmap, clustering


def one_vs_many_correlation(one, many, method='pearson'):
    if method == 'pearson':
        the_func = stats.pearsonr
    elif method == 'spearman':
        the_func = stats.spearmanr
    else:
        raise NotImplementedError("Unrecognised method %s" % method)

    cor = pd.Series(index=many.index)
    pval = pd.Series(index=many.index)

    for idx, vals in many.iterrows():
        cor[idx], pval[idx] = the_func(one, vals)

    return cor, pval


if __name__ == '__main__':
    cmyc_probe = 3115514
    cmyc_gene = 'MYC'
    corr_method = 'pearson'

    fn = os.path.join(DATA_DIR_NON_GIT, 'exon_array', 'GSE60982', 'raw', 'GSE60892_HuEx-ALL40-Upload-Transposed.txt.gz')
    rma_data = pd.read_csv(fn, sep='\t', comment='#', header=0, index_col=0)

    ann_fn = os.path.join(DATA_DIR_NON_GIT, 'exon_array', 'GSE60982', 'HuEx-1_0-st-v2.na36.hg19.probeset.csv.gz')
    ann = pd.read_csv(ann_fn, sep=',', comment='#', header=0, index_col=0, na_values='---')

    ann_fn2 = os.path.join(DATA_DIR_NON_GIT, 'exon_array', 'ash_cpt_project', 'GSE60892_annotation.xlsx')
    ann2 = pd.read_excel(ann_fn2)

    # annotation files
    # latest version from Affymetrix / Thermofisher
    ann = ann.loc[rma_data.index]
    ann = ann.loc[~ann.gene_assignment.isnull()]
    ann.loc[:, 'gene_assignment'] = ann.loc[:, 'gene_assignment'].str.replace(r' /+ ', '|')
    symb = ann.gene_assignment.str.split('|').apply(lambda x: x[1])

    # version used by Jennie originally
    ann2 = ann2.loc[~ann2.probeset_id.isnull()]
    ann2.loc[:, 'probeset_id'] = ann2.loc[:, 'probeset_id'].astype(int)
    ann2.set_index("probeset_id", inplace=True)
    ann2 = ann2.loc[rma_data.index]
    ann2 = ann2.loc[~ann2.gene_assignment.isnull()]
    ann2.loc[:, 'gene_assignment'] = ann2.loc[:, 'gene_assignment'].str.replace(r' /+ ', '|')
    symb2 = ann2.gene_assignment.str.split('|').apply(lambda x: x[1] if len(x) > 1 else None)

    the_symbols = symb2

    # aggregate all genes
    all_aggr_by_gene = rma_data.groupby(the_symbols, axis=0).mean()

    # Pearson correlation by gene against the cymc probe
    base = rma_data.loc[cmyc_probe]
    cor, pval = one_vs_many_correlation(base, rma_data, method=corr_method)

    # select probes exhibiting correlation
    probes = cor.index[(cor.abs() > 0.5) & (pval < 0.05)]
    genes = the_symbols.loc[probes].dropna()
    print "%d probes selected that correlate with MYC, corresponding to %d unique genes" % (probes.size, genes.unique().size)

    # go back to the genes, taking the median probe aggregation
    glist = genes.unique()
    plist = the_symbols[the_symbols.isin(glist)].index

    aggr_by_gene = rma_data.loc[plist].groupby(the_symbols[plist], axis=0).mean()

    # re-run correlation with MYC
    cor_gene, pval_gene = one_vs_many_correlation(aggr_by_gene.loc[cmyc_gene], aggr_by_gene, method=corr_method)

    ### alternative approach: just take the gene signature and replot the heatmap

    jennie_list = [
        'MYC',
        'CXCL2',
        'CXCL1',
        'TNFAIP3',
        'IL8',
        'C17orf47',
        'RAET1L',
        'TEX14',
        'SERPINE1',
        'CCL2',
        'LIF',
        'TMEM49',
        'IL1RL2',
        'ASTL',
        'MSC',
        'CXCL3',
        'NFKBIZ',
        'FOSB',
        'PLK3',
        'G0S2',
        'CHST7',
        'GADD45B',
        'MAP2K3',
        'IL28A',
        'CXCL6',
        'ZNF165',
        'C10orf55',
        'PCK1',
        'VCAM1',
        'PLAUR',
        'B3GNT5',
        'CCNL1',
        'ATP13A3',
        'TIFA',
        'ATF3',
        'SEC14L4',
        'EIF2AK3',
        'IL6',
        'PPP1R15A',
        'WTAP',
        'ADM',
        'ZC3H12A',
        'JUN',
        'ALOXE3',
        'NBAG3',
        'FOSL1',
        'BIRC3',
        'CCL20',
        'PLAU',
        'RRAD',
        'KLF4',
        'NFKBID',
        'RHOB',
        'SFRP4',
        'CTNNA3',
        'DUSP2',
        'MAFF',
        'HSPA5',
        'SLC4A7',
        'SPRY2',
        'STX5',
        'CREB5',
        'SLC27A6',
        'WDR74',
        'ANKRD1',
        'NAMPT',
        'IER5',
        'ICAM1',
        'CCLF1',
        'TFP12',
        'SOD2',
        'REL',
        'SLFN5',
        'C14orf43',
        'WEE1',
        'CH25H',
        'NFIL3',
        'NGEF',
        'EIF2A',
        'IRF1',
        'EDN1',
        'CD44',
        'CD69',
        'HSPA6',
        'IL13RA2',
        'EMP1',
        'TGIF1',
        'ELF3',
        'C15orf48',
        'CSRNP1',
        'NFKBIA',
        'DMAP1',
        'NCOA7',
        'IL28A',
        'IL1B',
        'CCL4',
        'SLC2A3',
        'TACSTD2',
        'MT1M',
        'DUSP4',
        'DNAJB1',
        'C1orf94',
        'BACH2',
        'PTMA',
        'NLRP3',
        'CCL4',
        'KLF10',
        'HSPA7',
        'OFD1',
        'SLC2A14',
        'PRMT3',
        'ILDR1',
        'LRRC2',
        'ME1',
        'PTGS2',
        'DDX26B',
        'MT2A',
        'KBTBD11',
        'LDHC',
        'UBC',
        'BHLHE41',
        'ULBP2',
        'CCL3',
        'P1GA',
        'FHL1',
        'TSPYL2',
        'NR4A1',
        'JUNB',
        'FAM46B',
        'HES6',
        'DMGDH',
        'ZFP36',
        'EIF1AX',
        'PMP22',
        'MT1X',
        'DUSP1',
        'MT1B',
        'GPR183',
        'HBEGF',
        'FLJ43826',
        'KCNJ5',
        'ENDOU',
        'STAT5A',
        'IER2',
        'VPS54',
        'MT1L',
        'ARRDC3',
        'TMEM',
        'KDM6A',
        'THOC6',
        'NCEH1',
        'TRIB1',
        'SYCP2',
        'KLF6',
        'AOX1',
        'PAMR1',
        'FERMT3',
        'DGKZ',
        'ARL14',
        'C6orf141',
        'RND1',
        'TXK',
        'EHD1',
        'CAST',
        'LAMC2',
        'GK',
        'TAGAP',
        'NAV2',
        'SKIL',
        'MAGEA3',
        'CCDC132',
        'CCDC106',
        'LRP5L',
        'NFIC',
        'GOLT1A',
        'RCVRN',
        'PEX19',
        'LRRC20',
        'ABHD1',
        'RWDD1',
        'LRRN2',
        'PPFIA3',
        'CAPRIN2',
        'HPSE',
        'HDGFRP2',
        'BTN2A1',
        'SFTPB',
        'ADD1',
        'HIC1',
        'NRG2',
        'SDCCAG8',
        'NOVA2',
        'SCAI',
        'CLPP',
        'C9orf86',
        'PPP6R1',
        'ICK',
        'ZNF442',
        'HCG8',
        'EPB41L4B',
        'SIGLEC1',
        'CNNM2',
        'APLN',
        'SERF2',
        'RAB6B',
        'FLJ42875',
        'MLL2',
        'DSC3',
        'SMUG1',
        'ZBTB40',
        'SCAP',
        'DNAJC16',
        'RAB33A',
    ]

    # bizarrely, not all of these genes have a match - in fact, some don't even exist online?!
    # anyway, do the best I can
    lookup = the_symbols[the_symbols.isin(jennie_list)]
    j_data = rma_data.loc[lookup.index].groupby(lookup.values, axis=0).mean()

    # couple of different options here

    cg = clustering.plot_clustermap(j_data, cmap='RdBu_r', metric='euclidean', method='ward')
    cg = clustering.plot_clustermap(j_data, cmap='RdBu_r', metric='correlation', method='ward')