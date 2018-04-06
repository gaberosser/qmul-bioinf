import pandas as pd
import numpy as np
from scipy import stats
import os
from settings import DATA_DIR_NON_GIT
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import heatmap, clustering
from utils import setops


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


def pairwise_correlation(many, method='pearson'):
    if method == 'pearson':
        the_func = stats.pearsonr
    elif method == 'spearman':
        the_func = stats.spearmanr
    else:
        raise NotImplementedError("Unrecognised method %s" % method)

    cor = pd.DataFrame(index=many.index, columns=many.index)
    pval = pd.DataFrame(index=many.index, columns=many.index)

    for i, row in enumerate(many.index):
        for j in range(i, len(many.index)):
            col = many.index[j]
            c, p = the_func(many.loc[row], many.loc[col])
            cor.loc[row, col] = c
            cor.loc[col, row] = c
            pval.loc[row, col] = p
            pval.loc[col, row] = p

    return cor, pval


def standardise(data, axis=1):
    """
    Standardise the RMA data in the exon array
    Defined as median centred and divided by the IQR.
    :param data:
    :param axis: The axis to run this - default is rows (probes)
    :return:
    """
    other_axis = 0 if axis == 1 else 1
    return data.subtract(data.median(axis=axis), axis=other_axis).divide(stats.iqr(data, axis=axis), axis=other_axis)


if __name__ == '__main__':
    myc_probe = 3115514
    myc_gene = 'MYC'
    corr_method = 'pearson'
    cross_corr_threshold = 0.5
    within_corr_threshold = 0.75
    alpha = 0.05

    fn = os.path.join(DATA_DIR_NON_GIT, 'exon_array', 'GSE60982', 'raw', 'GSE60892_HuEx-ALL40-Upload-Transposed.txt.gz')
    rma_data = pd.read_csv(fn, sep='\t', comment='#', header=0, index_col=0)

    ann_fn = os.path.join(DATA_DIR_NON_GIT, 'exon_array', 'GSE60982', 'HuEx-1_0-st-v2.na36.hg19.probeset.csv.gz')
    ann = pd.read_csv(ann_fn, sep=',', comment='#', header=0, index_col=0, na_values='---')

    # ann_fn2 = os.path.join(DATA_DIR_NON_GIT, 'exon_array', 'ash_cpt_project', 'GSE60892_annotation.xlsx')
    # ann2 = pd.read_excel(ann_fn2)

    # annotation files
    # latest version from Affymetrix / Thermofisher
    # every probe in this annotation is defined as 'core'
    ann = ann.loc[rma_data.index]
    ann = ann.loc[~ann.gene_assignment.isnull()]
    ann.loc[:, 'gene_assignment'] = ann.loc[:, 'gene_assignment'].str.replace(r' /+ ', '|')
    the_symbols = ann.gene_assignment.str.split('|').apply(lambda x: x[1])

    # version used by Jennie originally
    # ann2 = ann2.loc[~ann2.probeset_id.isnull()]
    # ann2.loc[:, 'probeset_id'] = ann2.loc[:, 'probeset_id'].astype(int)
    # ann2.set_index("probeset_id", inplace=True)
    # ann2 = ann2.loc[rma_data.index]
    # ann2 = ann2.loc[~ann2.gene_assignment.isnull()]
    # ann2.loc[:, 'gene_assignment'] = ann2.loc[:, 'gene_assignment'].str.replace(r' /+ ', '|')
    # symb2 = ann2.gene_assignment.str.split('|').apply(lambda x: x[1] if len(x) > 1 else None)

    # get a standardised dataset
    the_data = standardise(rma_data)

    # get all MYC probes
    myc_probes = the_symbols.index[the_symbols == myc_gene]
    n_myc_probes_orig = len(myc_probes)

    # reduce down to those that correlate with one another
    # form the square pairwise corr matrix, remove rows that don't correlate
    myc_probe_corr, myc_probe_pval = pairwise_correlation(the_data.loc[myc_probes], method=corr_method)
    drop_idx = ((myc_probe_corr < within_corr_threshold) & (myc_probe_pval > alpha)).sum(axis=0) == (n_myc_probes_orig - 1)
    myc_probes = myc_probes[~drop_idx]

    print "Of the original %d MYC probes, we retained %d that correlate well: %s." % (
        n_myc_probes_orig,
        len(myc_probes),
        ', '.join(myc_probes.astype(str))
    )

    # find all probes correlated with them (individually, taking the union)
    probes_corr_with_myc = []

    for p in myc_probes:
        print "Computing correlation with MYC probe %s" % str(p)
        base = the_data.loc[p]
        cor, pval = one_vs_many_correlation(the_data.loc[p], the_data, method=corr_method)
        keep_probes = cor.index[(cor.abs() > cross_corr_threshold) & (pval < alpha)]
        probes_corr_with_myc.append(keep_probes)

    #  out of interest, what is the overlap between these? (presumably quite high?)
    vs, vc = setops.venn_from_arrays(*probes_corr_with_myc)

    # union of probes
    probes_corr_with_myc_union = setops.reduce_union(*probes_corr_with_myc)

    print "After comparing all data against each MYC probe, we are left with %d correlated probes" % len(probes_corr_with_myc_union)

    # aggregate by gene (only within the pre-selected probes)
    genes_corr_with_myc = the_symbols.loc[probes_corr_with_myc_union]
    dat_corr_with_myc_aggr = the_data.loc[probes_corr_with_myc_union].groupby(genes_corr_with_myc, axis=0).mean()

    # re-run correlation analysis at the gene level
    base = dat_corr_with_myc_aggr.loc[myc_gene]
    cor_gene, pval_gene = one_vs_many_correlation(base, dat_corr_with_myc_aggr, method=corr_method)

    # reduce to significant and relevant
    keep_genes = cor_gene.index[(cor_gene.abs() > cross_corr_threshold) & (pval_gene < alpha)]

    # remove MYC itself when reporting
    print "Having aggregated these by gene, %d are correlated with %s" % (len(keep_genes.drop(myc_gene)), myc_gene)

    # cluster using this representation of the data
    # force the order of the columns to match the correlation with MYC
    keep_genes = cor_gene.loc[keep_genes].sort_values(ascending=False).index
    cg = clustering.plot_clustermap(
        dat_corr_with_myc_aggr.loc[keep_genes],
        cmap='RdBu_r',
        metric='euclidean',
        method='ward',
        row_cluster=False,
        vmin=-4.5, vmax=4.5
    )
    cg.gs.update(bottom=0.1)



    # aggregate all genes
    all_aggr_by_gene = rma_data.groupby(the_symbols, axis=0).mean()

    # Pearson correlation by gene against the cymc probe
    base = rma_data.loc[myc_probe]
    cor, pval = one_vs_many_correlation(base, rma_data, method=corr_method)

    # select probes exhibiting correlation
    probes = cor.index[(cor.abs() > 0.5) & (pval < 0.05)]
    genes = the_symbols.loc[probes].dropna()
    print "%d probes selected that correlate with MYC, corresponding to %d unique genes" % (probes.size, genes.unique().size)

    # drop probes that have no gene assignment
    probes_dropped = probes.difference(genes.index)
    print "Dropped %d probes as they had no gene assignment" % probes_dropped.size
    probes = genes.index

    # extract this block of data
    sub_rma_data = rma_data.loc[probes]

    # standardise
    sub_rma_data_ctr = standardise(sub_rma_data)

    # aggregate by gene
    aggr_rma_data_ctr_by_gene = sub_rma_data_ctr.groupby(genes, axis=0).mean()

    # re-run correlation with Myc
    cor_gene_ctr, pval_gene_ctr = one_vs_many_correlation(
        aggr_rma_data_ctr_by_gene.loc[myc_gene],
        aggr_rma_data_ctr_by_gene,
        method=corr_method
    )

    # go back to the genes, taking the average probe aggregation (for some value of average)
    glist = genes.unique()
    plist = the_symbols[the_symbols.isin(glist)].index

    aggr_by_gene = rma_data.loc[plist].groupby(the_symbols[plist], axis=0).mean()

    # re-run correlation with MYC
    cor_gene, pval_gene = one_vs_many_correlation(aggr_by_gene.loc[myc_gene], aggr_by_gene, method=corr_method)

    ### alternative approach: just take the gene signature and replot the heatmap
    # NB: some typos in here (wow - done by hand?!), now corrected


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
        'BAG3', # 'NBAG3',
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
        'CLCF1', # 'CCLF1',
        'TFPI2', # 'TFP12',
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
        'PIGA', # 'P1GA',
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
        'TMEM130', # 'TMEM',
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