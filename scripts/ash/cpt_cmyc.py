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

    fn = os.path.join(DATA_DIR_NON_GIT, 'exon_array', 'GSE60982', 'raw', 'GSE60892_HuEx-ALL40-Upload-Transposed.txt.gz')
    rma_data = pd.read_csv(fn, sep='\t', comment='#', header=0, index_col=0)

    ann_fn = os.path.join(DATA_DIR_NON_GIT, 'exon_array', 'GSE60982', 'HuEx-1_0-st-v2.na36.hg19.probeset.csv.gz')
    ann = pd.read_csv(ann_fn, sep=',', comment='#', header=0, index_col=0, na_values='---')

    ann_fn2 = os.path.join(DATA_DIR_NON_GIT, 'exon_array', 'ash_cpt_project', 'GSE60892_annotation.xlsx')
    ann2 = pd.read_excel(ann_fn2)

    # link probes to genes
    ann = ann.loc[rma_data.index]
    ann = ann.loc[~ann.gene_assignment.isnull()]
    ann.loc[:, 'gene_assignment'] = ann.loc[:, 'gene_assignment'].str.replace(r' /+ ', '|')
    symb = ann.gene_assignment.str.split('|').apply(lambda x: x[1])

    ann2 = ann2.loc[~ann2.probeset_id.isnull()]
    ann2.loc[:, 'probeset_id'] = ann2.loc[:, 'probeset_id'].astype(int)
    ann2.set_index("probeset_id", inplace=True)
    ann2 = ann2.loc[rma_data.index]
    ann2 = ann2.loc[~ann2.gene_assignment.isnull()]
    ann2.loc[:, 'gene_assignment'] = ann2.loc[:, 'gene_assignment'].str.replace(r' /+ ', '|')
    symb2 = ann2.gene_assignment.str.split('|').apply(lambda x: x[1] if len(x) > 1 else None)

    the_symbols = symb2

    # aggregate all genes
    aggr_by_gene = rma_data.groupby(the_symbols, axis=0).mean()

    # Pearson correlation by gene against the cymc probe
    base = rma_data.loc[cmyc_probe]
    pearson_cor, pearson_pval = one_vs_many_correlation(base, rma_data)

    all_others = rma_data.loc[rma_data.index != cmyc_probe]

    pearson_cor = pd.Series(index=all_others.index)
    pearson_pval = pd.Series(index=all_others.index)
    for pid, vals in all_others.iterrows():
        pearson_cor[pid], pearson_pval[pid] = stats.pearsonr(base, vals)

    # select probes exhibiting correlation
    probes = pearson_cor.index[(pearson_cor.abs() > 0.5) & (pearson_pval < 0.05)]
    genes = the_symbols.loc[probes].dropna()
    # genes2 = symb2.loc[probes].dropna()

    # go back to the genes, taking the median probe aggregation
    glist = genes.unique()
    plist = the_symbols[the_symbols.isin(glist)].index

    aggr_by_gene = rma_data.loc[plist].groupby(the_symbols[plist], axis=0).mean()

    # re-run correlation with MYC
    pearson_cor_gene, pearson_pval_gene = one_vs_many_correlation(aggr_by_gene.loc[cmyc_gene], aggr_by_gene)

