import os
import csv
import pandas as pd


def known_genes():
    in_file = 'reference_data/genenames.org.tsv'
    with open(in_file, 'rb') as f:
        df = pd.read_csv(f, delimiter='\t')
    return df


def gene_symbol_to_entrez(g):
    cat = known_genes()
    if isinstance(g, str) or not hasattr(g, '__iter__'):
        g = [g]
    res = [cat.loc[cat.loc[:, 'Approved Symbol'] == x, 'Entrez Gene ID'].values for x in g]
    for i in range(len(res)):
        if len(res[i]) == 0:
            res[i] = None
        elif len(res[i]) == 1:
            res[i] = res[i][0]
    return res
