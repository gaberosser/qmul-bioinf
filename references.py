import os
import csv
import pandas as pd
from settings import DATA_DIR


def known_genes():
    in_file = os.path.join(DATA_DIR, 'genenames', 'genenames.org.tsv')
    with open(in_file, 'rb') as f:
        df = pd.read_csv(f, delimiter='\t')
    return df


def _translate(cat, x, to_field, from_field):
    return cat.set_index(from_field).loc[x, to_field]


def gene_symbol_to_entrez(g):
    cat = known_genes()
    return _translate(cat, g, 'Entrez Gene ID', 'Approved Symbol')


def entrez_to_gene_symbol(e):
    cat = known_genes()
    return _translate(cat, e, 'Approved Symbol', 'Entrez Gene ID')


def ensembl_to_gene_symbol(e):
    cat = known_genes()
    return _translate(cat, e, 'Approved Symbol', 'Ensembl Gene ID')