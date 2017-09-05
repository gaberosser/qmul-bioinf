import os
import csv
import pandas as pd
from settings import DATA_DIR


def known_genes(tax_id=9606, index_field=None):
    """
    Get dataframe of known genes for the supplied taxonomy ID (default: human)
    Other options: mouse (10090)
    :param index_field: If supplied, reindex the dataframe by this field. Examples are 'Symbol', 'GeneID'
    """
    if tax_id == 9606:
        fn = os.path.join(DATA_DIR, 'ncbi_gene', 'Homo_sapiens.gene_info.gz')
    elif tax_id == 10090:
        fn = os.path.join(DATA_DIR, 'ncbi_gene', 'Mus_musculus.gene_info.gz')
    else:
        raise ValueError("Unrecognised tax_id %d.", tax_id)
    df = pd.read_csv(fn, sep='\t')
    if index_field is not None:
        df = df.set_index(index_field)
    return df


def conversion_table(type='protein_coding'):
    if type == 'protein_coding':
        in_file = os.path.join(DATA_DIR, 'genenames', 'protein_coding', 'genenames.org.2017.01.tsv')
    elif type == 'all':
        in_file = os.path.join(DATA_DIR, 'genenames', 'all', 'genenames.org.2017.09.tsv')
    else:
        raise ValueError("Unsupported type option '%s'" % type)
    # in_file = os.path.join(DATA_DIR, 'genenames', 'genenames.org.tsv')
    df = pd.read_csv(in_file, delimiter='\t')
    return df


def _translate(cat, x, to_field, from_field):
    return cat.set_index(from_field).loc[x, to_field]

def translate(x, to_field, from_field):
    cat = conversion_table()
    return cat.set_index(from_field).loc[x, to_field]


def gene_symbol_to_entrez(g):
    cat = conversion_table(type='all')
    return _translate(cat, g, 'Entrez Gene ID', 'Approved Symbol')


def entrez_to_gene_symbol(e):
    cat = conversion_table(type='all')
    return _translate(cat, e, 'Approved Symbol', 'Entrez Gene ID')


def gene_symbol_to_ensembl(g):
    cat = conversion_table(type='all')
    return _translate(cat, g, 'Ensembl Gene ID', 'Approved Symbol')


def ensembl_to_gene_symbol(e):
    cat = conversion_table(type='all')
    return _translate(cat, e, 'Approved Symbol', 'Ensembl Gene ID')


def ensembl_to_name(e):
    cat = conversion_table(type='all')
    return _translate(cat, e, 'Approved Name', 'Ensembl Gene ID')