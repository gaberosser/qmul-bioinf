import os
import pandas as pd
import numpy as np
from microarray import illumina, process
from settings import DATA_DIR


def load_from_r_processed(infile, sample_names, index_field):
    meta_cols = ['probeset_id', 'accession_id', 'gene_symbol', 'description', 'ensembl_id', 'entrez_id']
    cols = meta_cols + sample_names

    arr_data = pd.read_csv(infile, sep='\t', header=None, names=cols, skiprows=1)
    # drop other meta fields
    arr_data = arr_data.loc[:, [index_field] + sample_names]
    arr_data = process.aggregate_by_probe_set(arr_data, groupby=index_field)
    return arr_data


def load_annotated_microarray_gse54650(index_field='entrez_id'):
    infile = os.path.join(DATA_DIR, 'microarray_GSE54650', 'data.ann.txt')

    sample_names = [
        'Cer_CT18',
        'Cer_CT20',
        'Cer_CT22',
        'Cer_CT24',
        'Cer_CT26',
        'Cer_CT28',
        'Cer_CT30',
        'Cer_CT32',
        'Cer_CT34',
        'Cer_CT36',
        'Cer_CT38',
        'Cer_CT40',
        'Cer_CT42',
        'Cer_CT44',
        'Cer_CT46',
        'Cer_CT48',
        'Cer_CT50',
        'Cer_CT52',
        'Cer_CT54',
        'Cer_CT56',
        'Cer_CT58',
        'Cer_CT60',
        'Cer_CT62',
        'Cer_CT64',
    ]
    return load_from_r_processed(infile, sample_names, index_field=index_field)


def load_annotated_microarray_sb_data(index_field='entrez_id'):

    infile = os.path.join(DATA_DIR, 'sleeping_beauty_mouse_screen', 'data.ann.txt')
    sample_names = [
        "Wu050",
        "Wu053",
        "Wu054",
        "Wu051",
        "Wu052",
        "Wu055",
        "Wu056",
        "Wu057"
    ]
    arr_data = load_from_r_processed(infile, sample_names, index_field=index_field)
    # CHD7 status
    chd7 = pd.Series(data=[True] * 3 + [False] * 5, index=sample_names)
    return arr_data, chd7


def load_annotated_microarray_gse37382(index_field='entrez_id'):
    infile = os.path.join(DATA_DIR, 'microarray_GSE37382', 'data.ann.txt.gz')
    meta_fn = os.path.join(DATA_DIR, 'microarray_GSE37382', 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=0, sep=',')
    sample_names = list(meta.index)

    arr = load_from_r_processed(infile, sample_names, index_field=index_field)
    return arr, meta