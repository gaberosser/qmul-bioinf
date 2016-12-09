import pandas as pd
import os
import references
from settings import DATA_DIR

INDEX_FIELDS = (
    'Approved Symbol',
    'Entrez Gene ID',
    'RefSeq IDs',
    'Ensembl Gene ID'
)


def gse83696(index_by='Ensembl Gene ID'):
    """
    Data are initially indexed by Ensembl ID. Coversion is carried out using HGNC data, if required.
    Index field options are: Approved Symbol, Entrez Gene ID, RefSeq IDs
    :param index_by:
    :return:
    """
    # TODO: convert this into a generic loader for htseq-count outputs
    indir = os.path.join(DATA_DIR, 'rnaseq_GSE83696', 'htseq-count')
    samples = [
        ('XZ1', 'XZ-1.count'),
        ('XZ2', 'XZ-2.count'),
        ('XZ3', 'XZ-3.count'),
        ('XZ4', 'XZ-4.count'),
        ('XZ5', 'XZ-5.count'),
        ('XZ6', 'XZ-6.count'),
        ('XZ7', 'XZ-7.count'),
        ('XZ8', 'XZ-8.count'),
    ]
    df = pd.DataFrame()
    for sn, fn in samples:
        t = pd.read_csv(os.path.join(indir, fn), sep='\t', index_col=0, header=None).iloc[:, 0]
        df.loc[:, sn] = t

    if index_by is not None and index_by != 'Ensembl Gene ID':
        new_idx = references.translate(df.index, to_field=index_by, from_field='Ensembl Gene ID')
        new_idx.dropna(inplace=True)
        df = df.loc[new_idx.index]
        df.index = new_idx.values
    return df