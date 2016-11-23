import pandas as pd
import os
from settings import DATA_DIR


def gse83696(index_by=None):
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

    if index_by is None:
        return df
    else:
        pass