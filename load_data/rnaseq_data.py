import pandas as pd
import os
from settings import DATA_DIR


def gse83696():
    # TODO: convert this into a generic loader for htseq-count outputs
    indir = os.path.join(DATA_DIR, 'rnaseq_GSE83696', 'htseq-count')
    samples = [
        ('XZ1', 'xz1.counts'),
        ('XZ2', 'xz2.counts'),
    ]
    df = pd.DataFrame()
    for sn, fn in samples:
        t = pd.read_csv(os.path.join(indir, fn), sep='\t', index_col=0, header=None).iloc[:, 0]
        df.loc[:, sn] = t

