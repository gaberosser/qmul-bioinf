import pandas as pd
import os
import numpy as np
from settings import DATA_DIR
from microarray import process


def load_annotated_array_data():
    fn = os.path.join(DATA_DIR, 'sleeping_beauty_mouse_screen', 'Dubuc_BMi1_Tg Expression profile.csv')
    arr = pd.read_csv(fn, header=1, index_col=0)

    # remove mRNA accession labels
    arr = arr.loc[:, arr.columns.difference(['mRNA Accession'])]

    # strip whitespace from gene symbol
    arr.loc[:, 'Gene Symbol'] = arr.loc[:, 'Gene Symbol'].apply(lambda x: x.strip())

    # arr.loc[:, 'mRNA Accession'] = arr.loc[:, 'mRNA Accession'].apply(lambda x: x.strip())

    # replace --- with null
    arr.replace(to_replace='---', value=np.nan, inplace=True)

    # aggregate over redundant probe groups
    arr = process.aggregate_by_probe_set(arr, groupby='Gene Symbol', method='median')

    # shorten column names
    arr.columns = [t[:5] for t in arr.columns]

    # load CHD7 insertion marker
    chd7_lab = pd.read_csv(fn, header=None, nrows=1, usecols=range(1, 9))
    chd7_lab.columns = arr.columns

    return arr, chd7_lab


if __name__ == '__main__':
    pass