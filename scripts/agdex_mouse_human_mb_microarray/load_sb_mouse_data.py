import pandas as pd
import os
import numpy as np
from settings import DATA_DIR


def load_array_data():
    fn = os.path.join(DATA_DIR, 'sleeping_beauty_mouse_screen', 'Dubuc_BMi1_Tg Expression profile.csv')
    arr = pd.read_csv(fn, header=1, index_col=0)
    # strip whitespace from gene symbol
    arr.loc[:, 'Gene Symbol'] = arr.loc[:, 'Gene Symbol'].apply(lambda x: x.strip())
    arr.loc[:, 'mRNA Accession'] = arr.loc[:, 'mRNA Accession'].apply(lambda x: x.strip())
    # replace --- with null
    arr.replace(to_replace='---', value=np.nan, inplace=True)

    # load CHD7 insertion marker
    chd7_lab = pd.read_csv(fn, header=None, nrows=1, usecols=range(1, 9))
    chd7_lab.columns = arr.columns[1:9]

    return arr, chd7_lab


if __name__ == '__main__':
    pass