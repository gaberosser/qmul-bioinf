from scripts.agdex_mouse_human_mb_microarray import load_sb_mouse_data
import pandas as pd
import numpy as np
from settings import DATA_DIR
import os


if __name__ == '__main__':
    # load list of mouse housekeeping genes
    fn = os.path.join(DATA_DIR, 'mouse_housekeeping/nature05453_supptable1.csv')
    hkg_all = pd.read_csv(fn)
    hkg = hkg_all.loc[:, 'Near Ubiquitious']
    hkg = hkg[hkg.notnull()].values

    arr, chd7 = load_sb_mouse_data.load_array_data()
    # drop the mRNA column
    arr = arr.loc[:, arr.columns.difference(['mRNA Accession'])]
    # median of probes in each set
    arr = arr.groupby('Gene Symbol').median()


    # extract putative HKG from mouse array
    hkg_arr = arr.loc[arr.index.intersection(hkg)]
    # sort by descending median expression level
    hkg_arr = hkg_arr.loc[hkg_arr.median(axis=1).sort_values().index[::-1]]

    # play with this - a weighted avg might be nicer than choosing one?
    arr_norm = arr.divide(hkg_arr.iloc[0], axis=1)