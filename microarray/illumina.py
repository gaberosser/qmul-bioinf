import re
import os
import pandas as pd
import numpy as np
import csv


def load_full_probeset_definitions_txt(filename):
    with open(filename, 'rb') as f:
        # the first few lines do not contain data; skip until we find the header
        i = 0
        while True:
            line = f.readline()
            if re.match(r'ID', line):
                break
            i += 1
    res = pd.read_csv(filename, sep='\t', header=i, index_col=0)
    res.replace('---', np.nan, inplace=True)
    return res


def load_full_probeset_definitions_soft(filename):
    with open(filename, 'rb') as f:
        # the first few lines do not contain data; skip until we find the header
        i = 0
        while True:
            line = f.readline()
            if re.match(r'ID', line):
                break
            i += 1
    res = pd.read_csv(
        filename,
        sep='\t',
        header=0,
        index_col=0,
        comment='!',
        skip_blank_lines=True,
        skiprows=i,
        low_memory=False
    )
    res.replace('---', np.nan, inplace=True)
    return res


def load_full_probeset_definitions(filename, format='txt'):
    """
    Load the probe set definitions from a text file, as available from the GEO site.
    """
    if format.lower() == 'txt':
        return load_full_probeset_definitions_txt(filename)
    if format.lower() == 'soft':
        return load_full_probeset_definitions_soft(filename)


def load_data(directory, sample_names, extension=None, columns=('ID_REF', 'VALUE')):

    res = pd.DataFrame()

    for fn in sample_names:
        ff = os.path.join(directory, fn)
        ff += extension or ''
        this_df = pd.read_csv(
            ff,
            sep='\t',
            header=0,
            index_col=0,
            skip_blank_lines=True,
            usecols=columns,
            squeeze=True
        )
        res[fn] = this_df

    return res
