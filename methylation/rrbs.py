import pandas as pd
import os
import csv


def load_bismark_mbias_report(fn):
    res = {}
    with open(fn, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        this_res = None
        this_k = None
        for row in c:
            if len(row) < 5:
                # start a new section, saving the old one if applicable
                if this_res is not None:
                    res[this_k] = this_res
                    this_res = None
                    this_k = None
                if len(row) and row[0][0] != '=':
                    this_k = row[0]
            else:
                if this_res is None:
                    this_res = pd.DataFrame(columns=row[1:])
                else:
                    this_res.loc[int(row[0])] = [int(row[1]), int(row[2]), float(row[3]), int(row[4])]
    return res
