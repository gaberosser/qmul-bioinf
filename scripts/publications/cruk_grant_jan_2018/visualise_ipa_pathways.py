import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
import sys


if __name__ == '__main__':
    ## TODO: put this somewhere more portable!
    home = os.path.expanduser("~")

    # pathways by DE
    fn = os.path.join(home, 'Dropbox', 'research', 'qmul', 'data', 'hgic_project', 'ipa_de_all.xlsx')
    df = pd.read_excel(fn)

    # separate into patients
    pids = df.columns[::2].str.replace('GBM', '')
    res = {}
    for i, pid in enumerate(pids):
        sub_df = df.iloc[1:, (2 * i):(2 * i + 2)]
        sub_df.columns = ['pathway', '-log_pval']
        sub_df.set_index('pathway', inplace=True)
        res[pid] = sub_df.dropna()

    # top 30
    p30 = {}
    p30_all = set()
    for pid, df in res.items():
        p30[pid] = df.index[:30]
        p30_all.update(p30[pid])

    p30_all = sorted(p30_all)

    df_in30 = pd.DataFrame(False, index=p30_all, columns=pids)
    df_30 = pd.DataFrame(index=p30_all, columns=pids, dtype=float)

    for pid, d in res.items():
        df_30.loc[p30[pid], pid] = res[pid].loc[p30[pid]].values
        df_in30.loc[p30[pid], pid] = True

    # plot these with a sort of clustering output
    
