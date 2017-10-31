import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from settings import GIT_LFS_DATA_DIR
import os


if __name__ == "__main__":
    fn = os.path.join(GIT_LFS_DATA_DIR, 'ash_cpt', 'cpt.csv')
    # load data
    a = pd.read_csv(fn, header=0, index_col=0).dropna(how='all')
    arr = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            arr[i, j] = ((a.loc[:, 'C-Myc'] == (2 - i)) & (a.loc[:, 'P53 '] == j)).sum()

    sns.heatmap(arr, cmap='RdBu_r')
    plt.xlabel('P53 focal staining', fontsize=14)
    plt.ylabel('C-Myc focal staining', fontsize=14)
    ticklabels = [
        'Absent',
        'Diffuse',
        'Focal'
    ]
    ax = plt.gca()
    ax.set_xticklabels(ticklabels, fontsize=14)
    ax.set_yticklabels(ticklabels[::-1], fontsize=14)

    fig = ax.figure
    fig.savefig('cpt_staining.tiff', dpi=200)
    fig.savefig('cpt_staining.jpg', dpi=200)