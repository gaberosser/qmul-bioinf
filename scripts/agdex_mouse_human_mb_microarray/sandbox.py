import pandas as pd
import numpy as np
from settings import GIT_LFS_DATA_DIR
import os


if __name__ == '__main__':
    infile = os.path.join(GIT_LFS_DATA_DIR, 'sleeping_beauty_mouse_screen', 'Dubuc_BMi1_Tg Expression profile.csv')
    arr = pd.read_csv(infile, skiprows=1, header=0, index_col=0)
    arr = arr.iloc[:, :-2].astype(float)
    #

    