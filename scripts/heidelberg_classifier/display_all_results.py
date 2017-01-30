import pandas as pd
import zipfile
import os
from settings import DATA_DIR_NON_GIT
import re


if __name__ == "__main__":
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation')
    for t in os.walk(indir):
        if 'heidelberg_classifier' in t[1]:
            this_dir = os.path.join(t[0], 'heidelberg_classifier')
            for fn in os.listdir(this_dir):
                if fn[-4:] == '.zip':
                    ff = os.path.join(this_dir, fn)
                    print ff
                    with zipfile.ZipFile(ff, mode='r') as z:
                        flist = [x for x in z.namelist() if re.search(r'scores_cal.csv$', x)]
                        if len(flist) != 1:
                            raise Exception("Could not find a unique .scores_cal.csv member")
                        with z.open(flist[0]) as f:
                            dat = pd.read_csv(f, header=0, index_col=0)
                            dat.sort_values(by=dat.columns[0], ascending=False, inplace=True)
                            # start printing
                            for x, y in dat.iterrows():
                                if y[0] > 0.04:
                                    print '{file: <10}\t\t\t{key: <12}\t{val:.1f}'.format(
                                        file=fn.replace('.zip', ''), key=x, val=y[0] * 100.)
                    print "***"


