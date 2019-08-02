import pandas as pd
import zipfile
import os
import numpy as np
from settings import DATA_DIR
import re


def read_cal_scores(filename, min_score=None):
    """
    :param filename: Path to a zipped file containing the results of a classifier run
    :param min_score: Optionally specify a minimum score, below which matches will not be reported
    :return: A table of calibrated scores
    """
    dat = None
    with zipfile.ZipFile(filename, mode='r') as z:
        flist = [x for x in z.namelist() if re.search(r'scores_cal.csv$', x)]
        if len(flist) != 1:
            raise Exception("Could not find a unique .scores_cal.csv member")
        with z.open(flist[0]) as f:
            dat = pd.read_csv(f, header=0, index_col=0)
            dat.sort_values(by=dat.columns[0], ascending=False, inplace=True)
            if min_score is not None:
                dat = dat[dat.iloc[:, 0] > min_score]
    return dat


def get_top_n_results(indir, n_result=4, score_threshold=0.6):
    """
    Recurse through all zip files in the directory.
    Load the calibrated scores and generate a summary report of the top <=n matches
    :param indir:
    :param n_result:
    :param score_threshold: Minimum share of the top results required to declare a 'match'. Must be >0.5 to make sense.
    :return:
    """
    index = [
        'match',
    ]
    for i in range(n_result):
        index += ['family_%d' % (i + 1), 'family_%d_score' % (i + 1)]

    flist = [t for t in os.listdir(indir) if t[-4:].lower() == '.zip']
    data = []

    for fn in flist:
        name = re.sub(re.compile(r'\.zip', flags=re.IGNORECASE), '', fn)
        ff = os.path.join(indir, fn)
        this_res = read_cal_scores(ff).iloc[:n_result]
        midx = np.where(this_res.divide(this_res.sum(axis=0)) >= score_threshold)[0]

        match = 'None'
        if len(midx) == 1:
            match = this_res.index[midx[0]]

        remainder = list(reduce(lambda x, y: x+y, this_res.iloc[:, 0].iteritems()))
        s = pd.Series(data=[match] + remainder, index=index, name=name)
        data.append(s)

    return pd.DataFrame(data=data)


if __name__ == "__main__":
    indir = os.path.join(DATA_DIR, 'methylation')
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


