import os
import pandas as pd
from settings import OUTPUT_DIR
import re


if __name__ == '__main__':
    indir = os.path.join(OUTPUT_DIR, 'heidelberg_results.0')

    n = 4
    match_threshold = 0.3
    nm = 'No match'

    cols = ['Study ID', 'Batch', 'Result']
    for i in range(1, n + 1):
        cols += ['Reference %d' % i, '%% Match %d' % i]
    res = pd.DataFrame(columns=cols)

    i = 0

    for batch in os.listdir(indir):
        subdir = os.path.join(indir, batch)
        print batch
        for fn in os.listdir(subdir):
            if fn[0] == '.':
                # skip hidden file
                continue
            if 'calibrated_scores.csv' in fn:
                sname = fn.replace('.calibrated_scores.csv', '')
                print sname
                ff = os.path.join(subdir, fn)
                scores = pd.read_csv(ff, header=0)
                if scores.Score.iloc[0] >= match_threshold:
                    match = scores.iloc[0, 0]
                else:
                    match = nm
                # convert scores to pct
                scores.Score *= 100.

                this_res = [sname, batch, match] + scores.iloc[:4, [0, 2]].values.flatten().tolist()
                res.loc[i] = this_res
                i += 1

    # this conveniently fills in over half of the hGIC IDs. The NH-IDs still need to be converted, though.
    res.loc[:, 'hGIC ID'] = res.iloc[:, 0].apply(lambda x: re.sub(r'.*(DURA|GBM)(?P<id>[0-9]{3}).*', 'GBM\g<id>', x))
    print res.to_csv(index=False)