import os
import pandas as pd
from settings import DATA_DIR_NON_GIT


basedir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq')
dirname = 'wtchg_p170389'
indir = os.path.join(basedir, dirname)

flist = []
subdirs = []

for fn in os.listdir(indir):
    ff = os.path.join(indir, fn)
    if os.path.isdir(ff):
        src_ff = os.path.join(ff, 'sources.csv')
        if os.path.isfile(src_ff):
            flist.append(src_ff)
            subdirs.append(ff)

print "Found %d lanes." % len(subdirs)

reads = None
for fn in flist:
    this_meta = pd.read_csv(fn, header=0, index_col=None)
    try:
        this_meta.set_index('sample', inplace=True)
    except Exception:
        this_meta.set_index('filename', inplace=True)
    if reads is None:
        reads = this_meta.read_count
    else:
        reads += this_meta.read_count

print reads.to_csv()