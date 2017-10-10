import os
import pandas as pd
from settings import DATA_DIR_NON_GIT


basedir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq')
dirname = 'wtchg_p170446'
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
dupe = None

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

    # look for duplication levels
    idx = this_meta.columns.str.contains('dedupe')
    if idx.any():
        if dupe is None:
            dupe = this_meta.loc[:, idx]
        else:
            dupe += this_meta.loc[:, idx]

print reads.to_csv()

if dupe is not None:
    print "Duplicates"
    print (dupe / float(len(subdirs))).to_csv()