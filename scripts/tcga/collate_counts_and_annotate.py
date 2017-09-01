import os
import pandas as pd
from settings import DATA_DIR_NON_GIT
import re
import json
import references

## RNASEQ
base_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour')
indir = os.path.join(base_dir, 'htseq-count_fpkm')

dirlist = os.listdir(indir)

df = references.conversion_table()
df = df.loc[~df.loc[:, 'Ensembl Gene ID'].duplicated()].set_index('Ensembl Gene ID')

dat = None  # this will hold all our data
first = True

for d in dirlist:
    if re.match(r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}', d):
        the_dir = os.path.join(indir, d)
        the_file = os.listdir(the_dir)[0]  # assuming only one file per directory
        # load data and convert to a Series
        the_col = pd.read_csv(os.path.join(the_dir, the_file), sep='\t', index_col=0, header=None).iloc[:, 0]
        # only keep genes
        the_col = the_col.loc[the_col.index.str.contains('ENSG0')]
        # remove version numbers from the row names
        the_col.index = the_col.index.str.replace(r'\.[0-9]*', '')
        # rename
        col_name = the_file
        the_col.name = col_name

        # if this is the first file loaded, use it as a template
        if first:
            dat = pd.DataFrame(index=the_col.index)
            dat.loc[:, col_name] = the_col # the directory is the case ID
            first = False
        else:
            dat.loc[:, col_name] = the_col  # the directory is the case ID

# now we can use the annotation provided by Brennan et al. and the manifest to add subgroup information
# this is the old WHO subgrouping
# we do this via the
fn = os.path.join(
    base_dir,
    'brennan_s7.csv'
)
meta = pd.read_csv(fn, index_col=0, header=0)

fn = os.path.join(
    indir,
    'cases.json'
)
with open(fn, 'rb') as f:
    cases = json.load(f)
cases_mani = pd.Series(
    [t['submitter_id'] for t in cases],
    index=[t['case_id'] for t in cases]
)

fn = os.path.join(
    indir,
    'files.json'
)
with open(fn, 'rb') as f:
    files = json.load(f)
files_mani = pd.Series(
    [t['cases'][0]['case_id'] for t in files],
    index=[t['file_name'] for t in files],
)

# link filename -> case ID
aa = files_mani.loc[dat.columns[:-1]]

# I see duplicates here (2) - WHY? Probably the same tumour split into multiple samples?

# report them - show all cases
print "Duplicates detected: \n%s" % aa.loc[pd.Index(aa.values).duplicated(keep=False)].to_string()

# pick one in each case
# duplicated() doesn't mark the first entry by default
aa = aa[~pd.Index(aa.values).duplicated()]

dat2 = dat.loc[:, aa.index.intersection(dat.columns)]
cols = aa.loc[dat2.columns].values

bb = cases_mani.loc[cols]
# these should be unique...
if pd.Index(bb.values).duplicated().any():
    raise AttributeError("Some case IDs are duplicated in the final dataset")
dat2.columns = bb.values

# finally, only keep those cases that are also in the Brennan table
not_in_meta = pd.Index(bb.values).difference(meta.index)

if len(not_in_meta):
    this = bb.loc[pd.Index(bb.values).isin(not_in_meta)]
    print "Cases not in meta: \n%s" % this.to_string()

meta = meta.loc[meta.index.intersection(bb.values)]
dat2 = dat2.loc[:, meta.index]

# add gene symbols
gs = references.ensembl_to_gene_symbol(dat2.index)
dat2.loc[:, 'Approved Symbol'] = gs

# export
meta.to_csv(os.path.join(indir, 'sources.csv'))
dat2.to_csv(os.path.join(indir, 'counts.csv'))