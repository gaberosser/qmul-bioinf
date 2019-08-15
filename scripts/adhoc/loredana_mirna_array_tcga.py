from rnaseq import tcga
import pandas as pd
import os
from utils.output import unique_output_dir


if __name__ == "__main__":
    path_to_manifest = '/home/gabriel/Downloads/gdc_manifest.2017-06-30T15-30-57.570863.tsv'
    outdir = unique_output_dir("nih_legacy")
    tcga.download_from_manifest(path_to_manifest, outdir=outdir, legacy=True)

    flist = os.listdir(outdir)
    fin = os.path.join(outdir, flist[0])
    dat = pd.read_csv(fin, sep='\t', index_col=0, skiprows=1, header=0)
    res = pd.DataFrame(index=dat.index, columns=flist)

    for fn in flist:
        fin = os.path.join(outdir, fn)
        this = pd.read_csv(fin, sep='\t', index_col=0, skiprows=1, header=0)
        res.loc[:, fn] = this.iloc[:, 0]
        if res.loc[:, fn].isnull().any():
            print "Some were null - %s" % fn
            print this

    # add column for miRNA name
    annot_fn = '/home/gabriel/data/microarray/miRNA/tcga/annotation/h-mirna_8x15k_annotations.txt'
    meta = pd.read_csv(annot_fn, set='\t', index_col=1)
    meta = meta.loc[~pd.isnull(meta.index)]
