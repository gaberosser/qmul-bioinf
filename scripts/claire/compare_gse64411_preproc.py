import pandas as pd
import os
from settings import DATA_DIR_NON_GIT, GIT_LFS_DATA_DIR
from load_data import rnaseq_data
import numpy as np
from scipy import stats


if __name__ == '__main__':
    # create a lookup between RefSeq and Ensembl
    conv_fn = os.path.join(GIT_LFS_DATA_DIR, 'biomart', 'mm10', 'mm10_refseq_ensembl.txt')
    conv = pd.read_csv(conv_fn, header=0).dropna()
    conv.set_index(conv.columns[-1], inplace=True)
    conv = conv.loc[:, 'Gene stable ID']

    # load my version
    obj = rnaseq_data.gse64411(annotate_by='Ensembl Gene ID')
    my_idx = obj.data.index[obj.data.index.str.contains('ENS')]

    # load preprocessed datasets
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE64411', 'preprocessed')
    fns = ['%s_tagcounts.txt.gz' % t for t in (
        'GSM1570771_I_2', 'GSM1570773_I_5', 'GSM1570775_I_7', 'GSM1570777_I_9'
    )]
    accessions = [t[:10] for t in fns]
    preproc = {}
    the_idx = None
    for fn in fns:
        acc = fn[:10]
        the_dat = pd.read_csv(os.path.join(indir, fn), sep='\t', header=0, index_col=0)
        # convert to ENS
        idx = the_dat.index.str.replace(r'\.[0-9]', '')
        the_dat.index = idx
        idx_lookup = conv.loc[idx].dropna()
        to_keep = idx_lookup.index
        the_dat = the_dat.loc[to_keep]
        the_dat.index = to_keep
        # sum over any transcripts that correspond to the same gene
        the_dat = the_dat.groupby(idx_lookup).sum()
        the_dat.index.name = acc

        if the_idx is None:
            the_idx = the_dat.index.intersection(my_idx)

        the_dat = the_dat.loc[the_idx]

        preproc[acc] = the_dat

    pre = pd.DataFrame(index=the_idx, columns=preproc.keys())
    for k, v in preproc.items():
        pre.loc[the_idx, k] = v.values.flatten()

    # compare
    gsm_sample = obj.meta.set_index('accession')
    gsm_sample = gsm_sample.loc[:, 'sample']

    for acc in accessions:
        a = pre.loc[the_idx, acc]
        b = obj.data.loc[the_idx, gsm_sample.loc[acc]]
        al = np.log2(a.values + 1)
        bl = np.log2(b.values + 1)

        print acc
        print stats.linregress(al, bl)