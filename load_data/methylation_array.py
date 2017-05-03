import pandas as pd
import numpy as np
import os
from settings import DATA_DIR, DATA_DIR_NON_GIT


NORM_METHODS = {
    None,
    'raw',
    'bmiq',
    'swan',
    'pbc',
    'funnorm'
}


def load_illumina_methylationepic_annotation():
    fn = os.path.join(DATA_DIR, 'annotation', 'methylation', 'infinium-methylationepic-v1-0-b3-manifest-file-csv.zip')
    usecols = [
        'Name', 'CHR', 'MAPINFO', 'Strand', 'UCSC_RefGene_Name',
        'UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island'
    ]
    dtype = dict(
        Name=str,
        CHR=str,
        MAPINFO=str,  # should be int but there are some NA entries
        Strand=str,
        UCSC_RefGene_Name=str,
        UCSC_RefGene_Group=str,
        Relation_to_UCSC_CpG_Island=str
    )
    dat = pd.read_csv(
        fn, skiprows=7, usecols=usecols, dtype=dtype, header=0, index_col=0
    )
    # remove calibration probes
    dat = dat.loc[~dat.loc[:, 'MAPINFO'].isnull()]
    dat.loc[:, 'MAPINFO'] = dat.loc[:, 'MAPINFO'].astype(int)

    return dat


def gbm_nsc_methylationepic(norm_method='swan'):
    """
    Load beta values. These are all generated in one pass in R.
    :param norm_method: One of NORM_METHOD
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-12-19_ucl_genomics', 'beta')
    if norm_method is None:
        fn = os.path.join(indir, 'beta_raw.csv')
    elif norm_method.lower() in NORM_METHODS:
        fn = os.path.join(indir, 'beta_%s.csv' % norm_method.lower())
    else:
        raise AttributeError("Unrecognised norm_method %s. Options are (%s)." % (
            norm_method,
            ', '.join(str(t) for t in NORM_METHODS)
        ))

    b = pd.read_csv(fn, header=0, index_col=0)
    return b


