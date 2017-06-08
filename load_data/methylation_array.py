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


def load_beta_values(indir, metafile=None, norm_method='swan', samples=None):
    """
    Load beta values. These are all generated in one pass in R using `methylation/process.R`
    :param indir: Directory containing the data files.
    :param metafile: Optionally supply the path to the metafile
    :par0am norm_method: One of NORM_METHOD
    :param samples: If supplied, this is a list (or other iterable) of sample names to keep. These must match the
    column names in the input files
    :return: Single Pandas DataFrame (if metafile is None) or (data, meta) (if metafile is supplied).
    """
    if norm_method is None:
        fn = os.path.join(indir, 'beta_raw.csv.gz')
    elif norm_method.lower() in NORM_METHODS:
        fn = os.path.join(indir, 'beta_%s.csv.gz' % norm_method.lower())
    else:
        raise AttributeError("Unrecognised norm_method %s. Options are (%s)." % (
            norm_method,
            ', '.join(str(t) for t in NORM_METHODS)
        ))

    b = pd.read_csv(fn, header=0, index_col=0)
    meta = None
    if metafile is not None:
        meta = pd.read_csv(metafile, header=0, index_col=0)

    if samples is not None:
        if len(b.columns.intersection(samples)) != len(samples):
            missing = set(samples).difference(b.columns)
            raise KeyError("Some samples were not found: %s" % ', '.join(missing))
        b = b.loc[:, samples]
        if meta is not None:
            meta = meta.loc[samples, :]

    if meta is not None:
        return (b, meta)
    else:
        return b


def gbm_rtk1_and_paired_nsc(norm_method='swan'):
    samples1 = (
        'GBM018',
        'GBM019',
        'GBM031',
        'DURA018N2 NSC',
        'DURA019N8C NSC',
        'DURA031N44B NSC',
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-12-19_ucl_genomics', 'beta')
    metafile = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-12-19_ucl_genomics', 'sources.csv')
    b1, m1 = load_beta_values(indir, metafile=metafile, norm_method=norm_method, samples=samples1)

    samples2 = (
        'GBM018_P12',
        'GBM030_P5',  # probably RTK I: find out
        'DURA018N4NSC_P4',
        'DURA030N16B6NSC_P1', # probably RTK I: find out
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-05-12', 'beta')
    metafile = os.path.join(indir, '..', 'sources.csv')
    b2, m2 = load_beta_values(indir, metafile=metafile, norm_method=norm_method, samples=samples2)

    # combine
    meta = pd.concat((m1, m2), axis=0)
    data = pd.concat((b1, b2), axis=1)

    return data, meta


def hgic_methylationepic(norm_method='swan'):
    samples1 = (
        'GBM018',
        'GBM019',
        'GBM024',
        'GBM026',
        'GBM031',
        'DURA018N2 NSC',
        'DURA019N8C NSC',
        'DURA024N28 NSC',
        'DURA026N31D NSC',
        'DURA031N44B NSC',
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-12-19_ucl_genomics', 'beta')
    metafile = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-12-19_ucl_genomics', 'sources.csv')
    b1, m1 = load_beta_values(indir, metafile=metafile, norm_method=norm_method, samples=samples1)

    samples2 = (
        'GBM044_P4',
        'GBM026_P3_P4',
        'GBM018_P12',
        'GBM044_P8',
        'DURA044N8NSC_P2',
        'GIBCONSC_P4',
        'DURA044N17NSC_P3',
        'DURA018N4NSC_P4',
        'GBM030_P5',
        'DURA030N16B6NSC_P1',
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-05-12', 'beta')
    metafile = os.path.join(indir, '..', 'sources.csv')
    b2, m2 = load_beta_values(indir, metafile=metafile, norm_method=norm_method, samples=samples2)

    # combine
    meta = pd.concat((m1, m2), axis=0)
    data = pd.concat((b1, b2), axis=1)

    return data, meta
