import pandas as pd
import os
from settings import GIT_LFS_DATA_DIR, DATA_DIR
from utils.log import get_console_logger
from utils import setops
logger = get_console_logger(__name__)


NORM_METHODS = {
    None,
    'raw',
    'bmiq',
    'swan',
    'pbc',
    'funnorm'
}

project_dirs = {
    "2016-06-10_brandner": os.path.join(DATA_DIR, 'methylation', '2016-06-10_brandner'),
    "2016-09-21_dutt": os.path.join(DATA_DIR, 'methylation', '2016-09-21_dutt'),
    "2016-12-19_ucl_genomics": os.path.join(DATA_DIR, 'methylation', '2016-12-19_ucl_genomics'),
    "2017-01-17_brandner": os.path.join(DATA_DIR, 'methylation', '2017-01-17_brandner'),
    "2017-02-09_brandner": os.path.join(DATA_DIR, 'methylation', '2017-02-09_brandner'),
    "2017-05-12": os.path.join(DATA_DIR, 'methylation', '2017-05-12'),
    "2017-08-23": os.path.join(DATA_DIR, 'methylation', '2017-08-23'),
    "2017-09-19": os.path.join(DATA_DIR, 'methylation', '2017-09-19'),
    "2018-01-12": os.path.join(DATA_DIR, 'methylation', '2018-01-12'),
}

PATIENT_LOOKUP_FFPE = {}  # TODO?

PATIENT_LOOKUP_CELL = {
    '017': [
        ('GBM017_P3', "2017-09-19"),
        ('GBM017_P4', "2017-09-19"),
        ('DURA017_NSC_N3C5_P4', "2017-09-19"),
        ('DURA017_FB_P7', '2018-01-12'),
    ],
    '018': [
        ('GBM018_P12', '2017-05-12'),
        ('GBM018_P10', '2016-12-19_ucl_genomics'),
        ('DURA018_NSC_N4_P4', '2017-05-12'),
        ('DURA018_NSC_N2_P6', '2016-12-19_ucl_genomics'),
    ],
    '019': [
        ('GBM019_P4', '2016-12-19_ucl_genomics'),
        ('GBM019_P3n6', "2017-09-19"),
        ('DURA019_NSC_N8C_P2', '2016-12-19_ucl_genomics'),
        ('DURA019_NSC_N5C1_P2', '2018-01-12'),
        ('DURA019_FB_P7', '2018-01-12')
    ],
    '026': [
        ('GBM026_P8', '2016-12-19_ucl_genomics'),
        ('GBM026_P3n4', '2017-05-12'),
        ('DURA026_NSC_N31D_P5', '2016-12-19_ucl_genomics'),
    ],
    '030': [
        ('GBM030_P9', "2017-09-19"),
        ('GBM030_P5', '2017-05-12'),
        ('DURA030_NSC_N16B6_P1', '2017-05-12'),
        ('DURA030_NSC_N9_P2', '2018-01-12'),
        ('DURA030_FB_P8', '2018-01-12'),
    ],
    '031': [
        ('GBM031_P7', "2017-09-19"),
        ('GBM031_P4', '2016-12-19_ucl_genomics'),
        ('DURA031_NSC_N44B_P2', '2016-12-19_ucl_genomics'),
        ('DURA031_NSC_N44F_P3', '2018-01-12'),
        ('DURA031_FB_P7', '2018-01-12'),
    ],
    '044': [
        ('GBM044_P4', '2017-05-12'),
        ('GBM044_P8', '2017-05-12'),
        ('DURA044_NSC_N17_P3', '2017-05-12'),
        ('DURA044_NSC_N8_P2', '2017-05-12'),
    ],
    '049': [
        ('GBM049_P4', "2017-08-23"),
        ('GBM049_P6', "2017-08-23"),
        ('DURA049_NSC_N19_P4', "2017-08-23"),
        ('DURA049_NSC_N5_P2', "2017-08-23"),
        ('DURA049_IPSC_ N5_P10', '2018-01-12'),
    ],
    '050': [
        ('GBM050_P7n8', "2017-08-23"),
        ('GBM050_P9', "2017-08-23"),
        ('DURA050_NSC_N12_P3', "2017-08-23"),
        ('DURA050_NSC_N16_P4', "2017-08-23"),
        ('DURA050_IPSC_N12_P5', "2018-01-12"),
        ('DURA050_FB_P7', "2018-01-12"),
    ],
    '052': [
        ('GBM052_P6n7', "2017-09-19"),
        ('GBM052_P4n5', "2017-09-19"),
        ('DURA052_NSC_N4_P3', "2017-09-19"),
        ('DURA052_NSC_N5_P2', "2017-09-19"),
    ],
    '054': [
        ('GBM054_P4', "2017-08-23"),
        ('GBM054_P6', "2017-08-23"),
        ('DURA054_NSC_N3C_P2', "2017-08-23"),
        ('DURA054_NSC_N2E_P1', "2017-08-23"),
        ('DURA054_IPSC_N3C_P11', '2018-01-12'),
        ('DURA054_FB_P5', '2018-01-12'),
    ],
    '061': [
        ('GBM061_P3', "2017-08-23"),
        ('GBM061_P5', "2017-08-23"),
        ('DURA061_NSC_N4_P2', "2017-08-23"),
        ('DURA061_NSC_N6_P4', "2017-08-23"),
    ],
    'GIBCO': [
        ('GIBCONSC_P4', '2017-05-12'),
    ]
}


def load_illumina_methylationepic_annotation(split_genes=True):
    """

    :param split_genes: If True (default), the RefGene name column will be split into a set - useful for
    many downstream applications
    :return:
    """
    fn = os.path.join(GIT_LFS_DATA_DIR, 'annotation', 'methylation', 'infinium-methylationepic-v1-0-b3-manifest-file-csv.zip')
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

    if split_genes:
        dat.loc[:, 'UCSC_RefGene_Name'] = \
            dat.UCSC_RefGene_Name.str.split(';').apply(lambda x: set(x) if isinstance(x, list) else None)

    return dat


def load_illumina_methylation450_annotation():
    fn = os.path.join(GIT_LFS_DATA_DIR, 'annotation', 'methylation', 'GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz')
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


## TODO: unify the RNASeq batch combiner with this one?
class IlluminaHumanMethylationLoader(object):
    # for compatibility with RNASeq loader
    annotate_by = None
    def __init__(self, indir, meta_fn, norm_method='swan', samples=None):
        self.indir = indir
        self.meta_fn = meta_fn
        self.norm_method = norm_method
        self.samples = samples
        self.data = None
        self.meta = None
        self.load_data()

    def load_data(self):
        if self.norm_method is None:
            fn = os.path.join(self.indir, 'beta_raw.csv.gz')
        elif self.norm_method.lower() in NORM_METHODS:
            fn = os.path.join(self.indir, 'beta_%s.csv.gz' % self.norm_method.lower())
        else:
            raise AttributeError("Unrecognised norm_method %s. Options are (%s)." % (
                self.norm_method,
                ', '.join(str(t) for t in NORM_METHODS)
            ))

        self.data = pd.read_csv(fn, header=0, index_col=0)
        if self.meta_fn is not None:
            self.meta = pd.read_csv(self.meta_fn, header=0, index_col=0)

        if self.samples is not None:
            if len(self.data.columns.intersection(self.samples)) != len(self.samples):
                missing = set(self.samples).difference(self.data.columns)
                raise KeyError("Some samples were not found: %s" % ', '.join(missing))
            self.data = self.data.loc[:, self.samples]
            if self.meta is not None:
                self.meta = self.meta.loc[self.samples, :]
                # for compatibility with RNASeq batch loader need to duplicate sample names in a column called 'sample'
                if 'sample' not in self.meta.columns:
                    self.meta.insert(0, 'sample', self.meta.index)


def load_by_patient(patient_ids, norm_method='swan', include_control=True):
    """
    Load data for one or more patients, specified by ID.
    :param patient_ids: Iterable or single int or char
    :param include_control: If True (default), alos load Gibco NSC control data.
    """
    # ensure patient IDs are in correct form
    if patient_ids == 'all':
        patient_ids = [t for t in PATIENT_LOOKUP_CELL.keys() if t != 'GIBCO']
    elif hasattr(patient_ids, '__iter__'):
        patient_ids = [t if isinstance(t, str) else ('%03d' % t) for t in patient_ids]
    else:
        if isinstance(patient_ids, str):
            patient_ids = [patient_ids]
        else:
            patient_ids = ['%03d' % patient_ids]

    if include_control:
        patient_ids += ['GIBCO']

    # precompute the loaders required to avoid reloading multiple times
    # we'll also take a note of the order for later reordering
    sample_order = []
    by_loader = {}
    for pid in patient_ids:
        d = PATIENT_LOOKUP_CELL[pid]
        for s, ldr in d:
            by_loader.setdefault(ldr, []).append(s)
            sample_order.append(s)

    # loaders = []

    data = []
    metas = []

    for ldr, samples in by_loader.items():
        indir = project_dirs[ldr]
        beta_dir = os.path.join(indir, 'beta')
        meta_fn = os.path.join(indir, 'sources.csv')
        # loaders.append(IlluminaHumanMethylationLoader(beta_dir, meta_fn, norm_method=norm_method, samples=samples))
        b, m = load_beta_values(beta_dir, metafile=meta_fn, norm_method=norm_method, samples=samples)
        data.append(b)
        metas.append(m)

    # combine them, keeping only matching probes
    probes = setops.reduce_intersection(*[t.index for t in data])

    # report dropped probes
    dropped_probes = [t.index.difference(probes).size for t in data]
    if any([t > 0 for t in dropped_probes]):
        logger.warning(
            "Dropped some probes when combining the data from different batches: %s",
            ', '.join([str(t) for t in dropped_probes])
        )
    # TODO: handle sample name clash
    beta = pd.concat((t.loc[probes] for t in data), axis=1)
    meta_cols = reduce(lambda x, y: x.union(y), [t.columns for t in metas])
    meta = pd.DataFrame(columns=meta_cols, index=sample_order)
    for t in metas:
        meta.loc[t.index, t.columns] = t

    return beta, meta


def gse36278(dropna=True):
    """
    GBM dataset published by Sturm et al., supporting their study on classifying GBM subgroups.
    The data were provided as average beta values per probe. Different normalisation methods are therefore not
    available.
    :return:
    """
    indir = os.path.join(DATA_DIR, 'methylation', 'GSE36278')
    metafile = os.path.join(indir, 'sources.csv')
    data, meta = load_beta_values(indir, metafile=metafile, norm_method='raw')
    # data.columns = meta.loc[data.columns, 'title']
    if dropna:
        data = data.dropna()
    return data, meta


def gbm_rtk1_and_paired_nsc(norm_method='swan', ref=None):
    """

    :param norm_method:
    :param ref: Default None means no reference is loaded. If provided, this argument is a string containing a supported
     reference. Currently only 'gibco' is supported.
    :return:
    """
    samples1 = (
        'GBM018_P10',
        'GBM019_P4',
        'GBM031_P4',
        'DURA018_NSC_N2_P6',
        'DURA019_NSC_N8C_P2',
        'DURA031_NSC_N44B_P2',
    )
    indir = os.path.join(DATA_DIR, 'methylation', '2016-12-19_ucl_genomics', 'beta')
    metafile = os.path.join(DATA_DIR, 'methylation', '2016-12-19_ucl_genomics', 'sources.csv')
    b1, m1 = load_beta_values(indir, metafile=metafile, norm_method=norm_method, samples=samples1)

    samples2 = (
        'GBM018_P12',
        'GBM030_P5',  # probably RTK I: find out
        'DURA018_NSC_N4_P4',
        'DURA030_NSC_N16B6_P1', # probably RTK I: find out
    )
    if ref is not None:
        if ref == 'gibco':
            samples2 += ('GIBCONSC_P4', )
        else:
            raise NotImplementedError("Unrecognised reference %s" % ref)


    indir = os.path.join(DATA_DIR, 'methylation', '2017-05-12', 'beta')
    metafile = os.path.join(indir, '..', 'sources.csv')
    b2, m2 = load_beta_values(indir, metafile=metafile, norm_method=norm_method, samples=samples2)

    # combine
    meta = pd.concat((m1, m2), axis=0)
    data = pd.concat((b1, b2), axis=1).dropna()

    return data, meta


def hgic_methylationepic(norm_method='swan'):
    samples1 = (
        'GBM018_P10',
        'GBM019_P4',
        'GBM024_P9',
        'GBM026_P8',
        'GBM031_P4',
        'DURA018_NSC_N2_P6',
        'DURA019_NSC_N8C_P2',
        'DURA024_NSC_N28_P6',
        'DURA026_NSC_N31D_P5',
        'DURA031_NSC_N44B_P2',
    )
    indir = os.path.join(DATA_DIR, 'methylation', '2016-12-19_ucl_genomics', 'beta')
    metafile = os.path.join(DATA_DIR, 'methylation', '2016-12-19_ucl_genomics', 'sources.csv')
    b1, m1 = load_beta_values(indir, metafile=metafile, norm_method=norm_method, samples=samples1)

    samples2 = (
        'GBM044_P4',
        'GBM026_P3n4',
        'GBM018_P12',
        'GBM044_P8',
        'DURA044_NSC_N8_P2',
        'GIBCONSC_P4',
        'DURA044_NSC_N17_P3',
        'DURA018_NSC_N4_P4',
        'GBM030_P5',
        'DURA030_NSC_N16B6_P1',
    )
    indir = os.path.join(DATA_DIR, 'methylation', '2017-05-12', 'beta')
    metafile = os.path.join(indir, '..', 'sources.csv')
    b2, m2 = load_beta_values(indir, metafile=metafile, norm_method=norm_method, samples=samples2)

    # combine
    meta = pd.concat((m1, m2), axis=0)
    data = pd.concat((b1, b2), axis=1).dropna()

    return data, meta
