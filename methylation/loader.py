from load_data import loader
import pandas as pd
import os
import re
from settings import GIT_LFS_DATA_DIR, DATA_DIR_NON_GIT
from utils.log import get_console_logger
from utils import setops, string_manipulation


NORM_METHODS = {
    None,
    'raw',
    'bmiq',
    'swan',
    'pbc',
    'funnorm'
}

project_dirs = {
    "2016-06-10_brandner": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-06-10_brandner'),
    "2016-09-21_dutt": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-09-21_dutt'),
    "2016-12-19_ucl_genomics": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2016-12-19_ucl_genomics'),
    "2017-01-17_brandner": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-01-17_brandner'),
    "2017-02-09_brandner": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-02-09_brandner'),
    "2017-05-12": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-05-12'),
    "2017-08-23": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-08-23'),
    "2017-09-19": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2017-09-19'),
    "2018-01-12": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2018-01-12'),
    "2018-04-09": os.path.join(DATA_DIR_NON_GIT, 'methylation', '2018-04-09'),
    "gse38216": os.path.join(DATA_DIR_NON_GIT, 'methylation', 'GSE38216'),
    "gse65214": os.path.join(DATA_DIR_NON_GIT, 'methylation', 'GSE65214'),
    "gse67283": os.path.join(DATA_DIR_NON_GIT, 'methylation', 'GSE67283'),
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
        ('DURA061_NSC_N1_P3n4', "2018-04-09"),
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

    # correct gene symbols - look up what these should be - some kind of Excel fail?
    correction = {
        '1-Mar': 'MARCH1',
        '1-Sep': 'SEPT1',
        '10-Mar': 'MARCH10',
        '11-Mar': 'MARCH11',
        '11-Sep': 'SEPT11',
        '13-Sep': 'SEPT13',
        '2-Mar': 'MARCH2',
        '3-Mar': 'MARCH3',
        '4-Mar': 'MARCH4',
        '5-Sep': 'SEPT5',
        '6-Mar': 'MARCH6',
        '7-Mar': 'MARCH7',
        '8-Mar': 'MARCH8',
        '9-Sep': 'SEPT9',
    }
    regex = re.compile('|'.join(correction.keys()))
    corr_idx = dat.loc[:, 'UCSC_RefGene_Name'].dropna().str.contains(regex)
    corr_idx = corr_idx.index[corr_idx]
    dat.loc[corr_idx, 'UCSC_RefGene_Name'] = dat.loc[corr_idx, 'UCSC_RefGene_Name'].apply(lambda x: correction[x])

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


class IlluminaHumanMethylationLoader(loader.SingleFileLoader):
    def __init__(self, norm_method='swan', *args, **kwargs):
        self.norm_method = norm_method
        super(IlluminaHumanMethylationLoader, self).__init__(*args, **kwargs)

    def get_inputs(self):
        if self.norm_method is None:
            self.input_files = os.path.join(self.base_dir, 'beta_raw.csv.gz')
        elif self.norm_method.lower() in NORM_METHODS:
            self.input_files = os.path.join(self.base_dir, 'beta_%s.csv.gz' % self.norm_method.lower())
        else:
            raise AttributeError("Unrecognised norm_method %s. Options are (%s)." % (
                self.norm_method,
                ', '.join(str(t) for t in NORM_METHODS)
            ))

    def load_one_file(self, fn):
        return pd.read_csv(fn, header=0, index_col=0)


def load_by_patient(
        patient_ids,
        type='cell_culture',
        norm_method='swan',
        include_control=True
):
    """
    Load all RNA-Seq count data associated with the patient ID(s) supplied
    :param patient_ids: Iterable or single int or char
    :param source:
    :param include_control: If True (default) include Gibco reference NSC
    :return:
    """

    if type == "cell_culture":
        LOOKUP = PATIENT_LOOKUP_CELL
    elif type == "ffpe":
        LOOKUP = PATIENT_LOOKUP_FFPE
    else:
        raise NotImplementedError()

    cls = IlluminaHumanMethylationLoader

    # ensure patient IDs are in correct form
    if patient_ids == 'all':
        patient_ids = [t for t in LOOKUP.keys() if t != 'GIBCO']
    elif hasattr(patient_ids, '__iter__'):
        patient_ids = [t if isinstance(t, str) else ('%03d' % t) for t in patient_ids]
    else:
        if isinstance(patient_ids, str):
            patient_ids = [patient_ids]
        else:
            patient_ids = ['%03d' % patient_ids]

    if include_control and type == 'cell_culture':
        patient_ids += ['GIBCO']

    # precompute the loaders required to avoid reloading multiple times
    # we'll also take a note of the order for later reordering
    sample_order = []
    by_loader = {}
    for pid in patient_ids:
        d = LOOKUP[pid]
        for s, ldr in d:
            by_loader.setdefault(ldr, []).append(s)
            sample_order.append(s)

    objs = []
    for ldr, samples in by_loader.items():
        base_dir = os.path.join(project_dirs[ldr], 'beta')
        meta_fn = os.path.join(project_dirs[ldr], 'sources.csv')

        objs.append(
            cls(
                base_dir=base_dir,
                meta_fn=meta_fn,
                samples=samples,
                norm_method=norm_method,
                batch_id=ldr
            )
        )

    if len(objs) > 1:
        res = loader.MultipleBatchLoader(objs)
    else:
        res = objs[0]

    # apply original ordering
    res.meta = res.meta.loc[sample_order]
    res.data = res.data.loc[:, res.meta.index]

    return res


def load_reference(ref_names, norm_method='pbc', samples=None):
    ## TODO: support specifying samples in input?
    # ensure ref names are in a list (or iterable)
    if not hasattr(ref_names, '__iter__'):
        ref_names = [ref_names]

    objs = []
    for rid in ref_names:
        base_dir = project_dirs[rid]
        if not os.path.isdir(base_dir):
            raise ValueError("Directory %s for ref %s does not exist." % (base_dir, rid))

        beta_dir = os.path.join(base_dir, 'beta')
        meta_fn = os.path.join(base_dir, 'sources.csv')
        ldr = IlluminaHumanMethylationLoader(
            base_dir=beta_dir,
            meta_fn=meta_fn,
            batch_id=rid,
            norm_method=norm_method,
        )
        objs.append(ldr)

    if len(objs) > 1:
        res = loader.MultipleBatchLoader(objs)
    else:
        res = objs[0]

    if samples is not None:
        res = res.filter_by_sample_name(samples)

    return res


def gse38216(norm_method='bmiq', samples=None):
    base_dir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'GSE38216')
    beta_dir = os.path.join(base_dir, 'beta')
    meta_fn = os.path.join(base_dir, 'sources.csv')
    return IlluminaHumanMethylationLoader(
        base_dir=beta_dir,
        meta_fn=meta_fn,
        batch_id="GSE38216",
        norm_method=norm_method,
        samples=samples
    )


def gse67283(norm_method='bmiq', samples=None):
    base_dir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'GSE67283')
    beta_dir = os.path.join(base_dir, 'beta')
    meta_fn = os.path.join(base_dir, 'sources.csv')
    return IlluminaHumanMethylationLoader(
        base_dir=beta_dir,
        meta_fn=meta_fn,
        batch_id="GSE67283",
        norm_method=norm_method,
        samples=samples
    )


def gse65214(norm_method='bmiq', samples=None):
    base_dir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'GSE65214')
    beta_dir = os.path.join(base_dir, 'beta')
    meta_fn = os.path.join(base_dir, 'sources.csv')
    return IlluminaHumanMethylationLoader(
        base_dir=beta_dir,
        meta_fn=meta_fn,
        batch_id="GSE65214",
        norm_method=norm_method,
        samples = samples
    )


def encode_epic(norm_method='bmiq', samples=None):
    base_dir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'ENCODE_EPIC')
    beta_dir = os.path.join(base_dir, 'beta')
    meta_fn = os.path.join(base_dir, 'sources.csv')
    return IlluminaHumanMethylationLoader(
        base_dir=beta_dir,
        meta_fn=meta_fn,
        batch_id="Encode EPIC",
        norm_method=norm_method,
        samples=samples
    )


def encode_450k(norm_method='bmiq', samples=None):
    base_dir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'ENCODE_450k')
    beta_dir = os.path.join(base_dir, 'beta')
    meta_fn = os.path.join(base_dir, 'sources.csv')
    return IlluminaHumanMethylationLoader(
        base_dir=beta_dir,
        meta_fn=meta_fn,
        batch_id="Encode 450k",
        norm_method=norm_method,
        samples=samples
    )