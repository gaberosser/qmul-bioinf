import pandas as pd
import re
import os
import references
from load_data import loader
from load_data.loader import MultipleBatchLoader
from rnaseq.general import ensembl_transcript_quant_to_gene
from utils.log import get_console_logger
from utils import setops
from settings import RNASEQ_DIR
logger = get_console_logger(__name__)


INDEX_FIELDS = (
    'Approved Symbol',
    'Entrez Gene ID',
    'RefSeq IDs',
    'Ensembl Gene ID'
)

class RnaSeqFileLocations(object):
    def __init__(self, root_dir, alignment_subdir=None, batch_id = None, strandedness='r', tax_id=9606):
        self.root_dir = root_dir
        self.strandedness = strandedness
        self.alignment_subdir = alignment_subdir
        if batch_id is None:
            self.batch_id = os.path.split(self.root_dir)[-1]
        else:
            self.batch_id = batch_id
        self.tax_id = tax_id

        self.base_dir = self.root_dir if alignment_subdir is None else os.path.join(self.root_dir, self.alignment_subdir)
        self.meta_file = os.path.join(self.root_dir, 'sources.csv')

        common_kwds = {
            'strandedness': 'r',
            'tax_id': self.tax_id,
            'batch_id': self.batch_id,
            'meta_fn': self.meta_file
        }

        self.params = {}

        self.params['salmon'] = {
            'base_dir': os.path.join(self.base_dir, 'salmon')
        }

        self.params['star'] = {
            'base_dir': os.path.join(self.base_dir, 'star_alignment')
        }

        self.params['star_cufflinks'] = {
            'base_dir': os.path.join(self.base_dir, 'star_alignment', 'cufflinks')
        }

        for k, v in self.params.items():
            v.update(common_kwds)

    def loader_kwargs(self, typ):
        if typ == 'star':
            return self.params['star']
        elif typ == 'salmon':
            return self.params['salmon']
        elif typ == 'star/cufflinks':
            return self.params['star_cufflinks']
        else:
            raise NotImplementedError("Unrecognised type: %s" % typ)


wtchg_p160704 = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p160704'),
    alignment_subdir='human',
)

wtchg_p170218 = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p170218'),
    alignment_subdir='human',
)

wtchg_p170390 = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p170390'),
    alignment_subdir='human',
)

wtchg_p170503 = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p170503'),
    alignment_subdir='human',
)

wtchg_p160704_ribozero = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p160704_ribozero'),
    alignment_subdir='human',
)

wtchg_p160704_ribozero2 = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p160704_ribozero_rerun'),
    alignment_subdir='human',
)

wtchg_p170446 = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p170446'),
    alignment_subdir='human',
)

wtchg_p170446b = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p170446b'),
    alignment_subdir='human',
)

wtchg_p170582 = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p170582'),
    alignment_subdir='human',
)

wtchg_p170644 = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p170644'),
    alignment_subdir='human',
)

wtchg_p180199 = RnaSeqFileLocations(
    root_dir=os.path.join(RNASEQ_DIR, 'wtchg_p180199'),
    alignment_subdir='human',
)

PATIENT_LOOKUP_CC = {
    '017': [
        ('GBM017_P3', wtchg_p170390),
        ('GBM017_P4', wtchg_p170390),
        ('DURA017_NSC_N3C5_P4', wtchg_p170390),
        ('DURA017_FB_P7', wtchg_p180199),
    ],
    '018': [
        ('GBM018_P12', wtchg_p170218),
        ('GBM018_P10', wtchg_p160704),
        ('DURA018_NSC_N4_P4', wtchg_p170218),
        ('DURA018_NSC_N2_P6', wtchg_p160704),
        ('DURA018_FB_P6', wtchg_p180199),
    ],
    '019': [
        ('GBM019_P4', wtchg_p160704),
        ('GBM019_P3n6', wtchg_p170390),
        ('DURA019_NSC_N8C_P2', wtchg_p160704),
        ('DURA019_NSC_N5C1_P2', wtchg_p170582),
        ('DURA019_IPSC_N8_P13', wtchg_p180199),
    ],
    '024': [
        ('GBM024_P9', wtchg_p160704),
        ('DURA024_NSC_N28_P6', wtchg_p160704),
    ],
    '026': [
        ('GBM026_P8', wtchg_p160704),
        ('GBM026_P3n4', wtchg_p170218),
        ('DURA026_NSC_N31D_P5', wtchg_p160704),
        ('DURA026_FB_P8', wtchg_p180199),
    ],
    '030': [
        ('GBM030_P9n10', wtchg_p170390),
        ('GBM030_P5', wtchg_p170218),
        ('DURA030_NSC_N16B6_P1', wtchg_p170218),
        ('DURA030_NSC_N9_P2', wtchg_p170582),
        ('DURA030_IPSC_N16B6_P13', wtchg_p180199),
    ],
    '031': [
        ('GBM031_P7', wtchg_p170390),
        ('GBM031_P4', wtchg_p160704),
        ('DURA031_NSC_N44B_P2', wtchg_p160704),
        ('DURA031_NSC_N44_P3', wtchg_p170582),
        ('DURA031_IPSC_N44B_P10', wtchg_p180199),
    ],
    '044': [
        ('GBM044_P4', wtchg_p170218),
        ('GBM044_P8', wtchg_p170218),
        ('DURA044_NSC_N17_P3', wtchg_p170218),
        ('DURA044_NSC_N8_P2', wtchg_p170218),
    ],
    '047': [
        ('GBM047_P8', wtchg_p180199),
    ],
    '049': [
        # ('GBM049_P4', wtchg_p170503),  # IDH1 mutant switchover
        # ('GBM049_P6', wtchg_p170503),  # IDH1 mutant switchover
        ('GBM049_P7', wtchg_p180199),
        ('GBM049_P9', wtchg_p180199),
        ('DURA049_NSC_N19_P4', wtchg_p170503),
        ('DURA049_NSC_N5_P2', wtchg_p170503),
        ('DURA049_IPSC_N5_P10', wtchg_p170582)
    ],
    '050': [
        ('GBM050_P7n8', wtchg_p170503),
        ('GBM050_P9', wtchg_p170503),
        ('DURA050_NSC_N12_P3', wtchg_p170503),
        ('DURA050_NSC_N16_P4', wtchg_p170503),
        ('DURA050_IPSC_N12_P5', wtchg_p170582),
        ('DURA050_FB_P7', wtchg_p180199),
    ],
    '052': [
        ('GBM052_P6n7', wtchg_p170503),
        ('GBM052_P4n5', wtchg_p170503),
        ('DURA052_NSC_N4_P3', wtchg_p170503),
        ('DURA052_NSC_N5_P2', wtchg_p170503),
        ('DURA052_FB_P6', wtchg_p180199)
    ],
    '054': [
        ('GBM054_P4', wtchg_p170503),
        ('GBM054_P6', wtchg_p170503),
        ('DURA054_NSC_N3C_P2', wtchg_p170503),
        ('DURA054_NSC_N2E_P1', wtchg_p170503),
        ('DURA054_IPSC_N3C_P11', wtchg_p170582),
        ('DURA054_FB_P5', wtchg_p180199),
    ],
    '061': [
        ('GBM061_P3', wtchg_p170503),
        ('GBM061_P5', wtchg_p170503),
        ('DURA061_NSC_N4_P2', wtchg_p170503),
        ('DURA061_NSC_N6_P4', wtchg_p170503),
        ('DURA061_NSC_N1_P3', wtchg_p170644),
        ('DURA061_NSC_N1_P5', wtchg_p170644),
        ('DURA061_IPSC_N4_P5', wtchg_p170582),
        ('DURA061_FB_P7', wtchg_p180199),
    ],
    'GIBCO': [
        ('GIBCO_NSC_P4', wtchg_p170218),
    ],
    '3021': [
        ('3021_1Scr', wtchg_p170644),
        ('3021_1shB', wtchg_p170644),
        ('3021_1shC', wtchg_p170644),
        ('3021_1shB+C', wtchg_p170644),
        ('3021_2Scr', wtchg_p170644),
        ('3021_2shB', wtchg_p170644),
        ('3021_2shC', wtchg_p170644),
        ('3021_2shB+C', wtchg_p170644),
    ],
    'ICb1299': [
        ('shScr', wtchg_p170218),
        ('shBmi1', wtchg_p170218),
        ('shChd7', wtchg_p170218),
        ('shB+C', wtchg_p170218),
        ('1299_S', wtchg_p170582),
        ('1299_B', wtchg_p170582),
        ('1299_C', wtchg_p170582),
        ('1299_B+C', wtchg_p170582),
        ('ICb1299_3Scr', wtchg_p180199),
        ('ICb1299_3shBmi1', wtchg_p180199),
        ('ICb1299_3shChd7', wtchg_p180199),
        ('ICb1299_3shB+C', wtchg_p180199),
    ]
}

PATIENT_LOOKUP_FFPE = {
    '017': [
        ('NH15_1661DEF2C', wtchg_p170446)
    ],
    '018': [
        ('NH15_1877_SP1C', wtchg_p160704_ribozero)
    ],
    '019': [
        ('NH15_2101_DEF1A', wtchg_p160704_ribozero)
    ],
    '026': [
        ('NH16_270_DEF1A', wtchg_p160704_ribozero2),
        ('NH16_270_DEF1Ereplacement', wtchg_p170446b)
    ],
    '030': [
        ('NH16_616DEF1B', wtchg_p170446)
    ],
    '031': [
        ('NH16_677_SP1A', wtchg_p160704_ribozero)
    ],
    '044': [
        ('NH16_1574DEF1A', wtchg_p170446)
    ],
    '049': [
        ('NH16_1976DEF2A', wtchg_p170446),
        ('NH16_1976_DEF1Areplacement', wtchg_p170446b),
    ],
    '050': [
        ('NH16_2063DEF1B1', wtchg_p170446),
        ('NH16_2063_DEF1Areplacement', wtchg_p170446b),
    ],
    '052': [
        ('NH16_2214DEF1A', wtchg_p170446),
        ('NH16_2214_DEF1A', wtchg_p170446b),
    ],
    '054': [
        ('NH16_2255DEF1B2', wtchg_p170446)
    ],
    '061': [
        ('NH16_2806DEF3A1', wtchg_p170446)
    ],
}


class StarCountLoader(loader.MultipleFileLoader):
    file_pattern = '*ReadsPerGene.out.tab'
    extra_df_attributes = ('data_unassigned',)

    def __init__(self, strandedness='u', *args, **kwargs):
        if strandedness not in ('u', 'f', 'r'):
            raise ValueError("Variable 'strandedness' must be one of ('u', 'f', 'r')")
        self.strandedness = strandedness
        self.data_unassigned = None
        super(StarCountLoader, self).__init__(*args, **kwargs)

    def generate_input_path(self, fname):
        """
        Given the filename from the meta file, generate the path to the actual data (e.g. constructing subdir structure)
        """
        return os.path.join(self.base_dir, "%sReadsPerGene.out.tab" % fname)

    def generate_sample_name(self, file_path):
        ff = os.path.split(file_path)[-1]
        return ff.replace('ReadsPerGene.out.tab', '')

    def load_one_file(self, fn):
        dat = pd.read_csv(fn, sep='\t', index_col=0, header=None)
        # get correct strand data
        if self.strandedness == 'u':
            return dat.iloc[:, 0]
        elif self.strandedness == 'f':
            return dat.iloc[:, 1]
        else:
            return dat.iloc[:, 2]

    def post_process(self):
        super(StarCountLoader, self).post_process()
        unal = ['N_unmapped', 'N_multimapping', 'N_noFeature', 'N_ambiguous']
        self.data_unassigned = self.data.loc[self.data.index.isin(unal)]
        self.data = self.data.loc[~self.data.index.isin(unal)]

    def filter_samples(self, keep_idx):
        """
        Drop samples according to the supplied index. All data structures are updated in-place.
        :param drop_idx: Either a boolean index or an iterable that can be used as an indexer. The samples indicated
        will be retained.
        :return:
        """
        super(StarCountLoader, self).filter_samples(keep_idx)
        self.data_unassigned = self.data_unassigned.loc[:, self.meta.index]


class SalmonQuantLoader(loader.MultipleFileLoader):
    file_pattern = 'quant.sf'

    def __init__(self, units='tpm', aggregate_to_gene_level=True, *args, **kwargs):
        if units not in ('tpm', 'estimated_counts'):
            raise ValueError("Unsupported units %s" % units)
        self.units = units
        self.aggregate_to_gene_level = aggregate_to_gene_level
        self.transcript_lengths = None
        self.effective_transcript_lengths = None
        super(SalmonQuantLoader, self).__init__(*args, **kwargs)

    def generate_input_path(self, fname):
        """
        Given the filename from the meta file, generate the path to the actual data (e.g. constructing subdir structure)
        """
        return os.path.join(self.base_dir, fname, "quant.sf")

    def generate_sample_name(self, file_path):
        # first split off the quant.sf filename
        ff, _ = os.path.split(file_path)
        # now look at the subdirectory name
        return os.path.split(ff)[-1]

    def load_one_file(self, fn):
        dat = pd.read_csv(fn, sep='\t', index_col=0, header=0)
        if self.units == 'tpm':
            return dat.loc[:, 'TPM']
        elif self.units == 'estimated_counts':
            return dat.loc[:, 'NumReads']
        raise NotImplementedError("Unsupported units")

    def post_process(self):
        super(SalmonQuantLoader, self).post_process()
        # aggregate if required
        if self.aggregate_to_gene_level:
            self.data = ensembl_transcript_quant_to_gene(self.data, tax_id=self.tax_id)
        # reload one file to get additional data (transcript lengths)
        # only do this if we haven't aggregated
        if not self.aggregate_to_gene_level:
            if self.verbose:
                self.logger.info("Reloading input file %s to get transcript lengths", self.input_files.iloc[0])
            dat = pd.read_csv(self.input_files.iloc[0], sep='\t', index_col=0, header=0)
            self.transcript_lengths = dat.loc[:, 'Length']
            self.effective_transcript_lengths = dat.loc[:, 'EffectiveLength']


class StarCufflinksLoader(loader.MultipleFileLoader):
    file_pattern = "genes.fpkm_tracking"

    def __init__(self, remove_non_canonical=True, *args, **kwargs):
        self.gene_lengths = None
        self.remove_non_canonical = remove_non_canonical
        super(StarCufflinksLoader, self).__init__(*args, **kwargs)

    def generate_input_path(self, fname):
        """
        Given the filename from the meta file, generate the path to the actual data (e.g. constructing subdir structure)
        """
        res = os.path.join(self.base_dir, "%sAligned.sortedByCoord.out" % fname, "genes.fpkm_tracking")
        return res

    def generate_sample_name(self, file_path):
        # first split off the genes.fpkm_tracking filename
        ff, _ = os.path.split(file_path)
        # now look at the subdirectory name
        ff = os.path.split(ff)[-1]
        res = ff.replace("Aligned.sortedByCoord.out", "")
        return res

    def load_one_file(self, fn):
        dat = pd.read_csv(fn, sep='\t', index_col=0, header=0)
        # need to split the results into duplicated and unique genes, then resolve the duplicates
        dupe_idx = dat.index[dat.index.duplicated()].unique()
        dupe_fpkm = dat.loc[dupe_idx, 'FPKM']
        dupe_fpkm = dupe_fpkm.groupby(dupe_fpkm.index).max()
        unique_fpkm = dat.loc[~dat.index.duplicated(False), 'FPKM'].copy()

        return pd.concat((unique_fpkm, dupe_fpkm))

    def post_process(self):
        super(StarCufflinksLoader, self).post_process()
        if self.remove_non_canonical:
            # only keep known genes
            self.data = self.data.loc[~self.data.index.str.contains('CUFF.')]
        if self.verbose:
            self.logger.info("Reloading input file %s to get gene lengths", self.input_files.iloc[0])
        dat = pd.read_csv(self.input_files.iloc[0], sep='\t', index_col=0, header=0)
        self.gene_lengths = dat.loc[self.data.index, 'length']


def load_by_patient(
        patient_ids,
        type='cell_culture',
        source='star',
        include_control=True,
        **kwargs
):
    """
    Load all RNA-Seq count data associated with the patient ID(s) supplied
    :param patient_ids: Iterable or single int or char
    :param source:
    :param include_control: If True (default) include Gibco reference NSC
    :param kwargs: Passed to the loader
    :return:
    """

    if type == "cell_culture":
        LOOKUP = PATIENT_LOOKUP_CC
    elif type == "ffpe":
        LOOKUP = PATIENT_LOOKUP_FFPE
    else:
        raise NotImplementedError()

    if source == 'star':
        cls = StarCountLoader
    elif source == 'salmon':
        cls = SalmonQuantLoader
    elif source == 'star/cufflinks':
        cls = StarCufflinksLoader
    else:
        raise NotImplementedError()

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
        the_kwargs = dict(ldr.loader_kwargs(source))
        the_kwargs.update(kwargs)
        objs.append(
            cls(
                samples=samples,
                **the_kwargs
            )
        )

    if len(objs) > 1:
        res = loader.MultipleBatchLoader(objs)
    else:
        res = objs[0]
        # make samples column the meta index
        # res.meta.set_index('sample', inplace=True)

    # apply original ordering
    res.meta = res.meta.loc[sample_order]
    res.data = res.data.loc[:, res.meta.index]

    return res


def load_references(
    ref_names,
    source='star',
    tax_id=9606,
    batch_names=None,
    **kwargs
):
    """
    Load one or more reference samples
    Optionally add batch names to each
    """
    if source == 'star':
        cls = StarCountLoader
    elif source == 'salmon':
        cls = SalmonQuantLoader
    elif source == 'star/cufflinks':
        cls = StarCufflinksLoader
    else:
        raise NotImplementedError()

    if tax_id == 9606:
        alignment_subdir = 'human'
    elif tax_id == 10090:
        alignment_subdir = 'mouse'
    else:
        raise ValueError("Unrecognised tax_id %d" % tax_id)

    # ensure patient IDs are in correct form
    if not hasattr(ref_names, '__iter__'):
        ref_names = [ref_names]

    if batch_names is not None:
        if len(batch_names) != len(ref_names):
            raise ValueError("If batch_names are supplied, they must be the same length as ref_names.")

    objs = []
    for i, rid in enumerate(ref_names):
        bid = rid if batch_names is None else batch_names[i]
        loc = RnaSeqFileLocations(
            root_dir=os.path.join(RNASEQ_DIR, rid),
            alignment_subdir=alignment_subdir,
            batch_id=bid
        )
        if not os.path.isdir(loc.root_dir):
            raise ValueError("Directory %s for ref %s does not exist." % (loc.root_dir, rid))
        this_kwargs = dict(loc.loader_kwargs(source))
        this_kwargs.update(kwargs)
        this_kwargs['tax_id'] = tax_id
        this_obj = cls(**this_kwargs)
        if batch_names is not None:
            new_idx = ["%s (%s)" % (t, batch_names[i]) for t in this_obj.meta.index]
            this_obj.meta.index = new_idx
            this_obj.data.columns = new_idx
        objs.append(this_obj)

    if len(objs) > 1:
        res = loader.MultipleBatchLoader(objs)
    else:
        res = objs[0]

    return res
