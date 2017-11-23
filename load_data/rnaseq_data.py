import pandas as pd
import re
import os
import glob
import references
from rnaseq import normalisation, tcga
from utils.log import get_console_logger
from settings import GIT_LFS_DATA_DIR, DATA_DIR_NON_GIT, DATA_DIR_NON_GIT2
logger = get_console_logger(__name__)

INDEX_FIELDS = (
    'Approved Symbol',
    'Entrez Gene ID',
    'RefSeq IDs',
    'Ensembl Gene ID'
)

class RnaSeqStarFileLocations(object):
    def __init__(self, root_dir, lanes, alignment_subdir=None, strandedness='r'):
        self.root_dir = root_dir
        self.strandedness = strandedness
        self.alignment_subdir = alignment_subdir
        self.lane_dirs = [os.path.join(root_dir, l) for l in lanes]
        self.meta_files = [os.path.join(d, 'sources.csv') for d in self.lane_dirs]
        if alignment_subdir is None:
            self.count_dirs = [os.path.join(d, 'star_alignment') for d in self.lane_dirs]
        else:
            self.count_dirs = [os.path.join(d, alignment_subdir, 'star_alignment') for d in self.lane_dirs]


    @property
    def loader_kwargs(self):
        return {
            'count_dirs': self.count_dirs,
            'meta_fns': self.meta_files,
            'strandedness': self.strandedness
        }


wtchg_p160704 = RnaSeqStarFileLocations(
    root_dir=os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704'),
    lanes=['161222_K00198_0152_AHGYG3BBXX', '161219_K00198_0151_BHGYHTBBXX'],
)

wtchg_p170218 = RnaSeqStarFileLocations(
    root_dir=os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170218'),
    lanes=['170509_K00150_0192_BHJKCLBBXX', '170515_K00150_0196_BHJKC5BBXX_lane_2', '170515_K00150_0196_BHJKC5BBXX_lane_3'],
    alignment_subdir='human'
)

wtchg_p170390 = RnaSeqStarFileLocations(
    root_dir=os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170390'),
    lanes=[
        '170727_K00198_0222_AHKWW5BBXX',
        '170731_K00150_0226_AHL2CJBBXX_1',
        '170731_K00150_0226_AHL2CJBBXX_2'
    ],
    alignment_subdir='human'
)

wtchg_p170503 = RnaSeqStarFileLocations(
    root_dir=os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170503'),
    lanes=[
        '170929_K00150_0250_BHLGNHBBXX',
        '171003_K00198_0242_AHLGYVBBXX_1',
        '171003_K00198_0242_AHLGYVBBXX_2'
    ],
)

wtchg_p160704_ribozero = RnaSeqStarFileLocations(
    root_dir=os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704_ribozero'),
    lanes=[
        '170328_K00150_0177_BHJ2C2BBXX',
    ],
)

wtchg_p160704_ribozero2 = RnaSeqStarFileLocations(
    root_dir=os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704_ribozero_rerun'),
    lanes=[
        '170905_K00150_0238_BHLFMVBBXX',
        '170907_K00150_0239_AHLFL2BBXX',
    ],
)

wtchg_p170446 = RnaSeqStarFileLocations(
    root_dir=os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170446'),
    lanes=[
        '170905_K00150_0238_BHLFMVBBXX',
        '170907_K00150_0239_AHLFL2BBXX',
    ],
)

wtchg_p170582 = RnaSeqStarFileLocations(
    root_dir=os.path.join(DATA_DIR_NON_GIT2, 'rnaseq', 'wtchg_p170582'),
    lanes=[
        '171101_K00198_0249_BHLGFJBBXX',
        '171114_K00150_0261_AHM5WJBBXX',
    ],
    alignment_subdir='human'
)

PATIENT_LOOKUP_CC_STAR = {
    '017': [
        ('GBM017_P3', wtchg_p170390),
        ('GBM017_P4', wtchg_p170390),
        ('DURA017_NSC_N3C5_P4', wtchg_p170390),
    ],
    '018': [
        ('GBM018_P12', wtchg_p170218),
        ('GBM018_P10', wtchg_p160704),
        ('DURA018_NSC_N4_P4', wtchg_p170218),
        ('DURA018_NSC_N2_P6', wtchg_p160704),
    ],
    '019': [
        ('GBM019_P4', wtchg_p160704),
        ('GBM019_P3n6', wtchg_p170390),
        ('DURA019_NSC_N8C_P2', wtchg_p160704),
        ('DURA019_NSC_N5C1_P2', wtchg_p170582),
    ],
    '024': [
        ('GBM024_P9', wtchg_p160704),
        ('DURA024_NSC_N28_P6', wtchg_p160704),
    ],
    '026': [
        ('GBM026_P8', wtchg_p160704),
        ('GBM026_P3n4', wtchg_p170218),
        ('DURA026_NSC_N31D_P5', wtchg_p160704),
    ],
    '030': [
        ('GBM030_P9n10', wtchg_p170390),
        ('GBM030_P5', wtchg_p170218),
        ('DURA030_NSC_N16B6_P1', wtchg_p170218),
        ('DURA030_NSC_N9_P2', wtchg_p170582),
    ],
    '031': [
        ('GBM031_P7', wtchg_p170390),
        ('GBM031_P4', wtchg_p160704),
        ('DURA031_NSC_N44B_P2', wtchg_p160704),
        ('DURA031_NSC_N44_P3', wtchg_p170582),
    ],
    '044': [
        ('GBM044_P4', wtchg_p170218),
        ('GBM044_P8', wtchg_p170218),
        ('DURA044_NSC_N17_P3', wtchg_p170218),
        ('DURA044_NSC_N8_P2', wtchg_p170218),
    ],
    '049': [
        ('GBM049_P4', wtchg_p170503),
        ('GBM049_P6', wtchg_p170503),
        ('DURA049_NSC_N19_P4', wtchg_p170503),
        ('DURA049_NSC_N5_P2', wtchg_p170503),
    ],
    '050': [
        ('GBM050_P7n8', wtchg_p170503),
        ('GBM050_P9', wtchg_p170503),
        ('DURA050_NSC_N12_P3', wtchg_p170503),
        ('DURA050_NSC_N16_P4', wtchg_p170503),
    ],
    '052': [
        ('GBM052_P6n7', wtchg_p170503),
        ('GBM052_P4n5', wtchg_p170503),
        ('DURA052_NSC_N4_P3', wtchg_p170503),
        ('DURA052_NSC_N5_P2', wtchg_p170503),
    ],
    '054': [
        ('GBM054_P4', wtchg_p170503),
        ('GBM054_P6', wtchg_p170503),
        ('DURA054_NSC_N3C_P2', wtchg_p170503),
        ('DURA054_NSC_N2E_P1', wtchg_p170503),
    ],
    '061': [
        ('GBM061_P3', wtchg_p170503),
        ('GBM061_P5', wtchg_p170503),
        ('DURA061_NSC_N4_P2', wtchg_p170503),
        ('DURA061_NSC_N6_P4', wtchg_p170503),
    ],
    'GIBCO': [
        ('GIBCO_NSC_P4', wtchg_p170218),
    ]
}

PATIENT_LOOKUP_FFPE_STAR = {
    '017': [
        ('NH15_1661DEF2C', wtchg_p170446)
    ],
    '018': [
        ('NH15_1877_SP1C', wtchg_p160704_ribozero)
        # ('GBM018 FFPE', wtchg_p160704_ribozero)
    ],
    '019': [
        ('NH15_2101_DEF1A', wtchg_p160704_ribozero)
        # ('GBM019 FFPE', wtchg_p160704_ribozero)
    ],
    '026': [
        ('NH16_270_DEF1A', wtchg_p160704_ribozero2)
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
        ('NH16_1976DEF2A', wtchg_p170446)
    ],
    '050': [
        ('NH16_2063DEF1B1', wtchg_p170446)
    ],
    '052': [
        ('NH16_2214DEF1A', wtchg_p170446)
    ],
    '054': [
        ('NH16_2255DEF1B2', wtchg_p170446)
    ],
    '061': [
        ('NH16_2806DEF3A1', wtchg_p170446)
    ],
}


def strip_extension(s, file_ext):
    regex = file_ext.replace('.', '\.')
    return re.sub(r'%s$' % regex, '', s)


class CountDataMixin(object):

    def process(self):
        """
        Get processed data from raw
        The default use case involves no processing
        """
        return self.data

    def annotate(self, data, annotate_by=None, annotation_type=None):
        data = data if data is not None else self.data
        annotate_by = annotate_by or self.annotate_by
        annotation_type = annotation_type or self.annotation_type
        if annotate_by == 'Ensembl Gene ID':
            # this is the default indexing anyway
            annotate_by = None
        return annotate(
            data,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            tax_id=self.tax_id
        )

    @property
    def data_by_symbol(self):
        data = self.process()
        return self.annotate(
            data=data,
            annotate_by='Approved Symbol',
            annotation_type='protein_coding',
        )

    @property
    def data_by_ensembl(self):
        data = self.process()
        return self.annotate(
            data=data,
            annotate_by='Ensembl Gene ID',
            annotation_type='protein_coding',
        )

    @property
    def data_by_entrez(self):
        data = self.process()
        return self.annotate(
            data=data,
            annotate_by='Entrez Gene ID',
            annotation_type='protein_coding',
        )

    def get_counts(self):
        return self.data

    def get_fpkm(self):
        # if we already have FPKM computed, use these
        if getattr(self, '_fpkm', None) is not None:
            return self._fpkm
        gene_lengths = normalisation.gene_length_by_tax_id(self.tax_id)
        dat = self.data_by_ensembl
        # discard rows NOT in the gene lengths library
        # these are typically just those 4 rows added by the counter, e.g. 'N_unmapped'
        to_discard = dat.index.difference(gene_lengths.index)
        if len(to_discard):
            logger.warn("Discarding %d rows that were not found in the gene lengths library.", len(to_discard))
            dat = dat.drop(to_discard, axis=0)
        N = dat.sum(axis=0)
        self._fpkm = dat.divide(gene_lengths.iloc[:, 0], axis=0).divide(N, axis=1) * 1e9
        return self._fpkm

    def get_tpm(self):
        # if we already have TPM computed, use these
        if getattr(self, '_tpm', None) is not None:
            return self._tpm
        gene_lengths = normalisation.gene_length_by_tax_id(self.tax_id)
        dat = self.data_by_ensembl
        # discard rows NOT in the gene lengths library
        # these are typically just those 4 rows added by the counter, e.g. 'N_unmapped'
        to_discard = dat.index.difference(gene_lengths.index)
        if len(to_discard):
            logger.warn("Discarding %d rows that were not found in the gene lengths library.", len(to_discard))
            dat = dat.drop(to_discard, axis=0)
        rpb = dat.divide(gene_lengths.iloc[:, 0], axis=0)
        N = rpb.sum()
        self._tpm = rpb / N * 1e6
        return self._tpm

    def get_normed(self):
        if self.annotate_by == 'all':
            raise AttributeError("Cannot normalise data if annotation columns are present.")
        k = self.data.sum(axis=0)
        return self.data.divide(k, axis=1)


class CountDatasetLoader(CountDataMixin):
    cols_to_drop = tuple()
    default_file_ext = ''

    def __init__(
            self,
            meta_fn=None,
            samples=None,
            annotate_by=None,
            annotation_type='protein_coding',
            file_ext=None,
            tax_id=9606,
            *args,
            **kwargs):

        """
        Base class for loading a single dataset.
        :param annotate_by: If supplied, convert the index (initially Ensembl ID) to the requested annotation.
        If 'all' add all supported annotations (will result in additional columns).
        If None, add no extra annotations.
        :param annotation_type: Passed on to the `type` variable of the conversion table loader
        """

        self.meta_fn = meta_fn
        self.samples_to_keep = samples
        self.annotate_by = annotate_by
        self.annotation_type = annotation_type
        self.tax_id = tax_id

        self.logger = get_console_logger(self.__class__.__name__)

        self.raw_meta = None
        self.meta = None
        self.meta_has_sample_names = None
        self.raw_data = None
        self.data = None
        self.transcript_lengths = None
        self.num_reads = None

        # file extension appears in the column headings of the data but can be stripped for brevity
        self.file_ext = file_ext if file_ext is not None else self.__class__.default_file_ext

        self.load_meta()
        self.load_data()
        self.data = self.process()
        self.data = self.annotate(data=self.data)

    def load_meta(self):
        """
        Load the metadata, if a file has been specified
        :return:
        """
        if self.meta_fn is not None:
            self.raw_meta = pd.read_csv(self.meta_fn, header=0, index_col=0)
            if 'sample' not in self.raw_meta.columns:
                self.logger.warning("The meta data has no column 'sample', so sample names will not be used.")
                self.meta_has_sample_names = False
            else:
                self.meta_has_sample_names = True

    def load_data(self):
        """
        Load raw data
        :return:
        """
        raise NotImplementedError

    def process(self):
        """
        Performed after raw data is loaded
        :return: Processed data
        """
        # filter by sample if requested
        if self.samples_to_keep is not None:
            # we require the 'sample' column in metadata
            if not self.meta_has_sample_names:
                raise AttributeError("Sample names have been supplied but the meta data has no 'sample' column.")
            keep_idx = self.raw_meta.loc[:, 'sample'].isin(self.samples_to_keep)
            to_keep = self.raw_meta.loc[keep_idx].index
            # append file_ext if required
            to_keep = [t + self.file_ext for t in to_keep]

        else:
            # just drop the known columns if there are any
            keep_idx = self.raw_meta.index if self.raw_meta is not None else None
            to_keep = self.raw_data.columns.difference(self.cols_to_drop)

        if self.raw_meta is not None:
            self.meta = self.raw_meta.loc[keep_idx]
        data = self.raw_data.loc[:, to_keep]

        # set column names if we have metadata to do so
        # the file extension may already have been stripped by this point (e.g. HTSeqCountLoader)
        # but this will have no effect here
        col_idx = [strip_extension(t, self.file_ext) for t in data.columns]
        if self.meta_has_sample_names:
            data.columns = self.meta.loc[col_idx, 'sample'].values
        else:
            # just strip the extension where present
            data.columns = col_idx
        return data


class FeatureCountLoader(CountDatasetLoader):
    default_file_ext = '.bam'
    cols_to_drop = ['Chr', 'Start', 'End', 'Strand', 'Length']

    def __init__(self, count_file, *args, **kwargs):
        """
        If meta_fn is supplied, it is expected that the first column will match the columns in the data file.
        If samples is supplied, meta is required and the entries of the samples iterable should reference a column
        named `sample`
        :param count_file:
        :param file_ext: If supplied, this is appended to each meta index. In practice, featureCounts s
        :param args:
        :param kwargs:
        """
        self.data_files = count_file
        if kwargs.get('samples') is not None and 'meta_fn' not in kwargs:
            raise AttributeError("If samples are specified, must provide a meta filename.")
        super(FeatureCountLoader, self).__init__(*args, **kwargs)

    def load_data(self):
        self.raw_data = pd.read_csv(self.data_files, comment='#', header=0, index_col=0, sep='\t')

    def process(self):
        # we need to store a useful column now from the raw data
        self.transcript_lengths = self.raw_data.loc[:, 'Length']
        return super(FeatureCountLoader, self).process()

    ## TODO: can we delete this and use the parent method?
    def get_fpkm(self):
        if 'read_count' not in self.meta.columns:
            raise AttributeError("Cannot convert to FPKM without 'read_count' column in metadata.")
        nreads = self.meta.loc[:, 'read_count']
        return self.data.divide(nreads, axis=1).divide(self.transcript_lengths, axis=0) * 1e9

    ## TODO: can we delete this and use the parent method?
    def get_tpm(self):
        rpk = self.data.divide(self.transcript_lengths, axis=0)
        return rpk.divide(rpk.sum(axis=0), axis=1) * 1e6


class MultipleFileCountLoader(CountDatasetLoader):
    """
    Used for loading data when one sample corresponds to one file, as is the case for htseq-count and STAR.
    """
    def __init__(self, count_dir=None, count_files=None, *args, **kwargs):
        if (count_dir is None) == (count_files is None):
            raise AttributeError("Supply EITHER count_files OR count_dir")

        if count_dir is not None:
            # pick out all files in the directory
            file_ext = kwargs.get('file_ext', self.__class__.default_file_ext)
            if not os.path.exists(count_dir):
                raise ValueError("Supplied count_dir does not exist")
            flist = glob.glob(os.path.join(count_dir, '*%s' % file_ext))
            self.data_files = []
            for f in flist:
                if not os.path.isdir(f):
                    if os.path.getsize(f) > 0:
                        self.data_files.append(f)
                    else:
                        logger.warning(
                            "File %s should be included but has size 0. Did something go wrong when you created it?",
                            f
                        )
            if len(self.data_files) == 0:
                raise AttributeError("No files matching the file extension '%s' were found in %s" % (
                    file_ext,
                    count_dir
                ))
        else:
            if not hasattr(self.data_files, '__iter__'):
                self.data_files = [self.data_files]

        super(MultipleFileCountLoader, self).__init__(*args, **kwargs)


    def get_sample_names(self):
        """
        Generate an array of sample names for use as the columns of the counts matrix.
        We don't impose metadata sample names at this stage - that happens in the process() call.
        :return: Iterable of sample names
        """
        return [os.path.basename(t) for t in self.data_files]


class HTSeqCountLoader(MultipleFileCountLoader):
    default_file_ext = '.count'

    def load_data(self):
        sample_names = self.get_sample_names()

        self.raw_data = pd.DataFrame()
        for sn, fn in zip(sample_names, self.data_files):
            # logger.info("Loading file %s with sample name %s", fn, sn)
            self.raw_data.loc[:, sn] = pd.read_csv(fn, sep='\t', index_col=0, header=None).iloc[:, 0]


class StarCountLoader(MultipleFileCountLoader):
    default_file_ext = 'ReadsPerGene.out.tab'

    def __init__(self, strandedness='u', *args, **kwargs):
        if strandedness not in ('u', 'f', 'r'):
            raise ValueError("Variable 'strandedness' must be one of ('u', 'f', 'r')")
        self.strandedness = strandedness
        super(StarCountLoader, self).__init__(*args, **kwargs)

    def load_data(self):
        sample_names = self.get_sample_names()
        self.raw_data = pd.DataFrame()
        for sn, fn in zip(sample_names, self.data_files):
            dat = pd.read_csv(fn, sep='\t', index_col=0, header=None)
            # get correct strand
            if self.strandedness == 'u':
                self.raw_data.loc[:, sn] = dat.iloc[:, 0]
            elif self.strandedness == 'f':
                self.raw_data.loc[:, sn] = dat.iloc[:, 1]
            else:
                self.raw_data.loc[:, sn] = dat.iloc[:, 2]


class MultipleLaneCountDatasetLoader(CountDataMixin):
    loader_class = CountDatasetLoader

    def __init__(
            self,
            meta_fns=None,
            annotate_by='all',
            annotation_type='protein_coding',
            tax_id=9606,
            *args,
            **kwargs
    ):
        self.meta_fns = meta_fns
        self.loaders = None
        self.annotate_by = annotate_by
        self.annotation_type = annotation_type
        self.tax_id = tax_id
        self.load_lanes(annotate_by=annotate_by, annotation_type=annotation_type, tax_id=self.tax_id, *args, **kwargs)
        self.data = self.combine_counts()
        self.meta = self.combine_metas()

    def load_lanes(self, *args, **kwargs):
        raise NotImplementedError

    def combine_counts(self):
        ## FIXME: this results in us concatenating the gene symbol and summing the Entrez ID?
        data = self.loaders[0].data.copy()
        for l in self.loaders[1:]:
            data += l.data
        return data

    def combine_metas(self):
        meta = self.loaders[0].meta.copy()
        idx = meta.index
        b_sample = 'sample' in meta.columns
        if not b_sample:
            logger.warning("Sample names column not found so matching will assume that rows are matched.")
        b_readcount = 'read_count' in meta.columns

        for i, l in enumerate(self.loaders[1:]):
            this_meta = l.meta.copy()
            if b_sample:
                s = meta.loc[idx, 'sample']
                this_meta.set_index('sample', inplace=True)
                # check that entries match
                if len(this_meta.index.intersection(s)) != len(s):
                    raise ValueError("Meta file %s entries do not match meta file %s entries." % (
                        self.meta_fns[0],
                        self.meta_fns[i + 1],
                    ))
                # increment read count
                if b_readcount:
                    meta.loc[idx, 'read_count'] += this_meta.loc[s, 'read_count'].values

            else:
                s = meta.index
                # check that the number of entries match
                if len(this_meta.index) != len(s):
                    raise ValueError("Meta file %s entries do not match meta file %s entries." % (
                        self.meta_fns[0],
                        self.meta_fns[i + 1],
                    ))
                # increment read count
                if b_readcount:
                    meta.loc[:, 'read_count'] += this_meta.loc[:, 'read_count'].values

        return meta


    @property
    def data_by_symbol(self):
        data = self.loaders[0].data_by_symbol
        for l in self.loaders[1:]:
            data += l.data_by_symbol
        return data

    @property
    def data_by_ensembl(self):
        data = self.loaders[0].data_by_ensembl
        for l in self.loaders[1:]:
            data += l.data_by_ensembl
        return data

    @property
    def data_by_entrez(self):
        data = self.loaders[0].data_by_entrez
        for l in self.loaders[1:]:
            data += l.data_by_entrez
        return data

    def get_counts(self):
        return self.data

    # alternative approach to getting FPKM is to pass on to constituent loaders:

    # def get_fpkm(self):
    #     data = self.loaders[0].get_fpkm()
    #     for l in self.loaders[1:]:
    #         data += l.get_fpkm()
    #     return data
    #
    # def get_tpm(self):
    #     data = self.loaders[0].get_tpm()
    #     for l in self.loaders[1:]:
    #         data += l.get_tpm()
    #     # need to renormalise this now (sum to 1e6)
    #     return data


class MultipleLaneFeatureCountLoader(MultipleLaneCountDatasetLoader):
    loader_class = FeatureCountLoader

    def __init__(self, count_files, meta_fns=None, *args, **kwargs):
        if not hasattr(count_files, '__iter__'):
            raise ValueError("count_files must be iterable, otherwise you should use the single lane loader %s" %
                             self.loader_class.__name__)

        if meta_fns is None:
            meta_fns = [None] * len(count_files)

        if len(meta_fns) != len(count_files):
            raise ValueError("meta_fns and count_files must have the same length")
        self.data_files = count_files
        super(MultipleLaneFeatureCountLoader, self).__init__(meta_fns=meta_fns, *args, **kwargs)

    def load_lanes(self, *args, **kwargs):
        self.loaders = []
        for mfn, dfn in zip(self.meta_fns, self.data_files):
            self.loaders.append(self.loader_class(count_file=dfn, meta_fn=mfn, *args, **kwargs))


class MultipleLaneStarCountLoader(MultipleLaneCountDatasetLoader):
    loader_class = StarCountLoader

    def __init__(self, count_files=None, count_dirs=None, meta_fns=None, *args, **kwargs):
        """
        Load and combine count data from multiple lanes using STAR outputs.
        NB: exactly one of count_files and count_dirs must be supplied.
        :param count_files: If supplied, this is a list of lists (or other iterables). Each sublist is the files
        to use for one lane. This is a bit of a masochistic option, directories is much easier.
        :param count_dirs: If supplied, this is an iterable of directories which are searched for files matching the
        pattern ReadsPerGene.out.tab.
        :param meta_fns:
        :param args:
        :param kwargs:
        """
        if (count_dirs is None) == (count_files is None):
            raise AttributeError("Supply EITHER count_files OR count_dirs")

        if count_files is not None:
            self.data_files = count_files
            self.dirs_supplied = False
            if not hasattr(count_files, '__iter__'):
                raise ValueError("count_files must be iterable, otherwise you should use the single lane loader %s" %
                                 self.loader_class.__name__)

        if count_dirs is not None:
            self.data_files = count_dirs
            self.dirs_supplied = True
            if not hasattr(count_dirs, '__iter__'):
                raise ValueError("count_dirs must be iterable, otherwise you should use the single lane loader %s" %
                                 self.loader_class.__name__)

        if meta_fns is None:
            meta_fns = [None] * len(self.data_files)

        if len(meta_fns) != len(self.data_files):
            raise ValueError("meta_fns and count_files must have the same length")

        super(MultipleLaneStarCountLoader, self).__init__(meta_fns=meta_fns, *args, **kwargs)


    def load_lanes(self, *args, **kwargs):
        self.loaders = []
        for mfn, dfn in zip(self.meta_fns, self.data_files):
            # kwarg used depends on whether directories or files were supplied
            if self.dirs_supplied:
                self.loaders.append(self.loader_class(count_dir=dfn, meta_fn=mfn, *args, **kwargs))
            else:
                self.loaders.append(self.loader_class(count_file=dfn, meta_fn=mfn, *args, **kwargs))


class MultipleBatchLoader(CountDataMixin):
    def __init__(self, loaders, intersection_only=True):
        """
        Class to combine multiple loader objects.
        Each loader represents a separate batch. Inputs can include multiple lane loaders.
        :param loaders: Iterable of loader objects.
        :param intersection_only: If True (default), reduce counts to the indices (genes) that are present in all loaders.
        """
        if len(loaders) < 2:
            raise ValueError("Must supply 2 or more loaders to use a MultipleBatchLoader.")
        idx = loaders[0].data.index
        samples = loaders[0].meta.loc[:, 'sample'].tolist()
        meta_cols = set(loaders[0].meta.columns)

        # we may need to append a number to sample names
        sample_appendix = 1
        sample_names_modified = [list(samples)]
        sample_names_original = [list(samples)]

        for l in loaders[1:]:
            the_samples = l.meta.loc[:, 'sample'].values.tolist()
            sample_names_original.append(list(the_samples))

            # check for sample name clash
            first_warn = True
            renamed = False
            while any([t in samples for t in the_samples]):
                renamed = True
                if first_warn:
                    logger.warning("Found sample name clash. Modifying names to avoid errors.")
                    first_warn = False
                the_samples = ['%s_%d' % (t, sample_appendix) for t in l.meta.loc[:, 'sample'].values.tolist()]
                sample_appendix += 1

            # now we know the sample names are unique, append to the list
            samples.extend(the_samples)
            meta_cols.update(l.meta.columns)

            if intersection_only:
                idx = idx.intersection(l.data.index)
            else:
                idx = idx.union(l.data.index)

            sample_names_modified.append(the_samples)

        # we need to add a batch column to meta
        # however, we should first check it doesn't clash
        batch_colname = 'batch'
        if batch_colname in meta_cols:
            i = 0
            while True:
                batch_colname = 'batch_%d' % i
                i += 1
                if batch_colname in meta_cols:
                    continue
                else:
                    break
        self.batch_column = batch_colname
        meta_cols = list(meta_cols) + [batch_colname]

        self.data = pd.DataFrame(index=idx, columns=samples)
        self.meta = pd.DataFrame(index=samples, columns=meta_cols)
        self.tax_id = None

        add_anno = False
        annot_fields = set()
        for i, l in enumerate(loaders):
            this_samples_original = sample_names_original[i]
            this_samples_modified = sample_names_modified[i]
            this_index = l.data.index
            this_meta_cols = l.meta.columns

            # set taxonomy ID and check consistency between loaders
            if self.tax_id is None:
                self.tax_id = l.tax_id
            else:
                if self.tax_id != l.tax_id:
                    raise NotImplementedError("Taxonomy IDs differ between loaders: %s, %s. They must match." % (
                        str(self.tax_id),
                        str(l.tax_id)
                    ))

            # check whether annotations are included
            if l.annotate_by == 'all':
                add_anno = True
                this_annot_fields = l.data.drop(this_samples_original, axis=1).columns.intersection(INDEX_FIELDS)
                annot_fields.update(this_annot_fields)

            self.data.loc[this_index, this_samples_modified] = l.data.loc[this_index, this_samples_original].values
            self.meta.loc[this_samples_modified, this_meta_cols] = l.meta.values
            # add batch column data
            self.meta.loc[this_samples_modified, self.batch_column] = i + 1

        if add_anno:
            # maintain a record of the annotation status of this object, to allow further batching
            self.annotate_by = 'all'
            ann_data = pd.DataFrame(index=idx, columns=list(annot_fields))
            # add annotation columns back in
            for i, l in enumerate(loaders):
                # only looking up the indices and columns in this loader
                this_index = l.data.index
                this_annot_fields = l.data.columns.intersection(annot_fields)

                if len(this_annot_fields) == 0:
                    logger.error("No annotation data in loader %d. Skipping.")
                    continue

                try:
                    ann_data.loc[this_index, this_annot_fields] = l.data.loc[this_index, this_annot_fields]
                except KeyError:
                    logger.error("No annotation data in loader %d. Skipping.")
            self.data = pd.concat((self.data, ann_data), axis=1)
        else:
            abs = set([l.annotate_by for l in loaders])
            if len(abs) > 1:
                logger.error("The annotate_by variable differs between loaders: %s", ', '.join(list(abs)))
                logger.error("Picking one arbitrarily, but expect this to be broken!")
            self.annotate_by = loaders[0].annotate_by


def annotate(
    res,
    annotate_by='all',
    annotation_type='protein_coding',
    tax_id=9606,
):
    """
    Annotate the supplied dataframe
    :param res: The dataframe to annotate
    :param annotate_by: If supplied, convert the index (initially Ensembl ID) to the requested annotation.
    If 'all' add all supported annotations.
    If None, add no extra annotations.
    :param annotation_type: Passed on to the `type` variable of the conversion table loader
    :param tax_id: Default to human (9606). Mouse is 10090.
    """
    if annotate_by is not None:
        # load genenames data for annotation
        df = references.conversion_table(type=annotation_type, tax_id=tax_id)
        df.set_index('Ensembl Gene ID', inplace=True)
        # resolve duplicates (if present) by keeping the first
        df = df.loc[~df.index.duplicated(keep='first')]

        if annotate_by == 'all':
            # depending on the taxonomy and type, we may have different columns
            # therefore, if we are including all annotations, we should only use those that actually exist
            all_cols = [
                'Approved Symbol', 'Entrez Gene ID', 'RefSeq IDs'
            ]
            cols_to_use = df.columns.intersection(all_cols)
            annot = df.loc[res.index.intersection(df.index), cols_to_use]
        else:
            annot = df.loc[res.index.intersection(df.index), annotate_by]

        # add annotation columns
        res = pd.concat((res, annot), axis=1)

        # if one annotation was requested, set that as the index
        if annotate_by != 'all':
            # drop any rows that do not have an annotation
            res.dropna(axis=0, subset=[annotate_by], inplace=True)
            # set index
            res.set_index(annotate_by, inplace=True)
    return res


def htseqcounts(
        count_files,
        sample_names=None,
        metafile=None,
        units='counts',
        annotate_by='all',
        annotation_type='protein_coding',
):
    """
    Generic loader for htseq-count data
    :param count_files: Iterable containing paths to the count files
    :param sample_names: Iterable of same length as count_files containing the sample names. If missing, the
    filename is used.
    :param metafile: Path to the metadata; optional
    :param units: One of 'counts', 'fpkm', 'tpm'
    :param annotate_by: If supplied, convert the index (initially Ensembl ID) to the requested annotation.
    If 'all' add all supported annotations.
    If None, add no extra annotations.
    :param annotation_type: Passed on to the `type` variable of the conversion table loader
    :return:
    """
    ## TODO: units other than counts not supported at present
    if units not in ('counts',):
        raise ValueError("Unsupported units requested.")

    if sample_names is not None:
        if len(sample_names) != len(count_files):
            raise ValueError("Length of sample_names does not equal the length of count_files")
    else:
        sample_names = [os.path.basename(t) for t in count_files]

    dat = pd.DataFrame()
    for sn, fn in zip(sample_names, count_files):
        t = pd.read_csv(fn, sep='\t', index_col=0, header=None).iloc[:, 0]
        dat.loc[:, sn] = t

    dat = annotate(dat, annotate_by=annotate_by, annotation_type=annotation_type)

    if metafile is not None:
        meta = pd.read_csv(metafile, header=0, index_col=0)
        return dat, meta
    else:
        return dat


def featurecounts(
        count_files,
        metafiles,
        samples=None,
        units='counts',
        annotate_by='all',
        annotation_type='protein_coding'):
    """
    Iterate over one or more lanes of gene count data computed using featureCounts, summing data over the lanes.

    :param count_files: Iterable. Each entry is the path to a featureCounts output.
    :param metafiles: Iterable, same length as count_files. Each entry is the path to a CSV file. Each row details a
    sample. The column 'sample' gives the identifier. The first (index) column is the file stem in the counts file.
    :param samples: If supplied, these are the samples to keep. Everything else is dropped.
    :param units: One of 'counts', 'fpkm', 'tpm'
    :param annotate_by: If supplied, convert the index (initially Ensembl ID) to the requested annotation.
    If 'all' add all supported annotations.
    If None, add no extra annotations.
    :param annotation_type: Passed on to the `type` variable of the conversion table loader
    :return:
    """
    supported_annot = ('Approved Symbol', 'Entrez Gene ID', 'RefSeq IDs', 'all')
    if annotate_by is not None and annotate_by not in supported_annot:
        raise ValueError("Unrecognised annotation requested. Supported options are %s" % ', '.join(supported_annot))

    supported_units = ('counts', 'fpkm', 'tpm')
    if units not in supported_units:
        raise ValueError("Unrecognised units requested. Supported options are %s" % ', '.join(supported_units))

    first_run = True
    res = None
    lengths = None
    nreads = None
    meta = None
    for fn, mfn in zip(count_files, metafiles):
        meta = pd.read_csv(mfn, header=0, index_col=0)
        dat = pd.read_csv(fn, comment='#', header=0, index_col=0, sep='\t')

        # keep only counts for each requested sample
        if samples is not None:
            codes_to_keep = meta.loc[meta.loc[:, 'sample'].isin(samples)].index
        else:
            # keep all samples
            codes_to_keep = meta.index

        files_to_keep = ['%s.bam' % t for t in codes_to_keep]
        sample_names = meta.loc[codes_to_keep, 'sample'].values
        lengths = dat.Length

        dat = dat.loc[:, files_to_keep]
        dat.columns = sample_names

        if first_run:
            res = dat.copy()
            # get total reads and length of each transcript
            nreads = meta.loc[codes_to_keep, 'read_count']
            nreads.index = sample_names
            first_run = False
        else:
            # sum counts and total count tally
            res += dat
            nr = meta.loc[codes_to_keep, 'read_count']
            nr.index = sample_names
            nreads += nr

    # normalise (if required) AFTER combining all count files
    # we are assuming that the transcripts do not vary between files
    if units == 'fpkm':
        res = res.divide(nreads, axis=1).divide(lengths, axis=0) * 1e9
    elif units == 'tpm':
        rpk = res.divide(lengths, axis=0)
        res = rpk.divide(rpk.sum(axis=0), axis=1) * 1e6

    res = annotate(res, annotate_by=annotate_by, annotation_type=annotation_type)
    return res, meta


def gse83696(index_by='Ensembl Gene ID'):
    """
    Data are initially indexed by Ensembl ID. Coversion is carried out using HGNC data, if required.
    Index field options are: Approved Symbol, Entrez Gene ID, RefSeq IDs
    :param index_by:
    :return:
    """
    # TODO: convert this into a generic loader for htseq-count outputs
    indir = os.path.join(GIT_LFS_DATA_DIR, 'rnaseq_GSE83696', 'htseq-count')
    samples = [
        ('XZ1', 'XZ-1.count'),
        ('XZ2', 'XZ-2.count'),
        ('XZ3', 'XZ-3.count'),
        ('XZ4', 'XZ-4.count'),
        ('XZ5', 'XZ-5.count'),
        ('XZ6', 'XZ-6.count'),
        ('XZ7', 'XZ-7.count'),
        ('XZ8', 'XZ-8.count'),
    ]
    df = pd.DataFrame()
    for sn, fn in samples:
        t = pd.read_csv(os.path.join(indir, fn), sep='\t', index_col=0, header=None).iloc[:, 0]
        df.loc[:, sn] = t

    if index_by is not None and index_by != 'Ensembl Gene ID':
        new_idx = references.translate(df.index, to_field=index_by, from_field='Ensembl Gene ID')
        new_idx.dropna(inplace=True)
        df = df.loc[new_idx.index]
        df.index = new_idx.values
    return df


def gbm_paired_samples(units='counts', annotate_by='all', annotation_type='protein_coding'):
    # TODO: obsolete, remove
    obj = gbm_paired_samples_loader(annotate_by=annotate_by, annotation_type=annotation_type)

    if units == 'counts':
        return obj.data
    elif units == 'fpkm':
        return obj.get_fpkm()
    elif units == 'tpm':
        return obj.get_tpm()
    else:
        raise AttributeError("Units not recognised")


def gbm_paired_samples_loader(source='star', annotate_by='all', annotation_type='protein_coding'):
    # TODO: replace or rename?
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
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

    if source == 'star':
        count_dirs = [os.path.join(d, 'star_alignment') for d in (lane1dir, lane2dir)]
        obj = MultipleLaneStarCountLoader(
            count_dirs=count_dirs,
            meta_fns=metafiles,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            samples=samples,
            strandedness='r',
        )
    elif source == 'htseq-count':
        raise NotImplementedError
    elif source == 'featurecounts':
        count_files = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
        obj = MultipleLaneFeatureCountLoader(
            count_files=count_files,
            meta_fns=metafiles,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            samples=samples,
        )
    else:
        raise ValueError("Unrecognised source.")
    return obj


def gbm_astrocyte_nsc_samples_loader(source='star', annotate_by='all', annotation_type='protein_coding'):
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        'DURA018N2_NSC',
        'DURA019N8C_NSC',
        'DURA018N2_ASTRO_DAY12',
        'DURA019N8C_ASTRO_DAY12',
    )

    if source == 'star':
        count_dirs = [os.path.join(d, 'star_alignment') for d in (lane1dir, lane2dir)]
        obj = MultipleLaneStarCountLoader(
            count_dirs=count_dirs,
            meta_fns=metafiles,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            samples=samples,
            strandedness='r',
        )
    elif source == 'htseq-count':
        raise NotImplementedError
    elif source == 'featurecounts':
        count_files = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
        obj = MultipleLaneFeatureCountLoader(
            count_files=count_files,
            meta_fns=metafiles,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            samples=samples,
        )
    else:
        raise ValueError("Unrecognised source.")
    return obj


def all_samples_multilane_loader(
        countdirs,
        metafiles,
        samples=None,
        source='star',
        annotate_by='all',
        annotation_type='protein_coding',
        strandedness='r',
        tax_id=9606,
):
    """
    Load the full dataset from a multiple lane run
    :param lanedirs: The paths to the fastq files. It is expected that a subdirectory called star_alignment etc. exists
    here.
    :param metafiles: Paths to the metafiles
    :param source:
    :param annotate_by:
    :param annotation_type:
    :return:
    """
    if source == 'star':
        obj = MultipleLaneStarCountLoader(
            count_dirs=countdirs,
            meta_fns=metafiles,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            strandedness=strandedness,
            samples=samples,
            tax_id=tax_id
        )
    elif source == 'htseq-count':
        raise NotImplementedError
    elif source == 'featurecounts':
        count_files = [os.path.join(d, 'counts.txt') for d in countdirs]
        obj = MultipleLaneFeatureCountLoader(
            count_files=count_files,
            meta_fns=metafiles,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            samples=samples,
            tax_id=tax_id
        )
    else:
        raise ValueError("Unrecognised source.")
    return obj


def all_hgic_loader(source='star', annotate_by='all', annotation_type='protein_coding', include_derived=False):
    loaders = []

    # expt 1
    samples = (
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
    if include_derived:
        samples += ('DURA018_ASTRO_N2_DAY12', 'DURA019_ASTRO_N8C_DAY12')

    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lanedirs = [
        os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX'),
        os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'star_alignment') for d in lanedirs]
    loaders.append(
        all_samples_multilane_loader(
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
            samples=samples, strandedness='r'
        )
    )

    # expt 2
    samples = (
        'GBM044_P4',
        'GBM026_P3n4',
        'GBM018_P12',
        'GBM044_P8',
        'DURA044_NSC_N8_P2',
        'GIBCO_NSC_P4',
        'DURA044_NSC_N17_P3',
        'DURA018_NSC_N4_P4',
        'GBM030_P5',
        'DURA030_NSC_N16B6_P1',
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170218')
    lanedirs = [
        os.path.join(indir, '170509_K00150_0192_BHJKCLBBXX'),
        os.path.join(indir, '170515_K00150_0196_BHJKC5BBXX_lane_2'),
        os.path.join(indir, '170515_K00150_0196_BHJKC5BBXX_lane_3')
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'human', 'star_alignment') for d in lanedirs]
    loaders.append(
        all_samples_multilane_loader(
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
            samples=samples, strandedness='r'
        )
    )

    # expt 3
    samples = (
        'GBM030_P9n10',
        'GBM031_P7',
        'GBM019_P3n6',
        'GBM017_P3',
        'GBM017_P4',
        'DURA017_NSC_P4_N3C5',
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170390')
    lanedirs = [
        os.path.join(indir, '170727_K00198_0222_AHKWW5BBXX'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_1'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_2')
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'human', 'star_alignment') for d in lanedirs]
    loaders.append(
        all_samples_multilane_loader(
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
            samples=samples, strandedness='r'
        )
    )

    # expt 4
    samples = (
        'GBM049_P4',
        'GBM049_P6',
        'GBM050_P7',
        'GBM050_P9',
        'GBM054_P4',
        'GBM054_P6',
        'GBM061_P3',
        'GBM061_P5',
        'GBM052_P6n7',
        'GBM052_P4n5',
        'DURA049_NSC_N19_P4',
        'DURA049_NSC_N5_P2',
        'DURA050_NSC_N12_P3',
        'DURA050_NSC_N16_P4',
        'DURA054_NSC_N3C_P2',
        'DURA054_NSC_N2E_P1',
        'DURA061_N4_P2',
        'DURA061_NSC_N6_P4',
        'DURA052_NSC_N4_P3',
        'DURA052_NSC_N5_P2',
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170503')
    lanedirs = [
        os.path.join(indir, '170929_K00150_0250_BHLGNHBBXX'),
        os.path.join(indir, '171003_K00198_0242_AHLGYVBBXX_1'),
        os.path.join(indir, '171003_K00198_0242_AHLGYVBBXX_2')
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'star_alignment') for d in lanedirs]
    loaders.append(
        all_samples_multilane_loader(
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
            samples=samples, strandedness='r'
        )
    )

    return MultipleBatchLoader(loaders)


def rtki_hgic_loader(source='star', annotate_by='all', annotation_type='protein_coding', include_reference=True):
    """
    GBM018, GBM019, GBM030, GBM031
    Includes GBM cultures, iNSC
    :param source:
    :param annotate_by:
    :param annotation_type:
    :param include_reference: if True (default), also include the Gibco reference sample
    :return:
    """

    loaders = []

    # expt 1
    samples = (
        'GBM018_P10',
        'GBM019_P4',
        'GBM031_P4',
        'DURA018_NSC_N2_P6',
        'DURA019_NSC_N8C_P2',
        'DURA031_NSC_N44B_P2',
    )

    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lanedirs = [
        os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX'),
        os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'star_alignment') for d in lanedirs]
    loaders.append(
        all_samples_multilane_loader(
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
            samples=samples, strandedness='r'
        )
    )

    # expt 2
    samples = (
        'GBM018_P12',
        'DURA018_NSC_N4_P4',
        'GBM030_P5',
        'DURA030_NSC_N16B6_P1',
    )

    if include_reference:
        samples += ('GIBCO_NSC_P4',)

    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170218')
    lanedirs = [
        os.path.join(indir, '170509_K00150_0192_BHJKCLBBXX'),
        os.path.join(indir, '170515_K00150_0196_BHJKC5BBXX_lane_2'),
        os.path.join(indir, '170515_K00150_0196_BHJKC5BBXX_lane_3')
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'human', 'star_alignment') for d in lanedirs]
    loaders.append(
        all_samples_multilane_loader(
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
            samples=samples, strandedness='r'
        )
    )

    # expt 3
    samples = (
        'GBM030_P9n10',
        'GBM031_P7',
        'GBM019_P3n6',
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170390')
    lanedirs = [
        os.path.join(indir, '170727_K00198_0222_AHKWW5BBXX'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_1'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_2')
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'human', 'star_alignment') for d in lanedirs]
    loaders.append(
        all_samples_multilane_loader(
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
            samples=samples, strandedness='r'
        )
    )

    return MultipleBatchLoader(loaders)


def rtkii_hgic_loader(source='star', annotate_by='all', annotation_type='protein_coding', include_reference=True):
    """
    GBM017, GBM050, GBM054, GBM061
    Includes GBM cultures, iNSC
    :param source:
    :param annotate_by:
    :param annotation_type:
    :param include_reference: if True (default), also include the Gibco reference sample
    :return:
    """
    loaders = []
    # expt 1
    samples = (
        'GBM017_P3',
        'GBM017_P4',
        'DURA017_NSC_P4_N3C5',
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170390')
    lanedirs = [
        os.path.join(indir, '170727_K00198_0222_AHKWW5BBXX'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_1'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_2')
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'human', 'star_alignment') for d in lanedirs]
    loaders.append(
        all_samples_multilane_loader(
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
            samples=samples, strandedness='r'
        )
    )

    # expt 2
    samples = (
        'GBM050_P7',
        'GBM050_P9',
        'GBM054_P4',
        'GBM054_P6',
        'GBM061_P3',
        'GBM061_P5',
        'DURA050_NSC_N12_P3',
        'DURA050_NSC_N16_P4',
        'DURA054_NSC_N3C_P2',
        'DURA054_NSC_N2E_P1',
        'DURA061_N4_P2',
        'DURA061_NSC_N6_P4',
    )
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170503')
    lanedirs = [
        os.path.join(indir, '170929_K00150_0250_BHLGNHBBXX'),
        os.path.join(indir, '171003_K00198_0242_AHLGYVBBXX_1'),
        os.path.join(indir, '171003_K00198_0242_AHLGYVBBXX_2')
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'star_alignment') for d in lanedirs]
    loaders.append(
        all_samples_multilane_loader(
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
            samples=samples, strandedness='r'
        )
    )

    if include_reference:
        samples = ('GIBCO_NSC_P4',)

        indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170218')
        lanedirs = [
            os.path.join(indir, '170509_K00150_0192_BHJKCLBBXX'),
            os.path.join(indir, '170515_K00150_0196_BHJKC5BBXX_lane_2'),
            os.path.join(indir, '170515_K00150_0196_BHJKC5BBXX_lane_3')
        ]
        metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
        countdirs = [os.path.join(d, 'human', 'star_alignment') for d in lanedirs]
        loaders.append(
            all_samples_multilane_loader(
                countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type,
                samples=samples, strandedness='r'
            )
        )

    return MultipleBatchLoader(loaders)


def load_by_patient(
        patient_ids,
        type='cell_culture',
        source='star',
        annotate_by='all',
        annotation_type='protein_coding',
        include_control=True
):
    """
    Load all RNA-Seq count data associated with the patient ID(s) supplied
    :param patient_ids: Iterable or single int or char
    :param source:
    :param annotate_by:
    :param annotation_type:
    :param include_control: If True (default) include Gibco reference NSC
    :return:
    """

    if source == 'star':
        if type == "cell_culture":
            LOOKUP = PATIENT_LOOKUP_CC_STAR
        elif type == "ffpe":
            LOOKUP = PATIENT_LOOKUP_FFPE_STAR
        else:
            raise NotImplementedError()
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
        objs.append(
            MultipleLaneStarCountLoader(
                annotate_by=annotate_by,
                annotation_type=annotation_type,
                samples=samples,
                **ldr.loader_kwargs
            )
        )

    if len(objs) > 1:
        res = MultipleBatchLoader(objs)
    else:
        res = objs[0]
        # make samples column the meta index
        res.meta.set_index('sample', inplace=True)

    # apply original ordering
    res.meta = res.meta.loc[sample_order]
    ## FIXME: maintain all annotation columns if preset
    res.data = res.data.loc[:, res.meta.index]

    return res


def atcc_cell_lines(source='star', annotate_by='all', annotation_type='protein_coding'):
    """
    Two ATCC lines sequenced, but one is mouse contaminated
    :param source:
    :param annotate_by:
    :param annotation_type:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170390')
    lane_dirs = [
        os.path.join(indir, '170727_K00198_0222_AHKWW5BBXX'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_1'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_2'),
    ]
    metafiles = [os.path.join(d, 'sources.csv') for d in lane_dirs]
    samples = (
        'CRL3021_PRIMARYMB_P3',
    )

    if source == 'star':
        count_dirs = [os.path.join(d, 'human', 'star_alignment') for d in lane_dirs]
        obj = MultipleLaneStarCountLoader(
            count_dirs=count_dirs,
            meta_fns=metafiles,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            samples=samples,
            strandedness='r',
        )
    else:
        raise NotImplementedError("Unrecognised source.")
    return obj


def zhao_mb_cultures(source='star', annotate_by='all', annotation_type='protein_coding'):
    """
    NB: ICb1078 and ICb1487 are both mouse contaminated, so the results do not have any real meaning
    :param source:
    :param annotate_by:
    :param annotation_type:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        'ICb1078',
        'ICb1595',
        'ICb1487'
    )

    if source == 'star':
        count_dirs = [os.path.join(d, 'star_alignment') for d in (lane1dir, lane2dir)]
        obj = MultipleLaneStarCountLoader(
            count_dirs=count_dirs,
            meta_fns=metafiles,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            samples=samples,
            strandedness='r',
        )
    elif source == 'htseq-count':
        raise NotImplementedError
    elif source == 'featurecounts':
        count_files = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
        obj = MultipleLaneFeatureCountLoader(
            count_files=count_files,
            meta_fns=metafiles,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            samples=samples,
        )
    else:
        raise ValueError("Unrecognised source.")
    return obj


def brainrnaseq_preprocessed():
    """
    Load the pre-processed results from www.brainrnaseq.org
    These are in units of FPKM, annotated by gene. Annoyingly, the gene symbols are for mouse.
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE73721')
    infile = os.path.join(indir, 'fpkm', 'fpkm.csv')
    meta_fn = os.path.join(indir, 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=0)
    data = pd.read_csv(infile, header=0, index_col=0)

    return data, meta


def nih_gdc_gbm_preprocessed(indir, meta_fn=None):
    """
    Load the preprocessed gene expression data downloaded from NIH genomic data commons
    :param indir: Must be specified as different versions are available.
    :param meta_fn: Optionally supply the path to the Brenna metadata. This should be set anyway, so no need.
    :return: data, meta
    """
    infile = os.path.join(indir, 'data.csv')
    sources_fn = os.path.join(indir, 'sources.csv')

    if not os.path.exists(infile) or not os.path.exists(sources_fn):
        logger.info("Unable to find collated data file %s. Calculating from individual files now.", infile)
        dat, meta = tcga.prepare_data_and_meta(indir, meta_fn=meta_fn)

        logger.info("Saving data from %d samples to %s.", dat.columns.size, infile)
        dat.to_csv(infile)
        logger.info("Saving metadata from to %s.", sources_fn)
        meta.to_csv(sources_fn)
    else:
        dat = pd.read_csv(infile, header=0, index_col=0)
        meta = pd.read_csv(sources_fn, header=0, index_col=0)

    return dat, meta


def tcga_primary_gbm(units='counts'):
    """
    Load the preprocessed gene expression data downloaded from NIH genomic data commons
    :param units: Either 'counts' or 'fpkm' are available.
    :return: data, meta
    """
    basedir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour')
    meta_fn = os.path.join(basedir, 'brennan_s7.csv')
    if units == 'counts':
        indir = os.path.join(basedir, 'htseq-count')
    elif units == 'fpkm':
        indir = os.path.join(basedir, 'htseq-count_fpkm')
    else:
        raise ValueError("Unsupported units. Supported values are 'counts' and 'fpkm'.")

    return nih_gdc_gbm_preprocessed(indir, meta_fn=meta_fn)


def tcga_methylation_assayed(units='counts'):
    """
    Load the preprocessed gene expression data downloaded from NIH genomic data commons
    :param units: Either 'counts' or 'fpkm' are available.
    :return: data, meta
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'paired_rnaseq_methylation_cohort')
    return nih_gdc_gbm_preprocessed(indir, units=units)


def gse73721(source='star', annotate_by='all', annotation_type='protein_coding'):
    """
    Barres data on GEO.
    Astrocytes, oligodendrocytes, etc...
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE73721')
    metafn = os.path.join(indir, 'sources.csv')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            annotation_type=annotation_type
        )
    elif source == 'htseq-count':
        logger.warning("htseq-count fails with some samples in this dataset.")
        obj = HTSeqCountLoader(
            count_dir=os.path.join(indir, 'hisat2_alignment', 'htseq-count'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            annotation_type=annotation_type
        )
    elif source == 'featurecounts':
        obj = FeatureCountLoader(
            count_file=os.path.join(indir, 'hisat2_alignment', 'featureCounts'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            annotation_type=annotation_type
        )
    else:
        raise ValueError("Unrecognised source.")
    return obj


def gse61794(source='star', annotate_by='all', annotation_type='protein_coding', collapse_replicates=False):
    """
    2 samples similar to Gibco NSC line.
    These may be technical replicates as they are very highly correlated.
    :param collapse_replicates: If True, combine the counts to make one sample
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE61794')
    metafn = os.path.join(indir, 'sources.csv')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            strandedness='u',
        )
    elif source == 'htseq-count':
        raise NotImplementedError
        # obj = HTSeqCountLoader(
        #     count_dir=os.path.join(indir, 'hisat2_alignment', 'htseq-count'),
        #     meta_fn=metafn,
        #     annotate_by=annotate_by,
        #     annotation_type=annotation_type
        # )
    elif source == 'featurecounts':
        obj = FeatureCountLoader(
            count_file=os.path.join(indir, 'hisat2_alignment', 'featureCounts'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            annotation_type=annotation_type
        )
    else:
        raise ValueError("Unrecognised source.")

    if collapse_replicates:
        ## TODO?
        raise NotImplementedError("collapse_replicates")

    return obj


def pollard_nsc(source='star', annotate_by='all', annotation_type='protein_coding'):
    """
    2 biological replicates of NSC cells.
    :param source:
    :param annotate_by:
    :param annotation_type:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'E-MTAB-3867')
    metafn = os.path.join(indir, 'sources.csv')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            strandedness='u',
        )
    else:
        raise ValueError("Unrecognised source")

    return obj


def gbm_ribozero_samples_loader(source='star', annotate_by='all', annotation_type='protein_coding'):
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704_ribozero', '170328_K00150_0177_BHJ2C2BBXX')
    metafn = os.path.join(indir, 'sources.csv')
    samples = (
        'GBM018 FFPE',
        'GBM019 FFPE',
        'GBM031 FFPE',
        'GBM026 FF',
        'GBM031 FF',
    )

    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            strandedness='r',
        )
    else:
        raise NotImplementedError

    return obj


def mouse_nsc_validation_samples(source='star', annotate_by='all'):
    """
    3 mice, 4 conditions:
    eNSC in endogenous media
    eNSC in mouse induction media
    iNSC mouse protocol
    iNSC human protocol
    Three of the samples were re-run following a facility mess up (CRL3034 contamination). Use the newest versions.
    :param source:
    :param annotate_by:
    :return:
    """
    loaders = []

    # expt 1: original data
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170390')
    lanedirs = [
        os.path.join(indir, '170727_K00198_0222_AHKWW5BBXX'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_1'),
        os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_2')
    ]

    metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
    countdirs = [os.path.join(d, 'mouse', 'star_alignment') for d in lanedirs]
    samples = ['eNSC%dmouse' % i for i in (3, 5, 6)] \
    + ['mDura%smouse' % i for i in ('3N1', '5N24A', '6N6')] \
    + ['mDura%shuman' % i for i in ('3N1', '5N24A', '6N6')]

    obj = all_samples_multilane_loader(
        countdirs,
        metafiles,
        strandedness='r',
        source=source,
        annotate_by=annotate_by,
        samples=samples,
        tax_id=10090
    )
    loaders.append(obj)

    # expt 2: 3 x replacement runs

    samples = ['eNSC%dmed' % i for i in (3, 5, 6)]
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p170506', '170829_K00150_0236_AHL5YHBBXX')
    countdir = os.path.join(indir, 'mouse', 'star_alignment')
    metafn = os.path.join(indir, 'sources.csv')
    obj = StarCountLoader(
        count_dir=countdir,
        strandedness='r',
        samples=samples,
        meta_fn=metafn,
        source=source,
        annotate_by=annotate_by,
        tax_id=10090
    )
    loaders.append(obj)

    return MultipleBatchLoader(loaders)


def gse77920_loader(source='star', annotate_by='all', annotation_type='protein_coding'):
    """
    Single deep sequenced H9 ESC line
    :param source:
    :param annotate_by:
    :param annotation_type:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE77920')
    metafn = os.path.join(indir, 'sources.csv')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            annotation_type=annotation_type,
            strandedness='r',
        )
    else:
        raise ValueError("Unrecognised source")

    return obj

def gse24399_merged_loader(source='star', annotate_by='all', annotation_type='protein_coding'):
    """
    Merged technical triplicate H9 ESC cell line
    :param source:
    :param annotate_by:
    :param annotation_type:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE24399')
    metafn = os.path.join(indir, 'sources.csv')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            annotation_type=annotation_type
        )
    else:
        raise ValueError("Unrecognised source")

    return obj


def gse52564(source='star', annotate_by='all', **kwargs):
    """
    Mouse reference data. Duplicates of: astrocyte, neuron, OPC, microglia.
    Unstranded, based on looking at the gene counts
    :param source:
    :param annotate_by:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE52564')
    metafn = os.path.join(indir, 'sources.csv')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            strandedness='u',
            tax_id=10090,
            **kwargs
        )
    else:
        raise ValueError("Unrecognised source")

    return obj


def gse43916(source='star', annotate_by='all', **kwargs):
    """
    Mouse reference data. One astrocyte sample (poor quality), one NSC (good quality)
    Reverse stranded single end, based on looking at the gene counts
    :param source:
    :param annotate_by:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE43916')
    metafn = os.path.join(indir, 'sources.csv')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            strandedness='r',
            tax_id=10090,
            **kwargs
        )
    else:
        raise ValueError("Unrecognised source")

    return obj



def gse86248(source='star', annotate_by='all'):
    """
    Mouse reference data. 3 repeats of CL57/BL6 ESC cultured in LIF.
    Reverse stranded single end, based on looking at the gene counts
    :param source:
    :param annotate_by:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86248')
    metafn = os.path.join(indir, 'sources.csv')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            strandedness='r',
            tax_id=10090
        )
    else:
        raise ValueError("Unrecognised source")

    return obj


def gse36114(source='star', annotate_by='all', **kwargs):
    """
    Mouse reference data. E14 commercial ESC line at day 0 and 4 days under Activin A differentiation.
    Unstranded single end, based on looking at the gene counts
    :param source:
    :param annotate_by:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE36114')
    metafn = os.path.join(indir, 'sources.csv')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=os.path.join(indir, 'star_alignment'),
            meta_fn=metafn,
            annotate_by=annotate_by,
            strandedness='u',
            tax_id=10090,
            **kwargs
        )
    else:
        raise ValueError("Unrecognised source")

    return obj


def gse64411(source='star', annotate_by='all', trimmed=False, **kwargs):
    """
    Mouse reference data. ESC, neurons and astrocytes from C57BL/6.
    Unstranded single end, based on looking at the gene counts.
    :param source:
    :param annotate_by:
    :param trimmed: If True, use the data that was trimmed before alignment
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE64411')
    metafn = os.path.join(indir, 'sources.csv')
    if trimmed:
        count_dir = os.path.join(indir, 'trimgalore', 'star_alignment')
    else:
        count_dir = os.path.join(indir, 'star_alignment')
    if source == 'star':
        obj = StarCountLoader(
            count_dir=count_dir,
            meta_fn=metafn,
            annotate_by=annotate_by,
            strandedness='u',
            tax_id=10090,
            **kwargs
        )
    else:
        raise ValueError("Unrecognised source")

    return obj
