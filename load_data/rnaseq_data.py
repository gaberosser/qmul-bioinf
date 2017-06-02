import pandas as pd
import re
import os
import glob
import references
from utils.log import get_console_logger
from settings import DATA_DIR, DATA_DIR_NON_GIT
logger = get_console_logger(__name__)

INDEX_FIELDS = (
    'Approved Symbol',
    'Entrez Gene ID',
    'RefSeq IDs',
    'Ensembl Gene ID'
)


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
            annotation_type=annotation_type
        )

    @property
    def data_by_symbol(self):
        data = self.process()
        return self.annotate(data=data, annotate_by='Approved Symbol', annotation_type='protein_coding')

    @property
    def data_by_ensembl(self):
        data = self.process()
        return self.annotate(data=data, annotate_by='Ensembl Gene ID', annotation_type='protein_coding')

    @property
    def data_by_entrez(self):
        data = self.process()
        return self.annotate(data=data, annotate_by='Entrez Gene ID', annotation_type='protein_coding')

    def get_counts(self):
        return self.data

    def get_fpkm(self):
        raise NotImplementedError

    def get_tpm(self):
        raise NotImplementedError

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

    def get_fpkm(self):
        if 'read_count' not in self.meta.columns:
            raise AttributeError("Cannot convert to FPKM without 'read_count' column in metadata.")
        nreads = self.meta.loc[:, 'read_count']
        return self.data.divide(nreads, axis=1).divide(self.transcript_lengths, axis=0) * 1e9

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


class MultipleLaneCountDatasetLoader(object):
    loader_class = CountDatasetLoader

    def __init__(
            self,
            meta_fns=None,
            annotate_by='all',
            annotation_type='protein_coding',
            *args,
            **kwargs
    ):
        self.meta_fns = meta_fns
        self.loaders = None
        self.annotate_by = annotate_by
        self.annotation_type = annotation_type
        self.load_lanes(annotate_by=annotate_by, annotation_type=annotation_type, *args, **kwargs)
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

    def get_fpkm(self):
        raise NotImplementedError

    def get_tpm(self):
        raise NotImplementedError


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

        for l in loaders[1:]:
            samples.extend(l.meta.loc[:, 'sample'].values)
            meta_cols.update(l.meta.columns)

            if intersection_only:
                idx = idx.intersection(l.data.index)
            else:
                idx = idx.union(l.data.index)

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

        add_anno = False
        if any([l.annotate_by == 'all' for l in loaders]):
            add_anno = True
            ann_fields = [t for t in INDEX_FIELDS if t != 'Ensembl Gene ID']
            ann_data = pd.DataFrame(index=idx, columns=ann_fields)

        self.data = pd.DataFrame(index=idx, columns=samples)
        self.meta = pd.DataFrame(index=samples, columns=meta_cols)
        for i, l in enumerate(loaders):
            this_samples = l.meta.loc[:, 'sample'].values
            this_index = l.data.index
            this_meta_cols = l.meta.columns

            self.data.loc[this_index, this_samples] = l.data.loc[this_index, this_samples].values
            self.meta.loc[this_samples, this_meta_cols] = l.meta.values
            # add batch column data
            self.meta.loc[this_samples, self.batch_column] = i + 1

        if add_anno:
            # add annotation columns back in
            for i, l in enumerate(loaders):
                this_index = l.data.index
                try:
                    ann_data.loc[this_index, ann_fields] = l.data.loc[this_index, ann_fields]
                except KeyError:
                    logger.error("No annotation data in loader %d. Skipping.")
            self.data = pd.concat((self.data, ann_data), axis=1)


def annotate(
    res,
    annotate_by='all',
    annotation_type='protein_coding',
):
    """
    Annotate the supplied dataframe
    :param res: The dataframe to annotate
    :param annotate_by: If supplied, convert the index (initially Ensembl ID) to the requested annotation.
    If 'all' add all supported annotations.
    If None, add no extra annotations.
    :param annotation_type: Passed on to the `type` variable of the conversion table loader
    """
    if annotate_by is not None:
        # load genenames data for annotation
        df = references.conversion_table(type=annotation_type)
        df.set_index('Ensembl Gene ID', inplace=True)

        if annotate_by == 'all':
            annot = df.loc[res.index.intersection(df.index), ['Approved Symbol', 'Entrez Gene ID', 'RefSeq IDs']]
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
    indir = os.path.join(DATA_DIR, 'rnaseq_GSE83696', 'htseq-count')
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
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    count_files = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        'GBM018',
        'GBM019',
        'GBM024',
        'GBM026',
        'GBM031',
        'DURA018N2_NSC',
        'DURA019N8C_NSC',
        'DURA024N28_NSC',
        'DURA026N31D_NSC',
        'DURA031N44B_NSC',
    )

    obj = MultipleLaneFeatureCountLoader(
        count_files=count_files,
        meta_fns=metafiles,
        samples=samples,
        annotate_by=annotate_by,
        annotation_type=annotation_type
    )

    if units == 'counts':
        return obj.data
    elif units == 'fpkm':
        return obj.get_fpkm()
    elif units == 'tpm':
        return obj.get_tpm()

    return featurecounts(
        count_files,
        metafiles,
        samples=samples,
        units=units,
        annotate_by=annotate_by,
        annotation_type=annotation_type
    )


def gbm_paired_samples_loader(source='star', annotate_by='all', annotation_type='protein_coding'):
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        'GBM018',
        'GBM019',
        'GBM024',
        'GBM026',
        'GBM031',
        'DURA018N2_NSC',
        'DURA019N8C_NSC',
        'DURA024N28_NSC',
        'DURA026N31D_NSC',
        'DURA031N44B_NSC',
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



def gbm_astrocyte_nsc_samples(units='counts', annotate_by='all', annotation_type='protein_coding'):
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    count_files = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        'DURA018N2_NSC',
        'DURA019N8C_NSC',
        'DURA018N2_ASTRO_DAY12',
        'DURA019N8C_ASTRO_DAY12',
    )
    return featurecounts(
        count_files,
        metafiles,
        samples=samples,
        units=units,
        annotate_by=annotate_by,
        annotation_type=annotation_type
    )


def gbm_astrocyte_nsc_samples_loader(source='star', annotate_by='all', annotation_type='protein_coding'):
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        'DURA018N2_NSC',
        'DURA019N8C_NSC',
        'DURA018N2_ASTRO_DAY12',
        'DURA019N8C_ASTRO_DAY12',    )

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
        )
    else:
        raise ValueError("Unrecognised source.")
    return obj


def all_hgic_loader(source='star', annotate_by='all', annotation_type='protein_coding'):
    loaders = []

    # expt 1
    samples = (
        'GBM018',
        'GBM019',
        'GBM024',
        'GBM026',
        'GBM031',
        'DURA018N2_NSC',
        'DURA019N8C_NSC',
        'DURA024N28_NSC',
        'DURA026N31D_NSC',
        'DURA031N44B_NSC',
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
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type, samples=samples
        )
    )

    # expt 2
    samples = (
        "NH16_2383_P2",
        "GBM044_P4",
        "GBM026_P3_P4",
        "GBM018_P12",
        "GBM044_P8",
        "DURA044N8NSCP2",
        "GIBCONSC_P4",
        "DURA044N17_NSC_P3",
        "DURA018N4_NSC_P4",
        "GBM030_P5",
        "DURA030N16B6_NSC_P1",
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
            countdirs, metafiles, source=source, annotate_by=annotate_by, annotation_type=annotation_type, samples=samples
        )
    )

    return MultipleBatchLoader(loaders)



def mb_zhao_cultures(units='counts', annotate_by='all', annotation_type='protein_coding'):
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'wtchg_p160704')
    lane1dir = os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX')
    lane2dir = os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX')
    count_files = [os.path.join(d, 'featureCounts', 'counts.txt') for d in (lane1dir, lane2dir)]
    metafiles = [os.path.join(d, 'sources.csv') for d in (lane1dir, lane2dir)]
    samples = (
        '1078',
        '1595',
        '1487'
    )
    return featurecounts(
        count_files,
        metafiles,
        samples=samples,
        units=units,
        annotate_by=annotate_by,
        annotation_type=annotation_type
    )


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


def nih_gdc_gbm_preprocessed(units='counts'):
    """
    Load the preprocessed gene expression data downloaded from NIH genomic data commons
    :param units: Either 'counts' or 'fpkm' are available.
    :return: data, meta
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm')
    if units == 'counts':
        infile = os.path.join(indir, 'htseq_count', 'counts.csv')
    elif units == 'fpkm':
        infile = os.path.join(indir, 'htseq_count', 'fpkm.csv')
    else:
        raise ValueError("Unsupported units. Supported values are 'counts' and 'fpkm'.")

    meta_fn = os.path.join(indir, 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=0)

    if not os.path.exists(infile):
        logger.info("Unable to find summarised file %s. Calculating from individual files now.", infile)
        dat = []
        for cid in meta.index:
            if units == 'counts':
                this_infile = os.path.join(indir, 'htseq_count', 'counts', '%s.gz' % cid)
            else:
                this_infile = os.path.join(indir, 'htseq_count', 'fpkm', '%s.gz' % cid)
            try:
                this_dat = pd.read_csv(this_infile, header=None, index_col=0, sep='\t')
                # amend Ensembl IDs (remove version number)
                this_dat.index = this_dat.index.str.replace(r'\.[0-9]*$', '')
                dat.append(this_dat)
            except Exception:
                logger.exception("Failed to read file %s", this_infile)

        dat = pd.concat(dat, axis=1)
        dat.columns = meta.loc[:, 'sample']
        logger.info("Saving data from %d samples to %s.", dat.columns.size, infile)
        dat.to_csv(infile)
    else:
        dat = pd.read_csv(infile, header=0, index_col=0)

    meta.set_index('sample', inplace=True)
    return dat, meta


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


def gse61794(source='star', annotate_by='all', annotation_type='protein_coding'):
    """
    2 samples similar to Gibbco NSC line.
    These may be technical replicates as they are very highly correlated.
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
    indir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', '170328_K00150_0177_BHJ2C2BBXX')
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
