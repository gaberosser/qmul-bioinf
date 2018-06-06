import pandas as pd
import os
import copy
from utils import log
from utils import setops, string_manipulation
import fnmatch
import numpy as np


def recursive_file_search(root_dir, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches



def _filter_by_sample_name_index(sample_names, filt, exact=True):
    # def filter_by_sample_name(self, filt, exact=True):
    sn = pd.Index(sample_names)

    if hasattr(filt, "__iter__"):
        if not exact:
            filt = string_manipulation.construct_multiple_or_regex(filt)
    elif exact:
        filt = [filt]

    if not exact:
        idx = sn.str.contains(filt)
    else:
        idx = sn.isin(filt)

    if idx.sum() == 0:
        raise AttributeError("The requested filtering operation would remove all samples")

    return idx


class DatasetLoader(object):
    meta_col_sample_name = 'sample'
    meta_col_filename = 'filename'
    extra_df_attributes = tuple()

    def __init__(
            self,
            base_dir=None,
            meta_fn=None,
            samples=None,
            tax_id=9606,
            batch_id=None,
            verbose=True,
            *args,
            **kwargs):

        """
        Base class for loading a dataset.
        :param base_dir: Path to the root input directory. All data must be contained in this directory
        or below it.
        :param meta_fn: Path to the meta file.
        :param samples: If supplied, use this to filter the files loaded.
        :param tax_id: The taxonomy ID (default: 9606, human)
        :param batch_id: Optionally supply a name for this batch, useful when combining batches
        """
        self.base_dir = base_dir
        if not os.path.isdir(self.base_dir):
            raise ValueError("Supplied base_dir %s does not exist or is not a directory." % self.base_dir)

        self.meta_fn = meta_fn
        if self.meta_fn is not None:
            if not os.path.isfile(self.meta_fn):
                raise ValueError("Meta file %s does not exist." % self.meta_fn)

        self.meta_is_linked = None
        self.sample_names = None

        self.samples_to_keep = samples
        self.tax_id = tax_id
        self.batch_id = batch_id
        self.verbose = verbose

        self.logger = log.get_console_logger(self.__class__.__name__)

        self.meta = None
        self.input_files = None
        self.data = None

        self.load_meta()

        self.get_inputs()
        self.load_data()
        self.post_process()

    def load_meta(self):
        """
        Load the metadata, if a file has been specified
        :return:
        """
        if self.meta_fn is not None:
            self.meta = pd.read_csv(self.meta_fn, header=0, index_col=None, dtype=object)
            if self.meta_col_sample_name not in self.meta.columns:
                self.logger.warning(
                    "The meta data has no column '%s', so sample names will be inferred from the data.",
                    self.meta_col_sample_name
                )
                self.meta_is_linked = False
            else:
                self.meta.set_index(self.meta_col_sample_name, inplace=True)
                if self.meta_col_filename == self.meta_col_sample_name:
                    # the filename and sample name columns are shared, so we need to replace it
                    self.meta.insert(0, self.meta_col_filename, self.meta.index)
                self.meta_is_linked = True


    def get_inputs(self):
        """
        Set self.input_files based on the metadata (if loaded) and base_dir
        """
        raise NotImplementedError

    def load_one_file(self, fn):
        """
        Load a single datum or multiple data from a single file
        :param fn: Path to the file
        """
        raise NotImplementedError

    def load_data(self):
        """
        Load raw data
        :return:
        """
        raise NotImplementedError

    def post_process(self):
        """
        Carry out any processing required after loading, modifying meta and data accordingly
        """

    def filter_by_sample_name(self, filt, exact=True):
        if not self.meta_is_linked:
            raise NotImplementedError("Meta data are not linked to the data, so we can't filter. Run it manually on the data attribute.")
        idx = _filter_by_sample_name_index(self.meta.index, filt, exact=exact)
        self.meta = self.meta.loc[idx]
        self.data = self.data.loc[:, self.meta.index]

    def filter_samples(self, keep_idx):
        """
        Drop samples according to the supplied index. All data structures are updated in-place.
        :param drop_idx: Either a boolean index or an iterable that can be used as an indexer. The samples indicated
        will be retained.
        :return:
        """
        self.meta = self.meta.loc[keep_idx]
        self.data = self.data.loc[:, self.meta.index]


class SingleFileLoader(DatasetLoader):
    def load_data(self):
        self.data = self.load_one_file(self.input_files)

    def post_process(self):
        if self.meta_is_linked:
            # check for discrepancies and report
            s1 = self.meta.index.difference(self.data.columns)
            s2 = self.data.columns.difference(self.meta.index)
            if len(s1) > 0:
                self.logger.warn(
                    "%d samples were included in the meta data but not found in the data: %s",
                    len(s1),
                    ', '.join(s1)
                )
            if len(s2) > 0:
                self.logger.warn(
                    "%d samples were included in the data but not listed in the meta data: %s",
                    len(s2),
                    ', '.join(s2)
                )

        # apply sample filtering if requested
        if self.samples_to_keep is not None:
            self.data = self.data.loc[:, self.samples_to_keep]

        # ensure that meta has the same entries as data
        if self.meta_is_linked:
            self.meta = self.meta.loc[self.data.columns]
        super(SingleFileLoader, self).post_process()


# TODO idea: create a mixin that replaces load_data in the MultipleFileLoader with a more appropriate method when the
# data are NOT row-matched. For example. ChIP peaks are not joined on the index, so it's fairly pointless trying to
# stuff them into a single DataFrame?

def join_row_indexed_data(**kwargs):
    dat = None
    for sn, this_dat in kwargs.iteritems():
        # if each sample if represented by a simple series, the joining process is straightforward:
        if isinstance(this_dat, pd.Series):
            if dat is None:
                dat = pd.DataFrame(index=this_dat.index)
            dat.insert(dat.shape[1], sn, this_dat)

        # otherwise, we need to use a MultiIndex
        elif isinstance(this_dat, pd.DataFrame):
            # make a multiindex column
            this_dat = pd.DataFrame(
                this_dat.values,
                index=this_dat.index,
                columns=pd.MultiIndex.from_product([[sn], this_dat.columns], names=['sample', 'fields'])
            )
            if dat is None:
                dat = this_dat
            else:
                dat = dat.join(this_dat, how='outer')  # outer join to keep all values
    return dat


class MultipleFileLoader(DatasetLoader):
    file_pattern = ''
    row_indexed = True

    def generate_input_path(self, fname):
        """
        Given the filename from the meta file, generate the path to the actual data (e.g. constructing subdir structure)
        """
        return os.path.join(self.base_dir, fname)

    def generate_sample_name(self, file_path):
        """
        Given the input file path, generate the sample name. Used when there is no metatdata
        """
        # default behaviour: do nothing
        return file_path

    def load_meta(self):
        super(MultipleFileLoader, self).load_meta()
        if self.meta is not None:
            if self.meta_col_filename not in self.meta.columns:
                self.logger.warning(
                    "The meta data has no column '%s', so file names will not be used.",
                    self.meta_col_filename
                )
                self.meta_is_linked = False
            else:
                self.meta_is_linked = True

    def get_inputs(self):
        inputs = pd.Series()
        if self.meta_is_linked:
            # get the filenames directly from the meta data
            not_found = []
            for sn, f in self.meta.loc[:, self.meta_col_filename].iteritems():
                ff = self.generate_input_path(f)
                if os.path.isfile(ff):
                    inputs.loc[sn] = ff
                else:
                    not_found.append(sn)

            if inputs.size == 0:
                raise AttributeError("No files listed in the meta data were found.")

            if len(not_found):
                if self.verbose:
                    self.logger.warn(
                        "Unable to find %d datasets listed in the metadata: %s",
                        len(not_found),
                        ', '.join(not_found)
                    )
                    # remove from the samples_to_keep if supplied
                    if self.samples_to_keep is not None:
                        ## FIXME?
                        requested_but_missing = sorted(set(self.samples_to_keep).intersection(not_found))
                        if len(requested_but_missing) > 0:
                            self.logger.warn(
                                "Of the %d samples requested, %d cannot be loaded as they are missing: %s.",
                                len(self.samples_to_keep),
                                len(requested_but_missing),
                                ', '.join(requested_but_missing)
                            )
                            for t in self.samples_to_keep:
                                if t in requested_but_missing:
                                    self.samples_to_keep.remove(t)

        else:
            # no meta file names -> find the files by pattern
            flist = recursive_file_search(self.base_dir, self.file_pattern)
            for f in flist:
                if not os.path.isdir(f):
                    if os.path.getsize(f) > 0:
                        sn = self.generate_sample_name(f)
                        inputs.loc[sn] = f
                    else:
                        self.logger.warning(
                            "File %s should be included but has size 0. Did something go wrong at creation?",
                            f
                        )
            if len(inputs) == 0:
                raise AttributeError("No files matching the file pattern '%s' were found in %s" % (
                    self.file_pattern,
                    self.base_dir
                ))

        if self.samples_to_keep is not None:
            inputs = inputs.loc[self.samples_to_keep]

        self.input_files = inputs
        self.logger.info("Loading %d samples: %s", self.input_files.shape[0], ', '.join(self.input_files.index))

    def load_data(self):
        the_data_separates = {}
        for sn, fn in self.input_files.iteritems():
            the_data_separates[sn] = self.load_one_file(fn)
        if self.row_indexed:
            dat = join_row_indexed_data(**the_data_separates)
            # tolist() is important here for MultiIndexed data
            self.data = dat[self.input_files.index.tolist()]
        else:
            self.data = the_data_separates

    def post_process(self):
        """
        Carry out any processing required after loading, modifying meta and data accordingly
        """
        if self.meta_is_linked:
            # ensure that meta has the same entries as data
            self.meta = self.meta.loc[self.input_files.index]


class MultipleBatchLoader(object):
    ## FIXME: respect the batch ID somehow - so that MultipleBatch instances can also be combined
    def __init__(self, loaders, intersection_only=True):
        """
        Class to combine multiple loader objects.
        Each loader represents a separate batch. Inputs can include multiple lane loaders.
        :param loaders: Iterable of loader objects.
        :param intersection_only: If True (default), reduce counts to the indices (e.g. genes) that are present in all
        loaders.
        """
        self.logger = log.get_console_logger(self.__class__.__name__)

        if len(loaders) < 2:
            raise ValueError("Must supply 2 or more loaders to use a MultipleBatchLoader.")

        # we can only claim the meta data is linked here if all loaders have this property
        self.meta_is_linked = True
        for l in loaders:
            if not l.meta_is_linked:
                self.meta_is_linked = False

        # set the batch  column name avoiding clashes
        batch_col = 'batch'
        meta_cols = sorted(setops.reduce_union(*[t.meta.columns for t in loaders if t.meta is not None]))

        if batch_col in meta_cols:
            i = 1
            while batch_col in meta_cols:
                batch_col = "batch_%d" % i
                i += 1
        meta_cols += [batch_col]

        # check attributes that must match in all loaders
        if len(set([t.tax_id for t in loaders])) > 1:
            raise AttributeError(
                "The tax_id of the samples differ between loaders: %s" % ', '.join([str(t.tax_id) for t in loaders])
            )
        else:
            self.tax_id = loaders[0].tax_id

        if len(set([t.row_indexed for t in loaders])) > 1:
            raise AttributeError("row_indexed bool must be the same in all loaders")
        else:
            self.row_indexed = loaders[0].row_indexed

        extra_df_attributes = {}

        if self.row_indexed:
            row_indexed_dat_arr = {}
        else:
            dat = {}

        meta_values = []
        meta_index = []
        blank_meta_row = dict([(k, None) for k in meta_cols])

        # we may need to append a number to sample names
        sample_appendix = 0
        auto_batch = 1
        meta_auto_idx = 0
        samples_seen = set()

        for l in loaders:
            this_batch = l.batch_id
            if not hasattr(this_batch, '__iter__'):
                if l.batch_id is None:
                    this_batch = auto_batch
                    auto_batch += 1
                this_batch = pd.Series(this_batch, index=l.meta.index)

            try:
                this_samples = l.input_files.index.tolist()
            except AttributeError:
                # occurs when we are loading a single file
                # FIXME: find a better catch - this is too general
                if hasattr(l, 'input_files'):
                    # this occurs if l is a single file loader
                    this_samples = [l.input_files]
                else:
                    # this occurs if l is a batch loader
                    # FIXME: may not give us valid sample names?
                    this_samples = l.meta.index

            # get a copy of the data
            if self.row_indexed:
                this_dat = l.data.copy()
            else:
                this_dat = copy.copy(l.data)

            # get a copy of meta
            if l.meta is not None:
                this_meta = l.meta.copy()

            # resolve any sample clashes in the data (NOT the meta data)
            clash_resolved = False
            new_names = []

            while len(samples_seen.intersection(this_samples)) > 0:
                sample_appendix += 1
                # find the clash
                clashes = samples_seen.intersection(this_samples)
                self.logger.warning(
                    "Found sample name clash(es): %s. Modifying names to avoid errors.",
                    ', '.join(clashes)
                )
                for c in clashes:
                    new_names.append([
                        this_samples[this_samples.index(c)],
                        this_samples[this_samples.index(c)] + "_%d" % sample_appendix
                    ])
                    this_samples[this_samples.index(c)] += "_%d" % sample_appendix
                clash_resolved = True
            samples_seen.update(this_samples)

            if clash_resolved:
                # relabel metadata if linked
                if l.meta_is_linked:
                    # reorder first to be sure it's the same as data
                    this_meta = this_meta.loc[this_dat.columns]
                    this_meta.index = this_samples

                # relabel the data
                if self.row_indexed:
                    this_dat.columns = this_samples
                else:
                    for prev, new in new_names:
                        this_dat[new] = this_dat.pop(prev)

                # relabel the batch IDs
                this_batch.index = this_samples
                # relabel any other DF data if present
                for fld in l.extra_df_attributes:
                    x = getattr(l, fld)
                    x.columns = this_samples

            # data
            if self.row_indexed:
                if isinstance(this_dat.columns, pd.MultiIndex):
                    col_list = this_dat.columns.levels[0].tolist()
                else:
                    col_list = this_dat.columns.tolist()
                for c in col_list:
                    row_indexed_dat_arr[c] = this_dat[[c]]

            else:
                dat.update(this_dat)

            # other df attributes
            for fld in l.extra_df_attributes:
                if fld not in extra_df_attributes:
                    extra_df_attributes[fld] = getattr(l, fld).copy()
                else:
                    extra_df_attributes[fld] = pd.concat((extra_df_attributes[fld], getattr(l, fld)), axis=1)

            # rebuild meta
            if l.meta is not None:
                for i in this_meta.index:
                    this_row = dict(blank_meta_row)
                    this_row.update(this_meta.loc[i].to_dict())
                    this_row[batch_col] = this_batch[i]
                    meta_values.append(this_row)
                    if l.meta_is_linked:
                        meta_index.append(i)
                    else:
                        meta_index.append(meta_auto_idx)
                        meta_auto_idx += 1
            else:
                for c in this_dat.columns:
                    this_row = dict(blank_meta_row)
                    this_row[batch_col] = this_batch[c]
                    meta_values.append(this_row)
                    meta_index.append(meta_auto_idx)
                    meta_auto_idx += 1

        self.meta = pd.DataFrame(meta_values, index=meta_index, columns=meta_cols)

        if self.row_indexed:
            dat = pd.concat(
                [row_indexed_dat_arr[k] for k in self.meta.index],
                axis=1
            )

        self.data = dat
        self.batch_id = self.meta.loc[:, batch_col]

        self.extra_df_attributes = tuple()
        for fld in extra_df_attributes:
            setattr(self, fld, extra_df_attributes[fld])
            self.extra_df_attributes += (fld,)

    def filter_by_sample_name(self, filt, exact=True):
        if not self.meta_is_linked:
            raise NotImplementedError("Meta data are not linked to the data, so we can't filter. Run it manually on the data attribute.")
        idx = _filter_by_sample_name_index(self.meta.index, filt, exact=exact)
        self.meta = self.meta.loc[idx]
        self.data = self.data.loc[:, self.meta.index]
        self.batch_id = self.batch_id.loc[idx]

    def filter_samples(self, keep_idx):
        """
        Drop samples according to the supplied index. All data structures are updated in-place.
        :param drop_idx: Either a boolean index or an iterable that can be used as an indexer. The samples indicated
        will be retained.
        :return:
        """
        self.meta = self.meta.loc[keep_idx]
        self.data = self.data.loc[:, self.meta.index]
        self.batch_id = self.batch_id.loc[self.meta.index]
