import pandas as pd
import os
from utils import log
from utils import setops
import fnmatch
import numpy as np


def recursive_file_search(root_dir, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches


class DatasetLoader(object):
    meta_col_sample_name = 'sample'
    meta_col_filename = 'filename'

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


class MultipleFileLoader(DatasetLoader):
    file_pattern = ''

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
            for i, f in self.meta.loc[:, self.meta_col_filename].iteritems():
                ff = self.generate_input_path(f)
                if os.path.isfile(ff):
                    sn = i
                    inputs.loc[sn] = ff
                else:
                    not_found.append(f)
            if len(not_found):
                if self.verbose:
                    self.logger.warn(
                        "Unable to find %d datasets listed in the metadata: %s",
                        len(not_found),
                        ', '.join(not_found)
                    )
            if inputs.size == 0:
                raise AttributeError("No files listed in the meta data were found.")

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
        dat = None
        for sn, fn in self.input_files.iteritems():
            this_dat = self.load_one_file(fn)
            if dat is None:
                dat = pd.DataFrame(index=this_dat.index)
            dat.insert(dat.shape[1], sn, this_dat)
        self.data = dat

    def post_process(self):
        """
        Carry out any processing required after loading, modifying meta and data accordingly
        """
        if self.meta_is_linked:
            # ensure that meta has the same entries as data
            self.meta = self.meta.loc[self.data.columns]


class MultipleBatchLoader(object):
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

        # get the data index (e.g. gene labels) based on intersection_only
        if intersection_only:
            idx = setops.reduce_intersection(*[t.data.index for t in loaders])
        else:
            idx = setops.reduce_union(*[t.data.index for t in loaders])

        # set the batch  column name avoiding clashes
        batch_col = 'batch'
        meta_cols = sorted(setops.reduce_union(*[t.meta.columns for t in loaders if t.meta is not None]))

        if batch_col in meta_cols:
            i = 1
            while batch_col in meta_cols:
                batch_col = "batch_%d" % i
                i += 1
        meta_cols += [batch_col]

        if len(set([t.tax_id for t in loaders])) > 1:
            raise AttributeError(
                "The tax_id of the samples differ between loaders: %s" % ', '.join([str(t.tax_id) for t in loaders])
            )

        dat = pd.DataFrame(index=idx, dtype=float)
        meta_values = []
        meta_index = []
        blank_meta_row = dict([(k, None) for k in meta_cols])
        meta = pd.DataFrame(columns=meta_cols)

        # we may need to append a number to sample names
        sample_appendix = 0
        auto_batch = 1
        meta_auto_idx = 0
        samples_seen = set()

        for l in loaders:
            this_batch = l.batch_id
            if l.batch_id is None:
                this_batch = auto_batch
                auto_batch += 1

            this_samples = l.data.columns.tolist()
            # get a copy of the data
            this_dat = l.data.copy()
            # get a copy of meta
            if l.meta is not None:
                this_meta = l.meta.copy()


            # resolve any sample clashes in the data (NOT the meta data)
            clash_resolved = False
            while len(samples_seen.intersection(this_samples)) > 0:
                sample_appendix += 1
                # find the clash
                clashes = samples_seen.intersection(this_samples)
                self.logger.warning(
                    "Found sample name clash(es): %s. Modifying names to avoid errors.",
                    ', '.join(clashes)
                )
                for c in clashes:
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
                this_dat.columns = this_samples

            # data
            for c in this_dat.columns:
                dat.insert(dat.shape[1], c, this_dat.loc[idx, c])

            # rebuild meta

            if l.meta is not None:
                for i in this_meta.index:
                    this_row = dict(blank_meta_row)
                    this_row.update(this_meta.loc[i].to_dict())
                    this_row[batch_col] = this_batch
                    meta_values.append(this_row)
                    if l.meta_is_linked:
                        meta_index.append(i)
                    else:
                        meta_index.append(meta_auto_idx)
                        meta_auto_idx += 1
            else:
                for c in this_dat.columns:
                    this_row = dict(blank_meta_row)
                    this_row[batch_col] = this_batch
                    meta_values.append(this_row)
                    meta_index.append(meta_auto_idx)
                    meta_auto_idx += 1
                    # meta.loc[c, batch_col] = this_batch

        self.meta = pd.DataFrame(meta_values, index=meta_index, columns=meta_cols)
        self.data = dat
