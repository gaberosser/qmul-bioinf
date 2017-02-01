import os
import pandas as pd
import numpy as np
from microarray import process
from settings import DATA_DIR, DATA_DIR_NON_GIT
import re
from gzip import GzipFile
from log import get_console_logger
logger = get_console_logger(__name__)


AGGREGATION_FIELD_CHOICES = (
    'ENTREZID', 'SYMBOL', 'ENSEMBL'
)



def load_from_r_processed(infile, sample_names, aggr_field=None, aggr_method=None):
    """
    Load microarray data that has been created in R.
    This is (unfortunately) necessary when we want to apply certain pre-processing steps like RMA.
    :param infile: The input file name.
    :param sample_names: List (or iterable) of the sample names to keep. If None, all will be kept.
    NB: if None and the input file contains columns other than the 'standard' annotations, they will NOT be dropped.
    :param aggr_field:
    :return:
    """
    if (aggr_method is not None and aggr_field is None) or (aggr_method is None and aggr_field is not None):
        raise ValueError("Must either supply BOTH aggr_field and aggr_method or NEITHER.")
    if aggr_field is not None and aggr_field not in AGGREGATION_FIELD_CHOICES:
        logger.warning("Unrecognised aggregation field. Supported options are %s.", ', '.join(AGGREGATION_FIELD_CHOICES))

    arr_data = pd.read_csv(infile, sep='\t', header=0, index_col=0)
    if aggr_field is None:
        return arr_data

    # drop other meta fields
    if sample_names is None:
        # we'll have to guess what the other annotation fields are then drop them - this might not work
        remaining_annot = set(AGGREGATION_FIELD_CHOICES)
        try:
            remaining_annot.remove(aggr_field)
        except KeyError:
            pass
        arr_data.drop(remaining_annot, axis=1, errors='ignore', inplace=True)
    else:
        # Use sample names to ensure we only keep data columns
        arr_data = arr_data.loc[:, [aggr_field] + sample_names]

    if aggr_method is not None:
        arr_data = process.aggregate_by_probe_set(arr_data, groupby=aggr_field, method=aggr_method)
    return arr_data


def load_gse28192_probeset():
    # load library and annotate
    indir = os.path.join(DATA_DIR, 'microarray_GSE28192')
    annot_fn = os.path.join(indir, 'probe_set', 'GPL6102-11574.txt')
    ann = pd.read_csv(annot_fn, sep='\t', comment='#', header=0, index_col=0)
    return ann


def annotate_and_aggregate_gse28192(data, aggr_field=None, aggr_method=None):
    ann_cols = {
        'ENTREZID': 'Entrez_Gene_ID',
        'SYMBOL': 'Symbol',
    }

    if aggr_field == 'all' and aggr_method is not None:
        raise ValueError("Cannot specify an aggregation method when aggr_field=='all'.")
    if aggr_field is None and aggr_method is not None:
        raise ValueError("Must specify an aggr_field if aggr_method is not None")
    if aggr_field is not None and aggr_field != 'all' and aggr_field not in ann_cols:
        raise ValueError("Unrecognised aggr_field %s." % aggr_field)

    probeset = load_gse28192_probeset()
    common_probes = probeset.index.intersection(data.index)
    if aggr_field is None:
        return data
    if aggr_field == 'all':
        # add all relevant fields
        for a, b in ann_cols.items():
            data.loc[common_probes, a] = probeset.loc[common_probes, b]
    else:
        # include only the aggregation field
        data.loc[common_probes, aggr_field] = probeset.loc[common_probes, ann_cols[aggr_field]]
        # aggregate
        data = process.aggregate_by_probe_set(data, groupby=aggr_field, method=aggr_method)

    return data


def load_annotated_gse28192(aggr_field=None, aggr_method=None, sample_names=None, max_pval=None, log2=True):
    """
    Xiao-Nan Li data comprising several tumour/xenograft matched MB samples
    :param aggr_field:
    :param aggr_method:
    :param max_pval: If specified, any expression values with less statistical support are replaced with 0.
    In practice, this just imposes a minimum value cutoff. However, it does ensure no negative values are returned.
    :return:
    """

    indir = os.path.join(DATA_DIR, 'microarray_GSE28192')
    meta_fn = os.path.join(indir, 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=1, sep=',')
    if sample_names is None:
        sample_names = list(meta.index)
    else:
        meta = meta.loc[sample_names, :]

    arr = {}
    for sn in sample_names:
        ff = os.path.join(indir, "%s.gz" % sn)
        df = pd.read_csv(
            ff,
            sep='\t',
            header=0,
            index_col=0,
            usecols=[0, 1, 4],
        )
        df.columns = ['expr', 'pval']
        if max_pval is not None:
            idx = df.loc[df.loc[:, 'pval'] > max_pval].index
            df.loc[idx, 'expr'] = 0.
        arr[sn] = df.loc[:, 'expr']
    arr = pd.DataFrame.from_dict(arr, dtype=float)

    # do log2 transform first if required
    if log2:
        arr[arr < 1.] = 1.
        arr = np.log2(arr)

    arr = annotate_and_aggregate_gse28192(arr, aggr_field=aggr_field, aggr_method=aggr_method)

    return arr, meta


def load_annotated_microarray_gse54650(index_field='entrez_id'):
    """
    Circadian study (24 timepoints) in mouse cerebellum tissue.
    :param index_field:
    :return:
    """
    infile = os.path.join(DATA_DIR, 'microarray_GSE54650', 'data.ann.txt.gz')

    sample_names = [
        'Cer_CT18',
        'Cer_CT20',
        'Cer_CT22',
        'Cer_CT24',
        'Cer_CT26',
        'Cer_CT28',
        'Cer_CT30',
        'Cer_CT32',
        'Cer_CT34',
        'Cer_CT36',
        'Cer_CT38',
        'Cer_CT40',
        'Cer_CT42',
        'Cer_CT44',
        'Cer_CT46',
        'Cer_CT48',
        'Cer_CT50',
        'Cer_CT52',
        'Cer_CT54',
        'Cer_CT56',
        'Cer_CT58',
        'Cer_CT60',
        'Cer_CT62',
        'Cer_CT64',
    ]
    return load_from_r_processed(infile, sample_names, aggr_field=index_field)


def load_annotated_microarray_sb_data(aggr_field='ENTREZ', aggr_method='max'):
    """
    Data from Dubuc for 8 MB samples in mouse. 3 have an inactivating CHD7 insertion.
    :param index_field:
    :return:
    """
    infile = os.path.join(DATA_DIR, 'sleeping_beauty_mouse_screen', 'data.ann.txt.gz')
    sample_names = [
        "Wu050",
        "Wu053",
        "Wu054",
        "Wu051",
        "Wu052",
        "Wu055",
        "Wu056",
        "Wu057"
    ]
    arr_data = load_from_r_processed(infile, sample_names, aggr_field=aggr_field, aggr_method=aggr_method)
    # CHD7 status
    chd7 = pd.Series(data=[True] * 3 + [False] * 5, index=sample_names)
    return arr_data, chd7


def _annotated_microarray(indir, aggr_field=None, aggr_method=None):
    """
    Base loader for microarray data.
    :param aggr_field: If None, don't perform aggregation. Otherwise, use this field as the index.
    :param aggr_method: The method to use when aggregating over probe sets. If None, do not aggregate over probe sets.
    :param indir:
    :param aggr_field:
    :param aggr_method:
    :return:
    """
    infile = os.path.join(indir, 'expr.rma.csv.gz')
    meta_fn = os.path.join(indir, 'sources.csv')
    if os.path.isfile(meta_fn):
        meta = pd.read_csv(meta_fn, header=0, index_col=0, sep=',')
        sample_names = list(meta.index)
    else:
        logger.warning("Unable to find meta file %s", meta_fn)
        meta = None
        sample_names = None  # this will cause an error if aggr_field is not None

    arr = load_from_r_processed(infile, sample_names, aggr_field=aggr_field, aggr_method=aggr_method)
    return arr, meta


def load_annotated_microarray_gse37382(aggr_field=None, aggr_method=None):
    """
    Northcott human MB study comprising the subgroups C, D and SHH.
    :param aggr_field: If None, don't perform aggregation. Otherwise, use this field as the index. Options ae
    :param aggr_method: The method to use when aggregating over probe sets. If None, do not aggregate over probe sets.
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'microarray', 'GSE37382')
    infile = os.path.join(indir, 'expr.rma.csv.gz')
    meta_fn = os.path.join(indir, 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=0, sep=',')
    meta.index = meta.index.str.replace('_', '.')  # R changes column names
    sample_names = list(meta.index)
    arr = load_from_r_processed(infile, sample_names, aggr_field=aggr_field, aggr_method=aggr_method)
    return arr, meta


def load_annotated_microarray_gse10327(aggr_field=None, aggr_method=None):
    """
    Kool human MB dataset comprising 62 samples
    :param index_field:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'microarray', 'GSE10327')
    infile = os.path.join(indir, 'expr.rma.csv.gz')
    meta_fn = os.path.join(indir, 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=0, sep=',')
    sample_names = list(meta.index)
    arr = load_from_r_processed(infile, sample_names, aggr_field=aggr_field, aggr_method=aggr_method)
    return arr, meta


def load_annotated_microarray_gse37418(aggr_field=None, aggr_method=None):
    """
    Robinson human MB dataset comprising ~70 samples
    :param index_field:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'microarray', 'GSE37418')
    infile = os.path.join(indir, 'expr.rma.csv.gz')
    meta_fn = os.path.join(indir, 'sources.csv')
    meta = pd.read_csv(meta_fn, header=0, index_col=1, sep=',')  # NB index by accession here
    sample_names = list(meta.index)
    arr = load_from_r_processed(infile, sample_names, aggr_field=aggr_field, aggr_method=aggr_method)
    return arr, meta


def load_annotated_microarray_gse33199(aggr_field=None, aggr_method=None):
    """
    Roussel mouse MB dataset.  Filtered to include only the essential samples, as this dataset contains
    a large number of modified samples.
    :param aggr_field:
    :param aggr_method:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'microarray', 'GSE33199')
    return _annotated_microarray(indir, aggr_field=aggr_field, aggr_method=aggr_method)
    # infile = os.path.join(indir, 'expr.rma.csv.gz')
    # meta_fn = os.path.join(indir, 'sources.csv')
    # meta = pd.read_csv(meta_fn, header=0, index_col=1, sep=',')  # NB index by accession here
    # sample_names = list(meta.index)
    # arr = load_from_r_processed(infile, sample_names, aggr_field=aggr_field, aggr_method=aggr_method)
    # return arr, meta


def load_annotated_thompson2006(aggr_field=None, aggr_method=None):
    """
    Thompson 46 human MB samples. No metadata available.
    :param aggr_field:
    :param aggr_method:
    :return:
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'microarray', 'thompson2006')
    return _annotated_microarray(indir, aggr_field=aggr_field, aggr_method=aggr_method)


def load_metadata_from_series_matrix(infile):
    """
    Load metadata from the specified text matrix (available for many GEO datasets)
    :param infile: .txt or .txt.gz series matrix
    :return:
    """
    meta_map = {
        'Sample_title': 'title',
        'Sample_geo_accession': 'accession',
        'Sample_source_name_ch1': 'description',
    }
    nested_headers = (
        'Sample_characteristics_ch1',
    )

    if re.search(re.compile(r'\.gz$', re.IGNORECASE), infile):
        opener = lambda x: GzipFile(filename=x, mode='rb')
    else:
        opener = lambda x: open(x, 'rb')

    meta = pd.DataFrame()

    with opener(infile) as f:
        while True:
            line = f.readline().strip('\n')
            if len(line) == 0:
                continue
            header = re.match(r'^!(?P<hd>[^\t]*)', line).group('hd')
            if header in nested_headers:
                the_data = [t.strip('"') for t in line.split('\t')[1:]]
                hdr = re.match(r'^(?P<hd>[^:]*)', the_data[0]).group('hd')
                the_data = [re.sub(r'%s: ' % hdr, '', t) for t in the_data]
                meta[hdr] = the_data
            elif header in meta_map:
                the_data = [t.strip('"') for t in line.split('\t')[1:]]
                meta[meta_map[header]] = the_data
            if line == '!series_matrix_table_begin':
                break

    return meta