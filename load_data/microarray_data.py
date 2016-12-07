import os
import pandas as pd
from microarray import process
from settings import DATA_DIR, DATA_DIR_NON_GIT


AGGREGATION_FIELD_CHOICES = (
    'ENTREZID', 'SYMBOL', 'ENSEMBL'
)



def load_from_r_processed(infile, sample_names, aggr_field=None, aggr_method=None):
    """
    Load microarray data that has been created in R.
    This is (unfortunately) necessary when we want to apply certain pre-processing steps like RMA.
    :param infile: The input file name.
    :param sample_names: T
    :param aggr_field:
    :return:
    """
    if (aggr_method is not None and aggr_field is None) or (aggr_method is None and aggr_field is not None):
        raise ValueError("Must either supploy BOTH aggr_field and aggr_method or NEITHER.")
    if aggr_field is not None and aggr_field not in AGGREGATION_FIELD_CHOICES:
        raise ValueError("Unrecognised aggregation field. Supported options are %s." % ', '.join(AGGREGATION_FIELD_CHOICES))

    arr_data = pd.read_csv(infile, sep='\t', header=0, index_col=0)
    if aggr_field is None:
        return arr_data
    # drop other meta fields
    arr_data = arr_data.loc[:, [aggr_field] + sample_names]
    arr_data = process.aggregate_by_probe_set(arr_data, groupby=aggr_field, method=aggr_method)
    return arr_data


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


def load_annotated_microarray_sb_data(index_field='entrez_id'):
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
    arr_data = load_from_r_processed(infile, sample_names, aggr_field=index_field)
    # CHD7 status
    chd7 = pd.Series(data=[True] * 3 + [False] * 5, index=sample_names)
    return arr_data, chd7


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
    Kool human MB dataset comprising 62 samples
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

