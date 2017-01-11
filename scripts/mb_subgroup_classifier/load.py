from load_data import rnaseq_data, microarray_data
import pandas as pd
from microarray import process
from scripts.comparison_rnaseq_microarray import load_rnaseq_data


def load_xz_rnaseq(kind='cuff', yugene=True, gene_symbols=None):
    """
    Load RNA-Seq data from XZ samples
    :param kind: The source of the data, either 'cuff' or 'htseq'
    :param yugene: If True, apply YuGene normalisation
    :param gene_symbols: If supplied, this is a list containing the gene symbols. Any that are not present are filled
    with zeros
    :return:
    """
    if kind == 'cuff':
        X = load_rnaseq_data.load_rnaseq_cufflinks_gene_count_data(unit='fpkm')
        if yugene:
            X = process.yugene_transform(X)
        if gene_symbols is not None:
            X = pd.DataFrame(data=X, columns=X.columns, index=gene_symbols)
            X.fillna(0, inplace=True)
    elif kind == 'htseq':
        X = rnaseq_data.gse83696(index_by='Approved Symbol')
        if yugene:
            X = process.yugene_transform(X)
        if gene_symbols is not None:
            X = pd.DataFrame(data=X, columns=X.columns, index=gene_symbols)
            X.fillna(0, inplace=True)
    else:
        raise ValueError("Unrecognised kind '%s'" % kind)
    return X


def load_xiaonan_microarray(yugene=True, gene_symbols=None, sample_names=None):
    """
    Load the Xiao-Nan microarray data
    :param yugene: If True, apply YuGene normalisation
    :param gene_symbols: If supplied, this is a list containing the gene symbols. Any that are not present are filled
    with zeros
    :return:
    """
    X, meta = microarray_data.load_annotated_gse28192(
        aggr_field='SYMBOL',
        aggr_method='max',
        log2=True,
        sample_names=sample_names
    )
    if yugene:
        X = process.yugene_transform(X)
    if gene_symbols is not None:
        X = pd.DataFrame(data=X, columns=X.columns, index=gene_symbols)
        X.fillna(0, inplace=True)
    return X, meta


def load_sb_rnaseq(yugene=True, gene_symbols=None):
    """
    Load RNA-Seq from SB samples, counted using featureCounts
    :param yugene: If True, apply YuGene normalisation
    :param gene_symbols: If supplied, this is a list containing the gene symbols. Any that are not present are filled
    with zeros
    :return:
    """
    X = rnaseq_data.mb_zhao_cultures(units='fpkm', annotate_by='Approved Symbol')
    if yugene:
        X = process.yugene_transform(X)
    if gene_symbols is not None:
        X = pd.DataFrame(data=X, columns=X.columns, index=gene_symbols)
        X.fillna(0, inplace=True)
    return X