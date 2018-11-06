import os
import csv
import pandas as pd
from settings import GIT_LFS_DATA_DIR


def known_genes(tax_id=9606, index_field=None):
    """
    Get dataframe of known genes for the supplied taxonomy ID (default: human)
    Other options: mouse (10090)
    :param index_field: If supplied, reindex the dataframe by this field. Examples are 'Symbol', 'GeneID'
    """
    if tax_id == 9606:
        fn = os.path.join(GIT_LFS_DATA_DIR, 'ncbi_gene', 'Homo_sapiens.gene_info.gz')
    elif tax_id == 10090:
        fn = os.path.join(GIT_LFS_DATA_DIR, 'ncbi_gene', 'Mus_musculus.gene_info.gz')
    else:
        raise ValueError("Unrecognised tax_id %d.", tax_id)
    df = pd.read_csv(fn, sep='\t')
    if index_field is not None:
        df = df.set_index(index_field)
    return df


def conversion_table(type='protein_coding', tax_id=9606):
    ## TODO: add a version option
    if tax_id == 9606:
        if type == 'protein_coding':
            in_file = os.path.join(GIT_LFS_DATA_DIR, 'genenames', 'protein_coding', 'genenames.org.2017.01.tsv')
        elif type == 'all':
            # in_file = os.path.join(GIT_LFS_DATA_DIR, 'genenames', 'all', 'genenames.org.2017.09.tsv')
            # in_file = os.path.join(GIT_LFS_DATA_DIR, 'genenames', 'all', 'genenames.org.2018.03.tsv')
            in_file = os.path.join(GIT_LFS_DATA_DIR, 'genenames', 'all', 'genenames.org.2018.11.tsv')
        else:
            raise ValueError("Unsupported type option '%s'" % type)
        # in_file = os.path.join(DATA_DIR, 'genenames', 'genenames.org.tsv')
        df = pd.read_csv(in_file, delimiter='\t')
    elif tax_id == 10090:
        in_file = os.path.join(GIT_LFS_DATA_DIR, 'biomart', 'mm10', 'mm10.csv')
        df = pd.read_csv(in_file, header=0, index_col=None)
    else:
        raise NotImplementedError("Unknown taxonomy ID")
    return df


def _translate(cat, x, to_field, from_field):
    return cat.set_index(from_field).loc[x, to_field]


def translate(x, to_field, from_field, tax_id=9606, type='all'):
    cat = conversion_table(type=type, tax_id=tax_id)
    return cat.set_index(from_field).loc[x, to_field]


def gene_symbol_to_entrez(g, tax_id=9606):
    cat = conversion_table(type='all', tax_id=tax_id)
    return _translate(cat, g, 'Entrez Gene ID', 'Approved Symbol')


def entrez_to_gene_symbol(e, tax_id=9606):
    cat = conversion_table(type='all', tax_id=tax_id)
    return _translate(cat, e, 'Approved Symbol', 'Entrez Gene ID')


def entrez_to_ensembl(e, tax_id=9606):
    cat = conversion_table(type='all', tax_id=tax_id)
    return _translate(cat, e, 'Ensembl Gene ID', 'Entrez Gene ID')


def gene_symbol_to_ensembl(g, tax_id=9606):
    cat = conversion_table(type='all', tax_id=tax_id)
    return _translate(cat, g, 'Ensembl Gene ID', 'Approved Symbol')


def ensembl_to_gene_symbol(e, tax_id=9606):
    cat = conversion_table(type='all', tax_id=tax_id)
    return _translate(cat, e, 'Approved Symbol', 'Ensembl Gene ID')


def ensembl_to_name(e, tax_id=9606):
    cat = conversion_table(type='all', tax_id=tax_id)
    return _translate(cat, e, 'Approved Name', 'Ensembl Gene ID')


def translate_quantification_resolving_duplicates(dat, from_field, to_field, type='all', tax_id=9606, resolution='mean'):
    """
    Translate the index used for the input quantification data to some other
    :param dat: Pandas DataFrame containing some kind of quantification data
    :param to_field:
    :param from_field:
    :param type:
    :param tax_id:
    :param resolution: Approach used to resolve duplicates (TODO: add more as required).
    :return:
    """
    gs = translate(dat.index, to_field, from_field, tax_id=tax_id)
    gs = gs.loc[~gs.index.duplicated()]

    gs.dropna(inplace=True)
    dat = dat.loc[gs.index]
    dat.index = gs

    # this may leave some duplicate values to resolve
    dupe_idx = dat.index[dat.index.duplicated()]
    dupe_map = dat.index.isin(dupe_idx)
    if dupe_map.sum() != 0:
        dupes = dat.loc[dupe_map]
        dat = dat.loc[~dupe_map]
        if resolution == 'mean':
            dupes_aggr = dupes.groupby(dupes.index).mean()
        else:
            raise NotImplementedError("Unsupported resolution method %s." % resolution)
        dat = dat.append(dupes_aggr)
    return dat


def gene_symbol_to_x_robust(
        g,
        to_field,
        tax_id=9606,
        additional_search_fields=('Previous Symbols', 'Synonyms'),
        split_str=', ',
        type='all'
):
    """
    Translate from gene symbol to another compatible field.
    Where a search fails, we also look in the fields in additional_search_fields. If a match is found, this is used
    instead. Priority is given to fields mentioned earlier.
    :param g: Gene symbol or list of gene symbols
    :param to_field: Lookup column in the conversion table
    :param tax_id:
    :param additional_search_fields:
    :param split_str: String used to split the field (in the event that it contains multiple entries). If None, no
    splitting is attempted. This is necessary unless the entries are lists already.
    :param type:
    :return:
    """
    if not hasattr(g, '__iter__'):
        g = [g]
    cat = conversion_table(type=type, tax_id=tax_id).set_index('Approved Symbol')
    res = cat.loc[g, to_field]
    spares = pd.Index(g).difference(cat.index).tolist()
    for asf in additional_search_fields:
        if not len(spares):
            return res
        lookup = cat.loc[:, asf]
        if split_str is not None:
            lookup = lookup.fillna('').str.split(split_str)
        new_spares = []
        for t in spares:
            ix = lookup.apply(lambda x: t in x)
            if ix.sum() == 1:
                res.loc[t] = cat.loc[ix, to_field].values
            else:
                new_spares.append(t)
        spares = new_spares
    return res