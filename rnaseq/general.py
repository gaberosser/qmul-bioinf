import references
import pandas as pd
import os
from settings import LOCAL_DATA_DIR


def top_genes(
        data,
        n=100,
        convert_to_symbols=True,
        tax_id=9606,
):
    """
    Retrieve the top n genes from the data
    :param data: Indexed by ensembl_ID
    :param units:
    :param n:
    :return:
    """
    if convert_to_symbols:
        # get gene symbols and drop all NaN
        gs = references.ensembl_to_gene_symbol(data.index, tax_id=tax_id).dropna()
        gs = gs.loc[~gs.index.duplicated()]
        gs = gs.loc[~gs.duplicated()]
    res = {}
    for col in data.columns:
        t = data.loc[:, col].sort_values(ascending=False)[:n]
        if convert_to_symbols:
            new_idx = gs.loc[t.index]
            new_idx.loc[new_idx.isnull()] = t.index[new_idx.isnull()]
            t.index = new_idx
        res[col] = set(t.index)
    return res


def add_gene_symbols_to_ensembl_data(df, tax_id=9606):
    """
    Add gene symbols to the DataFrame df which is indexed by Ensembl IDs
    """
    gs = references.ensembl_to_gene_symbol(df.index, tax_id=tax_id)
    # resolve any duplicates arbitrarily (these should be rare)
    gs = gs.loc[~gs.index.duplicated()]
    df.insert(0, 'Gene Symbol', gs)


def add_fc_direction(df, logfc_field='logFC'):
    """
    Add direction column to DE data with the logFC in the field with name logfc_field
    """
    the_logfc = df.loc[:, logfc_field]
    direction = pd.Series(index=df.index, name='Direction')
    direction.loc[the_logfc < 0] = 'down'
    direction.loc[the_logfc > 0] = 'up'
    df.insert(df.shape[1], 'Direction', direction)


# These files are exported for a given reference from http://www.ensembl.org/biomart/martview
# Choose database 'Ensembl Genes 91', then desired reference, then export results with just two fields Gene Stable ID
# and Transcript Stable ID

GENE_TO_TRANSCRIPT_FILES = {
    9606: os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human',
        'ensembl',
        'GRCh38.p10.release90',
        'gene_to_transcript.txt'
    ),
    10090: os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'mouse',
        'ensembl',
        'GRCm38.p5.r90',
        'gene_to_transcript.txt'
    ),
}


def ensembl_transcript_quant_to_gene(dat, tax_id=9606, remove_ver=True):
    """
    Aggregate the supplied transcript-level quantification to gene level. Input index is Ensembl transcript ID.
    This is necessary for Salmon outputs, for example.
    :param dat: Pandas DataFrame containing transcript-level quantification
    :param tax_id: Default is human. Mouse is 10090.
    :param remove_ver: If True, remove accession version from dat.index (required for Salmon). Won't hurt if not needed.
    :return: Pandas DataFrame aggregated to gene level, indexed by Ensembl Gene ID.
    """
    if remove_ver:
        dat = dat.copy()
        dat.index = dat.index.str.replace(r'.[0-9]+$', '')

    if tax_id in GENE_TO_TRANSCRIPT_FILES:
        fn = GENE_TO_TRANSCRIPT_FILES[tax_id]
    else:
        raise AttributeError("Unsupported taxonomy ID %d" % tax_id)

    gene_transcript = pd.read_csv(fn, header=0, sep='\t').set_index('Transcript stable ID')

    # shouldn't be necessary, but remove transcripts that have no translation
    to_keep = dat.index.intersection(gene_transcript.index)
    if len(to_keep) != dat.shape[0]:
        to_drop = dat.index.difference(gene_transcript.loc[:, 'Transcript stable ID'])
        print "Discarding %d transcripts that have no associated gene: %s" % (
            len(to_drop), ', '.join(to_drop)
        )
        dat = dat.loc[to_keep]

    # gene list in same order as data
    genes = gene_transcript.loc[dat.index, 'Gene stable ID']

    return dat.groupby(genes).sum()
