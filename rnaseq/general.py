import references
import pandas as pd


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