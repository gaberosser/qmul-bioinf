import references

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
        gs = references.ensembl_to_gene_symbol(data.index).dropna()
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