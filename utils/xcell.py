import pandas as pd


def load_raw_results(
    score_fn,
    pval_fn,
    raw_score_fn=None
):
    """
    Load the raw results returned by the xCell server.
    :param score_fn: Path to the text file contianing normalised scores.
    :param pval_fn: Path to the text file contianing p values.
    :param raw_score_fn: (Optional) Path to the text file contianing raw scores.
    :return: Dictionary of DataFrames, which can be exported directly to Excel if required.
    """
    res = {}
    res['Proportions'] = pd.read_csv(score_fn, sep='\t', header=0, index_col=0)
    res['P values'] = pd.read_csv(pval_fn, sep='\t', header=0, index_col=0)
    res['P values'].columns = res['Proportions'].columns
    if raw_score_fn:
        res['Raw scores'] = pd.read_csv(raw_score_fn, sep='\t', header=0, index_col=0)
    return res