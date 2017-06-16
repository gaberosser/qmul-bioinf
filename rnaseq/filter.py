import pandas as pd


def filter_by_cpm(data, min_cpm=1, min_n_samples=3, unless_cpm_gt=None):
    """
    Filter the features in the supplied data (rows=features, cols=samples) by counts per million.
    :param min_n_samples: The minimum number of samples that must meet the criteria to have the feature accepted
    :param unless_cpm_gt: If any single reading exceeds this CPM, it is retained (optional)
    :return: New dataframe with (potentially) reduced number of features
    """
    n = data.sum(axis=0)
    cpm = data.divide(n, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= min_n_samples
    if unless_cpm_gt is not None:
        keep = keep | (cpm >= unless_cpm_gt).any(axis=1)
    return data.loc[keep]
