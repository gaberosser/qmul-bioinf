from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd


def svd_for_biplot(data, feat_axis=0, preserve_distance='samples', include_weighting=True, scale_preserved=1.):
    """
    Apply SVD to the supplied data, then extract the feature and sample representation that we require for a biplot.
    This is essentially just U and S \dot Vh (or vice versa, depending on the feature axis and in which axis we choose
    to preserve distances.
    Example: dat is a (p x n) matrix containing gene expression values from p genes measured over n samples.
    :param dat: Matrix containing multiple observations (samples) of a number of features.
    :param feat_axis: The axis corresponding to features.
    :param preserve_distance: Either 'samples' or 'features', this specifies which axis we are preserving. The other
    will be converted into loadings, which do NOT have preserved pairwise distances (?).
    :param include_weighting: If True (default), scale the preserved axis by S (the diagonal matrix). If this is False,
    pairwise distances will NOT be preserved.
    :param scale_preserved: Scalar parameter used to rescale the axis we are preserving.
    :return:
    """
    # ensure that features are on the rows
    if feat_axis == 1:
        data = data.transpose()

    n = float(data.shape[1])

    # standardise: mean centred data required for sensible decomposition
    # standardisation occurs along the FEATURES axis, which is dim 1
    scaler = StandardScaler(with_std=False)

    # features on the ROWS, mean centre by gene
    scaler = scaler.fit(data.transpose())
    X = scaler.transform(data.transpose()).transpose()

    # SVD
    u, s, vh = np.linalg.svd(X, full_matrices=False)

    # for info: checked this against the sklearn PCA code
    # n = float(data.shape[1])
    # explained_variance = (s ** 2) / n
    # explained_variance_ratio = explained_variance / explained_variance.sum()

    if preserve_distance == 'samples':
        # preserve inter-sample distances

        # project gene data into PCA
        # this matches the output of pca.transform() (except for possible sign switch)
        if include_weighting:
            # scaling by s: components scale by their relative explanatory power (not linearly)
            # the plot may appear 'squashed', depending on the weighting
            ss = np.diag(s)
        else:
            ss = np.eye(len(s))

        us = u.dot(ss)

        feat_dat = dict([
            (i + 1, scale_preserved * us[:, i]) for i in range(us.shape[1])
        ])
        sample_dat = dict([
            (i + 1, vh[i]) for i in range(vh.shape[0])
        ])
        # feat_dat = scale_preserved * us
        # sample_dat = vh
    else:
        # preserve inter-feature distances

        if include_weighting:
            ss = np.diag(s)
        else:
            ss = np.eye(len(s))
        # scaling by s: components scale by their relative explanatory power (not linearly)
        # the plot may appear 'squashed', depending on the weighting

        ## TODO: check this

        vs = vh.dot(np.diag(ss))

        sample_dat = dict([
            (i + 1, u[:, i]) for i in range(u.shape[1])
        ])
        feat_dat = dict([
            (i + 1, scale_preserved * vs[i]) for i in range(vs.shape[0])
        ])
        # sample_dat = u
        # feat_dat = scale_preserved * vs

    explained_variance = (s ** 2) / n
    explained_variance_ratio = explained_variance / explained_variance.sum()

    feat_dat = pd.DataFrame(feat_dat, index=data.index)
    sample_dat = pd.DataFrame(sample_dat, index=data.columns)

    return {
        'u': u,
        's': s,
        'vh': vh,
        'feat_dat': feat_dat,
        'sample_dat': sample_dat,
        'explained_variance': explained_variance,
        'explained_variance_ratio': explained_variance_ratio,
    }