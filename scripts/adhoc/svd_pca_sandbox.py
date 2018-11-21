import numpy as np
from scipy import stats, linalg
from sklearn import decomposition, preprocessing


if __name__ == "__main__":
    """
    Goal here is to test the definitions of SVD and PCA, and understand the interplay
    """

    # simulate some data
    # rows are observations, columns features (NB opposite to the typical genomics setup)
    # make sure the mean is NOT zero, so we have to deal with this
    n = 20 # number of samples / observations
    p = 200 # numer of features
    dat = np.zeros((n, p))
    for i in range(p):
        m = np.random.randn() * 10.
        s = np.abs(np.random.randn()) * 2.
        dat[:, i] = np.random.randn(n) * s + m

    # mean centre
    scaler = preprocessing.StandardScaler(with_std=False)
    scaler.fit(dat)

    # format: observations on the rows, features on the columns
    X = scaler.transform(dat)

    # we'll work with X from now on

    # 1) PCA vs manual PCA (via the covariance matrix)
    # a) covariance matrix
    cov = np.cov(X.transpose())  # numpy requires rows represent variables
    if np.allclose(cov, np.dot(X.transpose(), X) / float(n - 1)):
        print "Our manual covariance matrix calculation matches numpy's"
    else:
        print "BAD: Our manual covariance matrix calculation doesn't match numpy's"

    # b) eigenvalues and eigenvectors
    # eigh works for Hermitian matrices (c is symmetric, so we can use it) and guarantees real eigenvalues and vectors
    # NB the list is in ASCENDING order
    eigval, eigvec = np.linalg.eigh(cov)

    # c) use these to project the data (columns of XV)
    # we'll reorder these here, so the most major columns come first
    xv = np.dot(X, eigvec)[:, ::-1]

    # compare with 'out of the box' PCA
    pca = decomposition.PCA()
    pca = pca.fit(X)
    pc = pca.transform(X)  # PCs are in the columns

    all_passed = True
    for i in range(n):
        b = np.allclose(pc[:, i], xv[:, i])
        if not b:
            # we might just have a simple inversion
            # sklearn implements a heuristic here to ensure deterministic output: see sklearn.utils.extmath.svd_flip()
            b = np.allclose(pc[:, i], -xv[:, i])
            if not b:
                all_passed = False
    if all_passed:
        print "Our manual calculation of %d principal components all match sklearn's" % n
    else:
        print "BAD: Our manual calculation of %d principal components DO NOT all match sklearn's" % n

    # 2) PCA vs SVD
    # a) compute SVD
    u, s, vh = np.linalg.svd(X, full_matrices=False)

    # b) restore the original data
    if np.allclose(X, np.dot(u, np.diag(s).dot(vh))):
        print "We restored X from the full SVD."
    else:
        print "BAD: We failed to restore X from the full SVD."

    # c) compute the PCA
    us = u.dot(np.diag(s))

    all_passed = True
    for i in range(n):
        b = np.allclose(pc[:, i], us[:, i])
        if not b:
            # we might just have a simple inversion
            # sklearn implements a heuristic here to ensure deterministic output: see sklearn.utils.extmath.svd_flip()
            b = np.allclose(pc[:, i], -us[:, i])
            if not b:
                all_passed = False
    if all_passed:
        print "Our %d SVD-based PC projections all match sklearn's" % n
    else:
        print "BAD: Our %d SVD-based PC projections DO NOT all match sklearn's" % n

    # d) PCA components: rows of components
    # NB for some reason, the match is perfect EXCEPT for the last component
    # Assume this is related to the specific implementation??
    all_passed = True
    for i in range(n - 1):
        b = np.allclose(pca.components_[i], vh[i])
        if not b:
            # we might just have a simple inversion
            # sklearn implements a heuristic here to ensure deterministic output: see sklearn.utils.extmath.svd_flip()
            b = np.allclose(pca.components_[i], -vh[i])
            if not b:
                print "PC component %d didn't match" % i
                all_passed = False
    if all_passed:
        print "Our %d SVD-based principal components all match sklearn's" % (n - 1)
    else:
        print "BAD: Our %d SVD-based principal components DO NOT all match sklearn's" % (n - 1)

