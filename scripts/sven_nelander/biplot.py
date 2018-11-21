import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import pandas as pd

from rnaseq import loader, filter
from scripts.hgic_final import consts

if __name__ == "__main__":
    """
    Idea here: recreate the analysis Sven carried out, generating biplots for the RNA-Seq data.
    We can then extend this idea to methylation (?)

    Absolutely amazing resource here:
    https://www.fbbva.es/microsite/multivariate-statistics/biplots.html

    Code snippet inspiration here:
    https://stackoverflow.com/questions/39216897/plot-pca-loadings-and-loading-in-biplot-in-sklearn-like-rs-autoplot
    """
    eps = 1.  # offset applied during log transform

    # load data for iNSC and GBM (Salmon TPM)
    obj = loader.load_by_patient(consts.PIDS, source='salmon')
    obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES_GIC + consts.S1_RNASEQ_SAMPLES_INSC)

    scatter_colours = {
        'GBM': 'red',
        'iNSC': 'blue'
    }

    dat = filter.filter_by_cpm(obj.data, min_n_samples=2)
    # TODO: include VST or similar here
    dat = np.log(dat + eps)

    n = float(dat.shape[1])

    # standardise: mean centred data required for sensible decomposition
    # standardisation occurs along the FEATURES axis, which is dim 1
    scaler = StandardScaler(with_std=False)

    # features on the COLS, mean centre by PATIENT
    # scaler.fit(dat)
    # X = scaler.transform(dat).transpose()

    # features on the ROWS, mean centre by PATIENT
    # scaler.fit(dat)
    # X = scaler.transform(dat)

    # features on the ROWS, mean centre by gene
    scaler = scaler.fit(dat.transpose())
    X = scaler.transform(dat.transpose()).transpose()

    # features on the COLS, mean centre by gene
    # scaler.fit(dat.transpose())
    # X = scaler.transform(dat.transpose())

    # SVD
    u, s, vh = np.linalg.svd(X, full_matrices=False)

    preserve_dist_ax = 'samples'
    dims = (0, 1)  # dimensions / components to include in plot
    scale = .1

    if preserve_dist_ax == 'samples':
        # preserve inter-sample distances

        # project gene data into PCA
        # this matches the output of pca.transform() (except for possible sign switch)
        us = u.dot(np.diag(s))

        feat_x = scale * us[:, dims[0]]
        feat_y = scale * us[:, dims[1]]

        sample_x = vh[dims[0]]
        sample_y = vh[dims[1]]
    else:
        ## FIXME: finish this
        # preserve inter-feature distances
        xx = u[:,dims[0]]
        yy = u[:, dims[1]]

        vx = scale * s[dims[0]] * vh[dims[0]]
        vy = scale * s[dims[1]] * vh[dims[1]]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(feat_x, feat_y, c='gray', s=10, edgecolor='none', alpha=0.5)
    arr_ix, typ = obj.meta.type.factorize()
    for i, t in enumerate(typ):
        ax.scatter(sample_x[arr_ix == i], sample_y[arr_ix == i], c=scatter_colours[t], edgecolor='k', linewidth=1, label=t)
    # [ax.arrow(0, 0, dx, dy) for dx, dy in zip(vx, vy)]


    pass