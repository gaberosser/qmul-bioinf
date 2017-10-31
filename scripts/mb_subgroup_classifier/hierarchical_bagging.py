from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, fcluster
from scipy.spatial.distance import pdist
import pandas as pd
from load_data import microarray_data, allen_human_brain_atlas, rnaseq_data
from microarray import process
from matplotlib import pyplot as plt, gridspec
from matplotlib.colors import ListedColormap
import seaborn as sns
from settings import GIT_LFS_DATA_DIR
import numpy as np
import os
from scripts.comparison_rnaseq_microarray import consts


class Bagger(object):
    def __init__(self, data, axis=0):
        """
        Class for drawing random bootstrap samples with replacement
        :param data: The matrix of data to use
        :param axis: The axis to use when drawing random samples. axis=0 draws random rows.
        """
        self.data = data
        self.axis = axis
        self.m = self.data.shape[self.axis]

    def draw(self, n):
        """
        Draw random bootstrap sample of size n
        :return:
        """
        i = np.random.permutation(self.m)[:n]
        j = [slice(None)] * self.data.ndim
        j[self.axis] = i
        return (i, self.data[tuple(j)])

    def iterator(self, n=None):
        if n is None:
            # default behaviour: draw samples of same size as original data
            n = self.data.shape[self.axis]
        while True:
            yield self.draw(n)


def cluster_min_size(Z, min_count, min_n_cl):
    """
    Given the linkage matrix Z, find the flat clustering such that there are at least min_n_cl clusters
    with >= min_count members
    :param Z:
    :param min_count:
    :param min_n_cl:
    :return: (fcluster membership ID, cluster membership counts)
    """
    max_it = np.floor((Z.shape[0] + 1) / float(min_n_cl))
    i = 0
    cl = None
    cts = None
    while i < max_it:
        cl = fcluster(Z, min_count + i, criterion='maxclust')
        cts = np.histogram(cl, range(1, min_count + i + 2))[0]
        if (cts > min_n_cl).sum() == min_count:
            break
        else:
            i += 1
    return cl


if __name__ == '__main__':
    # ncott, ncott_meta = microarray_data.load_annotated_microarray_gse37382(
    #     aggr_field='SYMBOL',
    #     aggr_method='max'
    # )
    # sort_idx = ncott_meta.subgroup.sort_values().index
    # ncott_meta = ncott_meta.loc[sort_idx]
    # ncott = ncott.loc[:, sort_idx]

    # X = ncott.copy()
    # m = ncott_meta.copy()


    # load Kool dataset
    # kool, kool_meta = microarray_data.load_annotated_microarray_gse10327(
    #     aggr_field='SYMBOL',
    #     aggr_method='max',
    # )
    # sort_idx = kool_meta.subgroup.sort_values().index
    # kool_meta = kool_meta.loc[sort_idx]
    # kool = kool.loc[:, sort_idx]
    # kool_meta.loc[:, 'subgroup'] = (
    #     kool_meta.loc[:, 'subgroup'].str
    #         .replace('A', 'WNT')
    #         .replace('B', 'SHH')
    #         .replace('E', 'Group 3')
    #         .replace('C', 'Group 4')
    #         .replace('D', 'Group 4')
    # )
    # X = kool.copy()
    # m = kool_meta.copy()

    # load Robinson dataset (probe set level)
    robi, robi_meta = microarray_data.load_annotated_microarray_gse37418()
    # robi_meta = robi_meta.loc[~robi_meta.subgroup.isin(['U', 'SHH OUTLIER'])]
    sort_idx = robi_meta.subgroup.sort_values().index
    robi_meta = robi_meta.loc[sort_idx]
    robi = robi.loc[:, sort_idx]
    # robi_meta.loc[:, 'subgroup'] = robi_meta.subgroup.str.replace('G3', 'Group 3').replace('G4', 'Group 4')
    X = robi.copy()
    m = robi_meta.copy()

    # the samples must be represented in ROWS
    X = X.transpose()
    b = Bagger(X.values)
    it = b.iterator()

    n_it = 1000
    members = []
    res = []
    for i in range(n_it):
        ii, x = it.next()
        Z = linkage(x, method='average', metric='euclidean')
        res.append(cluster_min_size(Z, 4, 5))
        members.append(ii)


    # Z = linkage(X, method='average', metric='correlation')
    # clusters = []
    # counts = []
    #
    # k_min = 2
    # k_max = 5
    # for k in range(k_min, k_max + 1):
    #     i = 0
    #     THRESH = 5
    #     cl, cts = cluster_min_size(Z, k, THRESH)
    #     clusters.append(cl)
    #     counts.append(cts)
