from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
from load_data import microarray_data, allen_human_brain_atlas, rnaseq_data
from microarray import process
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.colors as colors
import seaborn as sns
from scripts.comparison_rnaseq_microarray import consts
from scripts.mb_subgroup_classifier.load import load_xz_rnaseq
from scripts.mb_subgroup_classifier.pca_exploratory import rotation, plot_ellipsoid, cluster_ellipsoid
from scipy.stats import chi2, multivariate_normal
from scipy.linalg import expm3, norm

plt.close('all')

# simulate some random data
covar = np.array([[ 4.86448242,  2.72598167,  1.7020008 ],
       [ 2.72598167,  2.93005138,  1.10588911],
       [ 1.7020008 ,  1.10588911,  5.46018852]])

data = multivariate_normal.rvs(cov=covar, size=100)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data[:, 0], data[:, 1], data[:, 2])
ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')

p = 0.9
chi2_cutoff = chi2.ppf(p, 3)  # 2nd arg is DOF

loc = data.mean(axis=0)
eigval, eigvec = np.linalg.eig(covar)
radii = np.sqrt(eigval * chi2_cutoff)

# plot the eigenvectors
for i in range(3):
    ax.plot3D(
        [0, eigvec[0, i] * eigval[i]],
        [0, eigvec[1, i] * eigval[i]],
        [0, eigvec[2, i] * eigval[i]],
        'k-'
    )

plot_ellipsoid(loc, radii, eigvec, ax=ax, cageAlpha=0.2, cageColor='g')

npt = 100
u = np.linspace(0.0, 2.0 * np.pi, npt)
v = np.linspace(0.0, np.pi, npt)

# cartesian coordinates that correspond to the spherical angles:
# x = radii[0] * np.outer(np.cos(u), np.sin(v))
# y = radii[1] * np.outer(np.sin(u), np.sin(v))
# z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
#
# # rotate accordingly
# for i in range(len(x)):
#     for j in range(len(x)):
#         [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]], eigvec) + loc
#
# ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color='b', alpha=0.3)

# start with the ellipsoid aligned along the axes
# plotEllipsoid(loc, radii, np.eye(3), ax=ax, cageAlpha=0.2, cageColor='y')
#
# # we now need to find the rotation matrix for the first eigenvector
# x, y, z = eigvec[:, 0]
#
# # two rotation operations: phi and theta
# # phi: rotate around z
# # TODO: why NEGATIVE rotation??
# phi = np.arctan2(y, x)
# rot1 = rotation(np.array([0, 0, -1]), phi)
#
# # plot it
# plotEllipsoid(loc, radii, rot1, ax=ax, cageAlpha=0.2, cageColor='g')
#
# # plot major ax
#
# # now rotate by theta around an orthogonal axis
# orth = np.dot([0, 1, 0], rot1)
# ax.plot3D(
#     [0, orth[0] * eigval[0]],
#     [0, orth[1] * eigval[0]],
#     [0, orth[2] * eigval[0]],
#     c='r',
#     lw=3.5
# )
#
# # theta: rotate around x (or negative y)
# # theta = np.arccos(z / (x ** 2 + y ** 2) ** 0.5)
# theta = np.arctan2((x ** 2 + y ** 2) ** 0.5, z)
# rot2 = rotation(orth, theta)
#
# rot = np.dot(rot1, rot2)
#
# plotEllipsoid(loc, radii, rot, ax=ax, cageAlpha=0.2, cageColor='b')

plt.show()