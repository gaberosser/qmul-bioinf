import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_ellipsoid(loc, radii, rot, ax=None, plotAxes=False, cageColor='b', cageAlpha=0.2, npt=100):
    """
    Modified from original at: https://github.com/minillinim/ellipsoid/blob/master/ellipsoid.py
    Plot an ellipsoid
    """
    make_ax = ax == None
    if make_ax:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    u = np.linspace(0.0, 2.0 * np.pi, npt)
    v = np.linspace(0.0, np.pi, npt)

    # cartesian coordinates that correspond to the spherical angles:
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

    # rotate accordingly
    for i in range(len(x)):
        for j in range(len(x)):
            [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]], rot) + loc

    if plotAxes:
        # make some purdy axes
        axes = np.array([[radii[0], 0.0, 0.0],
                         [0.0, radii[1], 0.0],
                         [0.0, 0.0, radii[2]]])
        # rotate accordingly
        for i in range(len(axes)):
            axes[i] = np.dot(axes[i], rot)

        # plot axes
        for p in axes:
            X3 = np.linspace(-p[0], p[0], 100) + loc[0]
            Y3 = np.linspace(-p[1], p[1], 100) + loc[1]
            Z3 = np.linspace(-p[2], p[2], 100) + loc[2]
            ax.plot(X3, Y3, Z3, color=cageColor)

    # plot ellipsoid
    ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color=cageColor, alpha=cageAlpha)

    return ax