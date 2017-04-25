from scipy.stats import chi2
from scipy.linalg import expm3, norm
import matplotlib.colors as colors
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import animation, cm
import numpy as np
import pandas as pd
from plotting.threed import plot_ellipsoid


def rotation(axis, theta):
    return expm3(np.cross(np.eye(3), axis/norm(axis)*theta))


def get_cmap(N, cmap='jet'):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap=cmap)
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color


def cluster_ellipse(data, p=0.99):
    """
    Compute the parameters required to plot a cluster ellipsoid, based on a multivariate normal assumption.
    :param data: Rows are observations, columns are PCA components
    :param p: The probability cutoff defining the centroid boundary.
    :return: (loc, width, height, angle)
    """
    chi2_cutoff = chi2.ppf(p, 2)  # 2nd arg is DOF

    loc = data.mean(axis=0)
    covar = np.cov(data.transpose())
    eigval, eigvec = np.linalg.eig(covar)
    ang = np.arctan2(eigvec[1, 0], eigvec[0, 0])  # y, x
    # define radii
    sx = np.sqrt(eigval[0] * chi2_cutoff)
    sy = np.sqrt(eigval[1] * chi2_cutoff)

    return loc, 2 * sx, 2 * sy, 180. * ang / np.pi


def cluster_ellipsoid(data, p=0.99):
    """
    Compute the parameters required to plot a cluster ellipsoid, based on a multivariate normal assumption.
    :param data: Rows are observations, columns are PCA components
    :param p: The probability cutoff defining the centroid boundary.
    :return: (loc, diameters, rotation matrix)
    """
    chi2_cutoff = chi2.ppf(p, 3)  # 2nd arg is DOF

    loc = data.mean(axis=0)
    covar = np.cov(data.transpose())
    eigval, eigvec = np.linalg.eig(covar)

    return loc, np.sqrt(eigval * chi2_cutoff), eigvec.transpose()


def autoscale(data, components, ax, **additional_data):
    """
    Apply custom autoscaling to the axes, based on the supplied data
    """
    nd = len(components)
    tmp_data = np.array(data, copy=True)
    for k, v in additional_data.items():
        tmp_data = np.concatenate((tmp_data, v), axis=0)
    ymax = np.ceil(tmp_data.max(axis=0) / 2. + 1) * 2.
    ymin = np.floor(tmp_data.min(axis=0) / 2. - 1) * 2.
    ax.set_xlim([ymin[components[0]], ymax[components[0]]])
    ax.set_ylim([ymin[components[1]], ymax[components[1]]])
    if nd == 3:
        ax.set_zlim([ymin[components[2]], ymax[components[2]]])


def pca_plot_by_group_2d(
        pca_data,
        subgroups,
        components=(0, 1),
        colour_map=None,
        marker_map=None,
        ellipses=True,
        ellipse_p=0.95,
        ax=None,
        legend=True,
        auto_scale=True,
        markersize=40,
        **additional_data
):
    """
    2D plot showing PCA respresentation of data, optionally overlaid with ellipses
    :param pca_data: Samples in rows, features in columns
    :param subgroups: The labelling for the pca_data
    :param components: An iterable of length 2 giving the 2 PCs to use
    :param colour_map: Optionally provide a dictionary, where key is the subgroup label and value is the plotting colour
    :param marker_map: Optionally provide a dictionary, where key is the subgroup label and value is the marker style
    :param ellipses: If True, plot ellipses
    :param ellipse_p: The 'probability' to include in the ellipse
    :param ax: If supplied, plot here, else create a new figure and axes
    :param legend: If True, add a legend
    :param auto_scale: If True, automatically set axis scaling
    :param additional_data: Dictionary containing optional addition data for the plot. If included, the key is the label
    used for the legend and the value is the same format as pca_data. To specify the colour, include the key in the
    colour_map, otherwise random colours will be selected
    :return:
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
    idx, labels = subgroups.factorize()

    if colour_map is None:
        cmap = get_cmap(len(labels) + len(additional_data))
        colour_map = dict([
            (k, cmap(i)) for i, k in enumerate(labels)
        ])
        for j, k in enumerate(additional_data.keys()):
            colour_map[k] = cmap(len(labels) + j)
    else:
        # check whether additional data have been included
        cmap = get_cmap(len(additional_data))
        for i, k in enumerate(additional_data.keys()):
            if k not in colour_map:
                colour_map[k] = cmap(i)

    if marker_map is None:
        marker_map = dict([(k, 'o') for k in labels])
        for k in additional_data:
            marker_map[k] = 'o'
    else:
        # check whether additional data have been included
        for k in additional_data:
            if k not in marker_map:
                marker_map[k] = 'o'

    for i, l in enumerate(labels):
        if l not in colour_map:
            continue
        c = colour_map[l]
        m = marker_map[l]
        j = idx == i
        ax.scatter(pca_data[j, components[0]], pca_data[j, components[1]], c=c, s=markersize, label=l, marker=m)
    for k, v in additional_data.items():
        c = colour_map[k]
        m = marker_map[k]
        # plot these on the top level
        ax.scatter(v[:, components[0]], v[:, components[1]], c=c, s=markersize, label=k, zorder=10, marker=m)
    if legend:
        plt.legend(frameon=False, loc='upper right')
    ax.set_xlabel("PCA component %s" % (components[0] + 1))
    ax.set_ylabel("PCA component %s" % (components[1] + 1))

    # (2D) ellipses
    if ellipses:
        for i, l in enumerate(labels):
            if l not in colour_map:
                continue
            c = colour_map[l]
            j = idx == i
            e = Ellipse(
                *cluster_ellipse(pca_data[j][:, components], p=ellipse_p),
                alpha=0.25,
                facecolor=c,
                edgecolor=c,
                lw=2
            )
            ax.add_artist(e)

    if auto_scale:
        autoscale(pca_data, components, ax, **additional_data)
    return ax


def pca_plot_by_group_3d(
    pca_data,
    subgroups,
    components=(0, 1, 2),
    colour_map=None,
    marker_map=None,
    plot_ellipsoids=True,
    ellipse_p=0.95,
    plot_pca_data=True,
    ax=None,
    markersize=40,
    auto_scale=True,
    **additional_data
):
    if ax is None:
        fig = plt.figure(figsize=(7.8, 6.6))
        ax = fig.add_subplot(111, projection='3d', aspect='equal')
    idx, labels = subgroups.factorize()
    if colour_map is None:
        cmap = get_cmap(len(labels) + len(additional_data))
        colour_map = dict([
            (k, cmap(i)) for i, k in enumerate(labels)
        ])
        for j, k in enumerate(additional_data.keys()):
            colour_map[k] = cmap(len(labels) + j)
    else:
        # check whether additional data have been included
        cmap = get_cmap(len(additional_data))
        for i, k in enumerate(additional_data.keys()):
            if k not in colour_map:
                colour_map[k] = cmap(i)

    if marker_map is None:
        marker_map = dict([(k, 'o') for k in labels])
        for k in additional_data:
            marker_map[k] = 'o'
    else:
        # check whether additional data have been included
        for k in additional_data:
            if k not in marker_map:
                marker_map[k] = 'o'

    if plot_pca_data:
        for i, l in enumerate(labels):
            if l not in colour_map:
                continue
            c = colour_map[l]
            m = marker_map[l]
            j = idx == i
            ax.scatter(
                pca_data[j, components[0]],
                pca_data[j, components[1]],
                pca_data[j, components[2]],
                c=c,
                s=markersize,
                marker=m,
                label=l
            )
    for k, v in additional_data.items():
        c = colour_map[k]
        m = marker_map[k]
        ax.scatter(
            v[:, components[0]],
            v[:, components[1]],
            v[:, components[2]],
            c=c,
            s=markersize,
            marker=m,
            label=k,
        )
    plt.legend(frameon=False, loc='upper right')

    if plot_ellipsoids:
        for i, l in enumerate(labels):
            if l not in colour_map:
                continue
            c = colour_map[l]
            j = idx == i
            loc, radii, rot = cluster_ellipsoid(pca_data[j][:, components], p=ellipse_p)
            plot_ellipsoid(loc, radii, rot, ax=ax, cageColor=c, cageAlpha=0.2, npt=100)

    if auto_scale:
        autoscale(pca_data, components, ax, **additional_data)

    ax.patch.set_facecolor('w')
    pane_c = (.8, .8, .8, 1.)
    pane_c_floor = (.8, .8, .8, .8)
    ax.w_xaxis.set_pane_color(pane_c)
    ax.w_yaxis.set_pane_color(pane_c)
    ax.w_zaxis.set_pane_color(pane_c_floor)
    ax.set_xlabel("PCA component %s" % (components[0] + 1))
    ax.set_ylabel("PCA component %s" % (components[1] + 1))
    ax.set_zlabel("PCA component %s" % (components[2] + 1))
    ax.view_init(65, 55)
    plt.tight_layout()

    return ax