from scipy.stats import chi2
from scipy.linalg import expm, norm
import matplotlib.colors as colors
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import animation, cm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
from plotting.threed import plot_ellipsoid


def rotation(axis, theta):
    # was using expm3, but this has been removed in a more recent version of Scipy
    # check that this still performs as expected (if anything, it should be more accurate)
    return expm(np.cross(np.eye(3), axis/norm(axis)*theta))


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
        data_supplied_as='pca',
        variance_explained=None,
        components=(0, 1),
        colour_map=None,
        marker_subgroups=None,
        marker_map=None,
        ellipses=True,
        ellipse_p=0.95,
        ax=None,
        legend=True,
        auto_scale=True,
        markersize=40,
        additional_data_dict=None,
        **additional_data
):
    """
    2D plot showing PCA respresentation of data, optionally overlaid with ellipses
    :param pca_data: Samples in rows, features in columns
    :param subgroups: The labelling for the pca_data
    :param data_supplied_as: Either 'pca' or 'raw'. If the latter, run PCA fitting first.
    :param components: An iterable of length 2 giving the 2 PCs to use
    :param colour_map: Optionally provide a dictionary, where key is the subgroup label and value is the plotting colour
    :param marker_map: Optionally provide a dictionary, where key is the subgroup label and value is the marker style
    TODO: The way subgroups is implemented here is a bit messy and could do with a refactor!
    :param ellipses: If True, plot ellipses
    :param ellipse_p: The 'probability' to include in the ellipse
    :param ax: If supplied, plot here, else create a new figure and axes
    :param legend: If True, add a legend
    :param auto_scale: If True, automatically set axis scaling
    :param variance_explained: If provided, this gives the % variance explained for each component. If None (default),
    do not plot, unless raw data are provided.
    :param additional_data_dict, additional_data: Dictionary containing optional addition data for the plot.
     Either specify as kwargs (**additional_data) or a single dict. Only one may be supplied.
    If included, the key is the label used for the legend and the value is the same format as pca_data.
    To specify the colour, include the key in the colour_map, otherwise random colours will be selected.
    :return:
    """
    if data_supplied_as != 'pca':
        p = PCA()
        pca_data = p.fit_transform(pca_data.transpose())
        variance_explained = p.explained_variance_ratio_ * 100.

    if additional_data_dict is not None:
        if len(additional_data) > 0:
            raise AttributeError("Must supply EITHER additional_data kwargs OR additional_data_dict, but not both.")
        additional_data = additional_data_dict

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

    if variance_explained is None:
        ax.set_xlabel("PCA component %s" % (components[0] + 1))
        ax.set_ylabel("PCA component %s" % (components[1] + 1))
    else:
        ax.set_xlabel("PCA component %s (%.1f%%)" % (components[0] + 1, variance_explained[components[0]]))
        ax.set_ylabel("PCA component %s (%.1f%%)" % (components[1] + 1, variance_explained[components[1]]))

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


def highlight_biplot_features(
    feat_x,
    feat_y,
    radius,
    ax
):
    """
    Highlight all features that are outside of the specified radius on the biplot
    :param feat_x:
    :param feat_y:
    :param radius:
    :param ax:
    :return:
    """
    selected = (feat_x ** 2 + feat_y ** 2) ** .5 > radius
    ax.scatter(
        feat_x[selected],
        feat_y[selected],
        facecolor='none',
        edgecolor='b',
        linewidths=1.,
        s=12.
    )
    return selected



def biplot(
    data,
    feat_axis=0,
    scale=1.,
    plot_dims=(0, 1),
    preserve_distance='samples',
    include_weighting=True,
    sample_colours=None,
    sample_markers=None,
    highlight_feature_radius=None
):
    """

    :param data:
    :param feat_axis:
    :param scale:
    :param include_weighting:
    :param sample_colours:
    :param sample_markers:
    :param highlight_feature_radius:
    :return:
    """
    if preserve_distance not in ('samples', 'features'):
        raise ValueError("preserve_distance parameter must be 'samples' or 'features'.")
    if feat_axis not in (0, 1):
        raise ValueError("Input data must be a 2D matrix, and feat_axis must be 0 or 1.")

    data = pd.DataFrame(data, copy=False)
    if feat_axis == 1:
        data = data.transpose()

    if sample_markers is None:
        sample_markers = dict([(t, 'o') for t in data.columns])

    # group sample markers, so we can run the minimum number of calls to plt.scatter()
    sample_markers_grouped = dict()
    sample_markers_grouped_ix = dict()
    for k, v in sample_markers.items():
        sample_markers_grouped.setdefault(v, []).append(k)
        sample_markers_grouped_ix.setdefault(v, []).append(np.where(data.columns == k)[0][0])

    # no need for the colours, because they can be passed in as a vector
    if sample_colours is None:
        sample_colours = dict([(t, 'k') for t in data.columns])

    n = float(data.shape[1])
    p = float(data.shape[0])

    # standardise: mean centred data required for sensible decomposition
    # standardisation occurs along the FEATURES axis, which is dim 1
    scaler = StandardScaler(with_std=False)

    # features on the ROWS, mean centre by gene
    scaler = scaler.fit(data.transpose())
    X = scaler.transform(data.transpose()).transpose()

    # SVD
    u, s, vh = np.linalg.svd(X, full_matrices=False)

    # checked this against the sklearn PCA code
    explained_variance = (s ** 2) / n
    explained_variance_ratio = explained_variance / explained_variance.sum()

    if preserve_distance == 'samples':
        # preserve inter-sample distances

        # project gene data into PCA
        # this matches the output of pca.transform() (except for possible sign switch)
        if include_weighting:
            # scaling by s: components scale by their relative explanatory power (not linearly)
            # the plot may appear 'squashed', depending on the weighting
            us = u.dot(np.diag(s))
            feat_x = scale * us[:, plot_dims[0]]
            feat_y = scale * us[:, plot_dims[1]]
        else:
            # alternatively, just plot the unscaled feature components (plot becomes more circular)
            feat_x = scale * u[:, plot_dims[0]]
            feat_y = scale * u[:, plot_dims[1]]

        sample_x = vh[plot_dims[0]]
        sample_y = vh[plot_dims[1]]
    else:
        # preserve inter-feature distances
        ## TODO: check this

        if include_weighting:
            # scaling by s: components scale by their relative explanatory power (not linearly)
            # the plot may appear 'squashed', depending on the weighting
            vs = vh.dot(np.diag(s))
            sample_x = vs[plot_dims[0]]
            sample_y = vs[plot_dims[1]]
        else:
            # alternatively, just plot the unscaled feature components (plot becomes more circular)
            sample_x = vh[plot_dims[0]]
            sample_y = vh[plot_dims[1]]

        feat_x = scale * u[:, plot_dims]
        feat_y = scale * u[:, plot_dims]

    # track legend entries to avoid double labelling
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # scatter features
    ax.scatter(feat_x, feat_y, c='gray', s=10, edgecolor='none', alpha=0.7)

    # scatter samples
    for mrk, arr in sample_markers_grouped.items():
        cs = [sample_colours[k] for k in arr]
        ax.scatter(
            sample_x[sample_markers_grouped_ix[mrk]],
            sample_y[sample_markers_grouped_ix[mrk]],
            c=cs,
            marker=mrk,
            edgecolor='k',
            linewidth=1.,
            zorder=10.
        )

    if highlight_feature_radius is not None:
        selected = highlight_biplot_features(feat_x, feat_y, highlight_feature_radius, ax=ax)

    ax.set_xlabel('PC%d (%.2f %%)' % (plot_dims[0] + 1, explained_variance_ratio[plot_dims[0]] * 100.))
    ax.set_ylabel('PC%d (%.2f %%)' % (plot_dims[1] + 1, explained_variance_ratio[plot_dims[1]] * 100.))

    return {
        'u': u,
        's': s,
        'vh': vh,
        'fig': fig,
        'ax': ax,
        'explained_variance': explained_variance,
        'explained_variance_ratio': explained_variance_ratio,
        'sample_data': (sample_x, sample_y),
        'feature_data': (feat_x, feat_y),
        'components': plot_dims
    }