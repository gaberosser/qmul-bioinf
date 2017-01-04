from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import os
from load_data import microarray_data, allen_human_brain_atlas, rnaseq_data
from microarray import process
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation, cm
import matplotlib.colors as colors
import seaborn as sns
from scripts.comparison_rnaseq_microarray import consts
from scripts.mb_subgroup_classifier.load import load_xz_rnaseq, load_xiaonan_microarray
from scipy.stats import chi2
from scipy.linalg import expm3, norm


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


def combine_expr(*args):
    common_genes = args[0].index
    sample_names = args[0].columns
    for t in args[1:]:
        common_genes = common_genes.intersection(t.index)
        sample_names = sample_names.append(t.columns)
    res = pd.DataFrame(index=common_genes, columns=sample_names)
    for t in args:
        res.loc[:, t.columns] = t.loc[common_genes]
    return res


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


def animation_func(ax, n_frames):
    # fixed elevation, vary azimuth
    elev = ax.elev
    def animate(i):
        ax.view_init(elev, i * (360. / float(n_frames)))
    return animate


if __name__ == '__main__':

    ELL_P = 0.95
    N_PC = 3
    ANIM = True
    # ANIM = False
    ANIM_NFRAMES = 300
    SAVEFIG = True
    # SAVEFIG = False

    if SAVEFIG:
        OUTDIR = 'pca_results.tmp.0'
        i = 1
        while os.path.exists(OUTDIR):
            OUTDIR = 'pca_results.tmp.%d' % i
            i += 1
        print "Creating temp output dir %s" % OUTDIR
        os.makedirs(OUTDIR)


    def plot_2d(y, lbl, colour_map, marker_map, title=None, **additional_data):
        # plots: PCA of classifier vs RNA-Seq
        fig, axs = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=False, figsize=(12, 5))
        for i, compo in enumerate([(0, 1), (1, 2), (0, 2)]):
            pca_plot_by_group_2d(
                y, lbl, components=compo, colour_map=colour_map,
                marker_map=marker_map,
                ellipse_p=ELL_P,
                ax=axs[i],
                legend=False,
                **additional_data
            )
        axs[-1].legend(loc='upper right')
        plt.tight_layout(pad=0.2, rect=[.02, .02, 1, 1])
        if title:
            fig.savefig(os.path.join(OUTDIR, "%s.png" % title), dpi=300)
            fig.savefig(os.path.join(OUTDIR, "%s.tiff" % title), dpi=200)
            fig.savefig(os.path.join(OUTDIR, "%s.pdf" % title))

        return fig, axs


    def plot_3d(
            y,
            lbl,
            colour_map,
            marker_map,
            title=None,
            anim=ANIM,
            n_frames=ANIM_NFRAMES,
            plot_y=True,
            **additional_data
    ):
        ax_3d = pca_plot_by_group_3d(
            y, lbl, colour_map=colour_map, ellipse_p=ELL_P,
            marker_map=marker_map, plot_pca_data=plot_y,
            **additional_data
        )
        ax_3d.view_init(35, 55)
        fig = ax_3d.get_figure()
        if title:
            fig.savefig(os.path.join(OUTDIR, "%s.png" % title), dpi=300)
            fig.savefig(os.path.join(OUTDIR, "%s.tiff" % title), dpi=200)
            fig.savefig(os.path.join(OUTDIR, "%s.pdf" % title))

        if anim:
            if not title:
                raise ValueError("If anim is True, must supply a valid title.")
            ax_3d.set_aspect('equal')
            animate = animation_func(ax_3d, n_frames)
            anim = animation.FuncAnimation(fig, animate, frames=n_frames, interval=20, repeat=False)
            anim.save(os.path.join(OUTDIR, '%s.mp4' % ttl), fps=30, writer='mencoder',
                      extra_args=['-ovc', 'lavc', '-lavcopts', 'vcodec=mpeg4:vbitrate=2400'])
            plt.close(fig)

        return fig, ax_3d


    # it's useful to maintain a list of known upregulated genes
    nano_genes = []
    for grp, arr in consts.NANOSTRING_GENES:
        if grp != 'WNT':
            nano_genes.extend(arr)
    nano_genes.remove('EGFL11')
    nano_genes.append('EYS')

    # load Ncott data (285 non-WNT MB samples)
    ncott, ncott_meta = microarray_data.load_annotated_microarray_gse37382(
        aggr_field='SYMBOL',
        aggr_method='max'
    )
    sort_idx = ncott_meta.subgroup.sort_values().index
    ncott_meta = ncott_meta.loc[sort_idx]
    ncott = ncott.loc[:, sort_idx]
    ncott = process.yugene_transform(ncott)

    # load Allen (healthy cerebellum)

    he, he_meta = allen_human_brain_atlas.cerebellum_microarray_reference_data(agg_field='gene_symbol', agg_method='max')
    he_meta.loc[:, 'subgroup'] = 'control'

    # load Kool dataset
    kool, kool_meta = microarray_data.load_annotated_microarray_gse10327(
        aggr_field='SYMBOL',
        aggr_method='max',
    )
    sort_idx = kool_meta.subgroup.sort_values().index
    kool_meta = kool_meta.loc[sort_idx]
    kool = kool.loc[:, sort_idx]
    kool_meta.loc[:, 'subgroup'] = (
        kool_meta.loc[:, 'subgroup'].str
            .replace('A', 'WNT')
            .replace('B', 'SHH')
            .replace('E', 'Group 3')
            .replace('C', 'Group 4')
            .replace('D', 'Group 4')
    )
    kool = process.yugene_transform(kool)

    # load Robinson dataset
    robi, robi_meta = microarray_data.load_annotated_microarray_gse37418(aggr_field='SYMBOL', aggr_method='max')
    robi_meta = robi_meta.loc[~robi_meta.subgroup.isin(['U', 'SHH OUTLIER'])]
    sort_idx = robi_meta.subgroup.sort_values().index
    robi_meta = robi_meta.loc[sort_idx]
    robi = robi.loc[:, sort_idx]
    robi_meta.loc[:, 'subgroup'] = robi_meta.subgroup.str.replace('G3', 'Group 3').replace('G4', 'Group 4')
    robi = process.yugene_transform(robi)

    # combine them all - this ensures a common list of genes
    combs = []
    grps = []
    # combs.append(he); grps.extend(he_meta.subgroup.values)
    combs.append(ncott); grps.extend(ncott_meta.subgroup.values)
    combs.append(kool); grps.extend(kool_meta.subgroup.values)
    combs.append(robi); grps.extend(robi_meta.subgroup.values)

    all_expr = combine_expr(*combs)

    ###################################
    # extract only the desired samples
    # can switch classifier here.
    ###################################
    title = 'robi'
    X = all_expr.loc[:, robi_meta.index].transpose()
    m = robi_meta.copy()

    # title = 'kool'
    # X = all_expr.loc[:, kool_meta.index].transpose()
    # m = kool_meta.copy()

    # title = 'northcott'
    # X = all_expr.loc[:, ncott_meta.index].transpose()
    # m = ncott_meta.copy()

    # first fit the pca with lots of components to investigate the explained variance
    pca = PCA(n_components=10)
    pca.fit(X)
    y = pca.transform(X)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(range(1, 11), np.cumsum(pca.explained_variance_), '-o')
    plt.axis([0.8, 10.2, 0, 100.])
    ax.set_xticks(range(1, 11))
    plt.xlabel('Principal component')
    plt.ylabel('Cumulative % variance explained')
    if SAVEFIG:
        fig.savefig(os.path.join(OUTDIR, "pca_variance_explained.png"), dpi=200)
        fig.savefig(os.path.join(OUTDIR, "pca_variance_explained.pdf"), dpi=200)


    # now fit again, this time with only the necessary number of components
    # PCA fitting requires sample to be represented in the ROWS
    pca = PCA(n_components=N_PC)
    pca.fit(X)
    y = pca.transform(X)

    # load RNA-Seq data
    X_htseq = load_xz_rnaseq(kind='htseq', yugene=True, gene_symbols=X.columns).transpose()
    X_cuff = load_xz_rnaseq(kind='cuff', yugene=True, gene_symbols=X.columns).transpose()
    y_htseq = pca.transform(X_htseq)
    y_cuff = pca.transform(X_cuff)

    # load Xiao-Nan data
    xnan_sample_names = (
        'ICb1299-III',
        'ICb1299-IV',
        'ICb1487-I',
        'ICb1487-III',
        'ICb1595-I',
        'ICb1595-III',
    )
    xnan, xnan_meta = load_xiaonan_microarray(yugene=True, gene_symbols=X.columns, sample_names=xnan_sample_names)
    xnan = xnan.transpose()
    y_xnan = pca.transform(xnan)

    y_1299 = y_xnan[0:2]
    y_1487 = y_xnan[2:4]
    y_1595 = y_xnan[4:]

    # get labels and plot by subgroup
    idx, labels = m.subgroup.factorize()

    # define colours and labels
    late_pass_lbl = 'ICb1299 (this study)'
    lbl_1299 = 'ICb1299 (Zhao et al.)'
    lbl_1487 = 'ICb1487 (Zhao et al.)'
    lbl_1595 = 'ICb1595 (Zhao et al.)'

    colour_map = {
        'Group 3': '#F2EA00',
        'Group 4': '#2A8C43',
        'WNT': '#2D438E',
        'SHH': '#E5161C',
        'control': 'gray',
        'Late passage': 'k',
        late_pass_lbl: 'k',
        'Early passage': 'c',
        lbl_1299: 'c',
        lbl_1487: 'm',
        lbl_1595: 'k'
    }

    marker_map = dict([(k, 'o') for k in colour_map])
    marker_map[lbl_1299] = 's'
    marker_map[lbl_1487] = 'D'
    marker_map[lbl_1595] = '^'

    # plots: PCA of classifier vs RNA-Seq
    ttl = ("pca_%s-rnaseq_2d" % title) if SAVEFIG else None
    ad = {late_pass_lbl: y_cuff}
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plots: PCA of classifier vs Xiao-Nan
    ttl = ("pca_%s-1299_2d" % title) if SAVEFIG else None
    ad = {lbl_1299: y_1299}
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plots: PCA of classifier vs Xiao-Nan AND RNA-Seq
    ttl = ("pca_%s-1299-rnaseq_2d" % title) if SAVEFIG else None
    ad = {lbl_1299: y_1299, late_pass_lbl: y_cuff}
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 3D scatterplot of classifier data with ellipsoids and RNA-Seq sample overlaid
    ttl = ("pca_%s-rnaseq_3d" % title) if SAVEFIG else None
    ad = {late_pass_lbl: y_cuff}
    fig, ax_3d = plot_3d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 3D scatterplot of DIFFERENT microarray data with classifier ellipsoids (no classifier data)
    X2 = all_expr.loc[:, ncott_meta.index].transpose()
    m2 = ncott_meta.copy()
    idx2, labels2 = m2.subgroup.factorize()
    y2 = pca.transform(X2)


    ttl = ("pca_%s-ncott_3d" % title) if SAVEFIG else None
    ad = {}
    for i, l in enumerate(labels2):
        j = idx2 == i
        ad[l] = y2[j]

    fig, ax_3d = plot_3d(y, m.subgroup, colour_map, marker_map, title=ttl, plot_y=False, **ad)

    # plot: 3D scatterplot of classifier data with ellipsoids and Xiao-Nan data overlaid

    ttl = ("pca_%s-1299_3d" % title) if SAVEFIG else None
    ad = {lbl_1299: y_1299}
    fig, ax_3d = plot_3d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 3D scatterplot of classifier data with ellipsoids and both RNA-Seq AND Xiao-Nan data overlaid
    ttl = ("pca_%s-1299-rnaseq_3d" % title) if SAVEFIG else None
    ad = {late_pass_lbl: y_cuff, lbl_1299: y_1299}
    fig, ax_3d = plot_3d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 2D centroids + Xiao-Nan 1299, 1487, 1595
    ttl = ("pca_%s-1299-1487-1595_2d" % title) if SAVEFIG else None
    ad = {
        lbl_1299: y_1299,
        lbl_1487: y_1487,
        lbl_1595: y_1595
    }

    fig, axs = fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 3D centroids + Xiao-Nan 1299, 1487, 1595
    ttl = ("pca_%s-1299-1487-1595_3d" % title) if SAVEFIG else None
    fig, ax_3d = plot_3d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)