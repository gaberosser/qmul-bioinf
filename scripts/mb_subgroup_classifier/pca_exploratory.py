from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import os
from load_data import microarray_data, allen_human_brain_atlas, rnaseq_data
from microarray import process
from matplotlib import pyplot as plt
from matplotlib import animation

from plotting.pca import pca_plot_by_group_2d, pca_plot_by_group_3d, cluster_ellipsoid
from scripts.comparison_rnaseq_microarray import consts
from scripts.mb_subgroup_classifier.load import load_xz_rnaseq, load_xiaonan_microarray, load_sb_rnaseq


def animation_func(ax, n_frames):
    # fixed elevation, vary azimuth
    elev = ax.elev
    def animate(i):
        ax.view_init(elev, i * (360. / float(n_frames)))
    return animate


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


if __name__ == '__main__':

    ELL_P = 0.95
    N_PC = 3
    # ANIM = True
    ANIM = False
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

    # load XZ RNA-Seq count data
    # X_xz = load_xz_rnaseq(kind='htseq', yugene=True, gene_symbols=X.columns).transpose()
    X_xz = load_xz_rnaseq(kind='cuff', yugene=True, gene_symbols=X.columns).transpose()
    y_xz = pca.transform(X_xz)

    # load SB RNA-Seq count data
    # NB: have checked and using TPM rather than FPKM makes no difference, as expected
    X_sb = load_sb_rnaseq(yugene=True, gene_symbols=X.columns).transpose()
    y_sb = pca.transform(X_sb)

    # extract samples and reshape to make sure they are 2D arrays
    y_sb_1078 = y_sb[0].reshape((1, N_PC))
    y_sb_1595 = y_sb[1].reshape((1, N_PC))
    y_sb_1487 = y_sb[2].reshape((1, N_PC))

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
    lbl_1078_this = 'ICb1078 (this study)'
    lbl_1299_this = 'ICb1299 (this study)'
    lbl_1487_this = 'ICb1487 (this study)'
    lbl_1595_this = 'ICb1595 (this study)'
    lbl_1078_zhao = 'ICb1078 (Zhao et al.)'
    lbl_1299_zhao = 'ICb1299 (Zhao et al.)'
    lbl_1487_zhao = 'ICb1487 (Zhao et al.)'
    lbl_1595_zhao = 'ICb1595 (Zhao et al.)'

    colour_map = {
        'Group 3': '#F2EA00',
        'Group 4': '#2A8C43',
        'WNT': '#2D438E',
        'SHH': '#E5161C',
        'control': 'gray',
        lbl_1299_zhao: 'c',
        lbl_1487_zhao: 'm',
        lbl_1595_zhao: 'k',
        lbl_1078_zhao: '#ff6600',
    }
    colour_map[lbl_1299_this] = colour_map[lbl_1299_zhao]
    colour_map[lbl_1078_this] = colour_map[lbl_1078_zhao]
    colour_map[lbl_1487_this] = colour_map[lbl_1487_zhao]
    colour_map[lbl_1595_this] = colour_map[lbl_1595_zhao]

    marker_map = dict([(k, 'o') for k in colour_map])
    marker_map[lbl_1299_zhao] = 'D'
    marker_map[lbl_1487_zhao] = 'D'
    marker_map[lbl_1595_zhao] = 'D'
    marker_map[lbl_1078_this] = 's'
    marker_map[lbl_1299_this] = 's'
    marker_map[lbl_1487_this] = 's'
    marker_map[lbl_1595_this] = 's'


    # plots: PCA of classifier vs RNA-Seq
    ttl = ("pca_%s-rnaseq_2d" % title) if SAVEFIG else None
    ad = {lbl_1299_this: y_xz}
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plots: PCA of classifier vs Xiao-Nan
    ttl = ("pca_%s-1299_2d" % title) if SAVEFIG else None
    ad = {lbl_1299_zhao: y_1299}
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plots: PCA of classifier vs Xiao-Nan AND RNA-Seq
    ttl = ("pca_%s-1299-rnaseq_2d" % title) if SAVEFIG else None
    ad = {lbl_1299_zhao: y_1299, lbl_1299_this: y_xz}
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 3D scatterplot of classifier data with ellipsoids and RNA-Seq sample overlaid
    ttl = ("pca_%s-rnaseq_3d" % title) if SAVEFIG else None
    ad = {lbl_1299_this: y_xz}
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
    ad = {lbl_1299_zhao: y_1299}
    fig, ax_3d = plot_3d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 3D scatterplot of classifier data with ellipsoids and both RNA-Seq AND Xiao-Nan data overlaid
    ttl = ("pca_%s-1299-rnaseq_3d" % title) if SAVEFIG else None
    ad = {lbl_1299_this: y_xz, lbl_1299_zhao: y_1299}
    fig, ax_3d = plot_3d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 2D centroids + Xiao-Nan 1299, 1487, 1595
    ttl = ("pca_%s-1299-1487-1595_2d" % title) if SAVEFIG else None
    ad = {
        lbl_1299_zhao: y_1299,
        lbl_1487_zhao: y_1487,
        lbl_1595_zhao: y_1595
    }

    fig, axs = fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 3D centroids + Xiao-Nan 1299, 1487, 1595
    ttl = ("pca_%s-1299-1487-1595_3d" % title) if SAVEFIG else None
    fig, ax_3d = plot_3d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 3D ellipsoids + paired RNA-Seq and Xiao-Nan data overlaid
    ttl = ("pca_%s-early_late_all_3d" % title) if SAVEFIG else None
    ad = {
        lbl_1299_this: y_xz,
        lbl_1299_zhao: y_1299,
        lbl_1487_this: y_sb_1487,
        lbl_1487_zhao: y_1487,
        lbl_1595_this: y_sb_1595,
        lbl_1595_zhao: y_1595,
        lbl_1078_this: y_sb_1078
    }
    fig, ax_3d = plot_3d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    ttl = ("pca_%s-early_late_all_2d" % title) if SAVEFIG else None
    fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    # plot: 2D centroids + 1595 old and new
    ttl = ("pca_%s-1595_2d" % title) if SAVEFIG else None
    ad = {
        lbl_1595_zhao: y_1595,
        lbl_1595_this: y_sb_1595,
    }

    fig, axs = fig, axs = plot_2d(y, m.subgroup, colour_map, marker_map, title=ttl, **ad)

    ###### experimental: plotting with Plotly ######
    if False:

        import plotly.plotly as py
        from plotly import graph_objs as go

        scatter = []
        wireframe = []

        # input args
        pca_data = y
        components = (0, 1, 2)
        ellipse_p = 0.99
        npt = 80
        markersize = 5

        # formatting
        plotly_symbol_map = {
            'o': 'circle',
            'D': 'diamond',
            's': 'square'
        }

        idx, labels = m.subgroup.factorize()

        for i, l in enumerate(labels):
            if l not in colour_map:
                continue
            c = colour_map[l]
            mk = plotly_symbol_map[marker_map[l]]
            j = idx == i
            marker = {
                'size': markersize,
                'color': c,
                # 'line': {'color': c, 'width': 0.5},
                'opacity': 0.8,
                'symbol': mk
            }
            trace = go.Scatter3d(
                x=pca_data[j, components[0]],
                y=pca_data[j, components[1]],
                z=pca_data[j, components[2]],
                mode='markers',
                marker=marker,
                name=l,
                showlegend=True,
            )
            scatter.append(trace)

        for k, v in ad.items():
            c = colour_map[k]
            mk = plotly_symbol_map[marker_map[k]]
            marker = {
                'size': markersize,
                'color': c,
                # 'line': {'color': c, 'width': 0.5},
                'opacity': 0.8,
                'symbol': mk
            }
            trace = go.Scatter3d(
                x=v[:, components[0]],
                y=v[:, components[1]],
                z=v[:, components[2]],
                mode='markers',
                marker=marker,
                name=k,
                showlegend=True,
            )
            scatter.append(trace)

        # plot_ellipsoids
        for i, l in enumerate(labels):
            if l not in colour_map:
                continue
            c = colour_map[l]
            j = idx == i
            loc, radii, rot = cluster_ellipsoid(pca_data[j][:, components], p=ellipse_p)
            u = np.linspace(0.0, 2.0 * np.pi, npt)
            v = np.linspace(0.0, np.pi, npt)

            # cartesian coordinates that correspond to the spherical angles:
            xx = radii[0] * np.outer(np.cos(u), np.sin(v))
            yy = radii[1] * np.outer(np.sin(u), np.sin(v))
            zz = radii[2] * np.outer(np.ones_like(u), np.cos(v))

            # rotate accordingly
            for i in range(len(xx)):
                for j in range(len(xx)):
                    [xx[i, j], yy[i, j], zz[i, j]] = np.dot([xx[i, j], yy[i, j], zz[i, j]], rot) + loc

            line = {
                'color': c,
                'width': 1.0,
            }
            # Mesh3d looks nice, but messes up colours inside closed objects
            # trace = go.Mesh3d(
            #     x=xx.flatten(),
            #     y=yy.flatten(),
            #     z=zz.flatten(),
            #     alphahull=0,
            #     opacity=0.4,
            #     color=c,
            #     name=l,
            #     showlegend=True
            # )
            # wireframe.append(trace)

            # build the ellipsoids out of lots of 3D lines
            leg = True
            for x1, y1, z1 in zip(xx, yy, zz):
                trace = go.Scatter3d(
                    x=x1, y=y1, z=z1, mode='lines', line=line, showlegend=leg, name=l
                )
                # only show legend on the first
                leg=False
                wireframe.append(trace)

        # TODO: disable hover projection (can do this manually until then)
        layout = go.Layout(
            title='Wireframe Plot',
            scene=dict(
                xaxis=dict(
                    gridcolor='rgb(255, 255, 255)',
                    zerolinecolor='rgb(255, 255, 255)',
                    showbackground=True,
                    backgroundcolor='rgb(230, 230,230)'
                ),
                yaxis=dict(
                    gridcolor='rgb(255, 255, 255)',
                    zerolinecolor='rgb(255, 255, 255)',
                    showbackground=True,
                    backgroundcolor='rgb(230, 230,230)'
                ),
                zaxis=dict(
                    gridcolor='rgb(255, 255, 255)',
                    zerolinecolor='rgb(255, 255, 255)',
                    showbackground=True,
                    backgroundcolor='rgb(230, 230,230)'
                )
            ),
            showlegend=True,
        )

        fig = go.Figure(data=scatter + wireframe, layout=layout)
        # publish to the web
        py.plot(fig, filename='wire_and_scatter')
