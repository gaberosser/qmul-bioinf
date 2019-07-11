import numpy as np
from matplotlib import pyplot as plt, patches, collections as plt_collections, colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from utils import genomics, setops, dictionary, log
from plotting import common
logger = log.get_console_logger()

PATCH_KWS = {
    'edgecolor': 'k',
    'linewidth': 0.75,
    'zorder': 15,
    'facecolor': 'w'
}


FEATURE_FORMATTER = dict([
    (k, dict(PATCH_KWS)) for k in ['default', 'exon', 'five_prime_utr', 'three_prime_utr']
])
FEATURE_FORMATTER['exon']['facecolor'] = 'k'
FEATURE_FORMATTER['three_prime_utr']['zorder'] = 16
FEATURE_FORMATTER['five_prime_utr']['zorder'] = 16


DIRECTION_COLOURS = {
    'down': '#89CD61',
    'up': '#FF381F',
}


def direction_colour_getter(direction_colours=DIRECTION_COLOURS, vmin=None, vmax=None):
    """
    Generate a direction colour getters, which returns a hexadecimal colour code when called with a direction.
    :param direction_colours: Various inputs supported
    string: The name of a matplotlib colourmap. Must specify vmin and vmax.
    dict: Must be keyed by 'up' and 'down', this is used to create a binary colourmap. No need to specify vmin or vmax.
    array of hex strings: Used to generate a discrete colourmap between the specified vmin and vmax values.
    :param vmin, vmax: The min and max values used for the scalar colour mapping.
    :return: Function that takes a single scalar input and returns a hex colour.
    """
    if isinstance(direction_colours, str):
        if vmin is None or vmax is None:
            raise AttributeError("When supplying the name of a matplotlib colourmap, must also supply vmin and vmax.")
        return common.continuous_cmap(vmin, vmax, cmap=direction_colours)
    elif isinstance(direction_colours, dict):
        if len({'up', 'down'}.intersection(direction_colours.keys())) != 2:
            raise ValueError("When supplying a dictionary, it must contain the keys 'up' and 'down'.")
        return lambda x: direction_colours['up' if x > 0 else 'down']
    elif hasattr(direction_colours, '__iter__'):
        if vmin is None or vmax is None:
            raise AttributeError("When supplying an iterable of colours, must also supply vmin and vmax.")
        norm = common.Normalize(vmin=vmin, vmax=vmax)
        lsm = colors.LinearSegmentedColormap.from_list(direction_colours)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=lsm)
        return lambda x: colors.to_hex(sm.to_rgba(x))
    elif callable(direction_colours):
        return direction_colours
    else:
        ## TODO: support ColorMap objects from matplotlib??
        raise NotImplementedError("Unsupported direction_colours object.")


DIRECTION_COLOUR_GETTER = direction_colour_getter(DIRECTION_COLOURS)


def plot_gene_transcripts(feature_res, bar_height=0.8, edge_buffer_bp=500, ax=None, patch_kwargs=None, add_arrow=True):
    """

    :param feature_res: Output from utils.genomics.get_features_from_gtf. Nested dictionary. Keys are features. The
    outer level corresponds to genes. The inner level corresponds to transcripts. The bottom level is a list of
    features to plot.
    :return:
    """
    if ax is None:
        ax = plt.gca()

    if patch_kwargs is None:
        patch_kwargs = {}

    # plot the tracks
    xmin = 1e12
    xmax = -1

    y = 0
    gene_y = []
    chroms = set()
    the_strand = None

    for gene_ftr, d1 in feature_res.items():
        chroms.add(gene_ftr.chrom)
        if len(chroms) > 1:
            raise AttributeError("The supplied features inhabit different chromosomes (%s)." % ', '.join(chroms))
        if the_strand is None:
            the_strand = gene_ftr.strand
        elif the_strand != gene_ftr.strand:
            logger.warn("Supplied features are present on both strands. Axis labels will be confusing.")

        gene_y.append(y)
        if y != 0:
            # start new gene section
            ax.axhline(y, c='k', linewidth=2.)
            y += 1
        for transcript_ftr, arr in d1.items():
            if len(arr) == 0:
                continue
            arr = sorted(arr, key=lambda x: x.start)
            this_patches = []
            for ftr in arr:
                the_patch_kws = dict(FEATURE_FORMATTER.get(ftr.featuretype, FEATURE_FORMATTER['default']))
                the_patch_kws.update(patch_kwargs)
                the_patch = patches.Rectangle(
                    (ftr.start - 0.5, y - bar_height * 0.5),
                    ftr.stop - ftr.start + 0.5,
                    bar_height,
                    **the_patch_kws
                )
                this_patches.append(the_patch)
                xmin = min(xmin, ftr.start - edge_buffer_bp)
                xmax = max(xmax, ftr.stop + edge_buffer_bp)
            x0 = transcript_ftr.start
            x1 = transcript_ftr.stop
            ax.add_collection(plt_collections.PatchCollection(this_patches, match_original=True, zorder=20))
            ax.plot([x0, x1], [y, y], color='k', linestyle='-', linewidth=2., zorder=15)
            y += 1
    gene_y.append(y)
    # to reach this point without exception, we must have a single chromosome
    the_chrom = list(chroms)[0]

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([-bar_height, y - 1 + bar_height])

    if add_arrow:
        if the_strand == '+':
            loc_str = 'bottom right'
        else:
            loc_str = 'bottom left'
        common.arrowed_spines(ax, lw=0.3, locations=(loc_str,))

    ax.set_xlabel('Chromosome %s' % the_chrom)
    ax.yaxis.set_visible(False)

    return ax


class MexLocusPlot(object):

    def __init__(self, row_names, fig_kws=None, edge_buffer=500):
        if fig_kws is None:
            fig_kws = {'figsize': (8, 6)}
        self.fig_kws = fig_kws
        self.fig = None
        self.gs = None
        self.track_ax = None
        self.m_axs = {}
        self.de_axs = {}
        self.coord_min = None
        self.coord_max = None
        self.track_arrow = None

        self.row_names = row_names
        self.nrows = len(row_names) + 1  # additional row for the track plot

        self.mdat_min = None
        self.mdat_max = None
        self.de_min = None
        self.de_max = None
        self.coord_min = None
        self.coord_max = None
        self.strand = None

        self.edge_buffer = edge_buffer

        self.setup_figure_axes()

    def set_de_on(self, toggle=True):
        if toggle:
            self.gs.set_width_ratios((10, 1))
        else:
            self.gs.set_width_ratios((10000, 1))

    def setup_figure_axes(self):
        self.fig = plt.figure(**self.fig_kws)
        self.gs = plt.GridSpec(nrows=self.nrows, ncols=2, width_ratios=(10000, 1))
        self.track_ax = self.fig.add_subplot(self.gs[-1, 0])
        dm_sharex = None
        de_sharey = None

        self.m_axs = {}
        self.de_axs = {}
        for i, nm in enumerate(self.row_names):
            self.m_axs[nm] = self.fig.add_subplot(self.gs[i, 0], sharex=dm_sharex)
            dm_sharex = self.m_axs[nm]
            self.de_axs[nm] = self.fig.add_subplot(self.gs[i, 1], sharey=de_sharey)
            de_sharey = self.de_axs[nm]

    def check_setup(self):
        if self.fig is None:
            raise AttributeError("Must run setup_figure_axes() before plotting")

    def plot_tracks(self, feat_res):
        self.check_setup()
        plot_gene_transcripts(feat_res, ax=self.track_ax, add_arrow=False)

        if self.coord_max is None:
            self.coord_min, self.coord_max = self.track_ax.get_xlim()

        # get strand - assume this is the same for all features
        self.strand = feat_res.keys()[0].strand

    def plot_dmrs(self, dmr_res, clusters, colours=DIRECTION_COLOUR_GETTER):
        """
        :param dmr_res: Nested dictionary. First level is keyed by cluster ID. Second level is keyed by comparison.
        :param clusters: Dictionary keyed by cluster ID. Each entry is an instance of ProbeCluster
        :return:
        """
        self.check_setup()

        for cid, d1 in dmr_res.items():
            pc = clusters[cid]
            for nm, d2 in d1.items():
                if nm in self.row_names:
                    this_ax = self.m_axs[nm]
                    # bar showing DMR coverage and delta
                    this_ax.barh(
                        0,
                        pc.coord_range[1] - pc.coord_range[0],
                        left=pc.coord_range[0],
                        height=d2['median_change'],
                        color=colours(d2['median_change']),
                        align='edge',
                        zorder=15,
                        edgecolor='k',
                        linewidth=0.75
                    )

    def plot_de(self, de_res, fdr_cutoff=0.01, colours=DIRECTION_COLOUR_GETTER):
        """
        :param de_res:
        :param fdr_cutoff:
        :return:
        """
        self.check_setup()
        de_min = 0
        de_max = 0

        for nm in self.row_names:
            this_de_ax = self.de_axs[nm]
            this_res = de_res[nm]
            if this_res is not None:
                # bar chart to represent DE
                this_de_ax.barh(
                    [0],
                    [this_res['logFC']],
                    edgecolor='k' if this_res['FDR'] < fdr_cutoff else 'none',
                    linewidth=1.,
                    color=colours(this_res['logFC'])
                )
                de_min = min(de_min, this_res['logFC'])
                de_max = max(de_max, this_res['logFC'])

        self.set_de_on(True)
        self.de_min = de_min
        self.de_max = de_max

    def plot_m_values(
            self,
            mdat,
            probe_locations,
            comparisons,
            colours='default',
            markers='default',
            zorder='default',
            alpha='default',
            size='default',
    ):
        """

        :param mdat: pd.DataFrame containing the data to plot. Columns are samples, rows are probes
        :param probe_locations: pd.Series containing the probe IDs to include and their genomic coordinates
        :param comparisons: Dictionary keyed by comparison (equivalent to row_names). Each entry is a dictionary keyed
        by group name (e.g. 'Disease' / 'Healthy') and with values giving the samples in that group. The sample names
        must be in the columns of `mdat`.
        :param colours: Dictionary keyed by group name (e.g. 'Disease') giving the colour to use for that group.
        Defaults are used if not supplied. To disable colours, set to None.
        :param markers: Dictionary keyed by group name giving the marker to use for that group.
        Defaults are used if not supplied. To use circle markers for everything, set to None.
        :param zorder: Dictionary keyed by group name giving the zorder to use for that group.
        Defaults are used if not supplied. To use matplotlib defaults for everything, set to None.
        :param alpha: Dictionary keyed by group name giving the alpha to use for that group.
        Defaults are used if not supplied. To use matplotlib defaults for everything, set to None.
        :return:
        """
        all_groups = sorted(setops.reduce_union(*(t.keys() for t in comparisons.values())))
        n_groups = len(all_groups)

        def set_property(x, default, default_static):
            if x == 'default':
                out = dict(zip(all_groups, default))
            elif x is None:
                out = dict([(k, default_static) for k in all_groups])
            elif not hasattr(x, 'get'):
                # single value supplied
                out = dict([(k, x) for k in all_groups])
            else:
                out = x
            return out

        colours = set_property(
            colours,
            common.get_best_cmap(n_groups),
            '0.5'
        )
        markers = set_property(
            markers,
            common.get_best_marker_map(n_groups),
            'o'
        )
        zorder = set_property(
            zorder,
            range(20, 20 + n_groups),
            20
        )
        # default alpha will be based on zorder
        a = sorted([(k, zorder[k]) for k in all_groups], key=lambda x: x[1])
        a_ix = dict([(t[0], i) for i, t in enumerate(a)])
        alpha_values = np.linspace(0.4, 0.6, n_groups)
        alpha_default = [alpha_values[a_ix[k]] for k in all_groups]

        alpha = set_property(
            alpha,
            alpha_default,
            '0.6'
        )

        # default size will be based on zorder
        s_values = range(20, 20 + n_groups)
        s_default = [s_values[a_ix[k]] for k in all_groups]
        size = set_property(
            size,
            s_default,
            20
        )

        # scatter plot individual probes
        ymin = 0
        ymax = 0
        for nm in self.row_names:
            grp_dict = comparisons[nm]
            this_ax = self.m_axs[nm]
            for grp_nm, grp_samples in grp_dict.items():
                the_colour = colours.get(grp_nm)
                the_marker = markers.get(grp_nm)
                the_z = zorder.get(grp_nm)
                the_alpha = alpha.get(grp_nm)
                the_s = size.get(grp_nm)
                for col, x in mdat.loc[probe_locations.index, grp_samples].iteritems():
                    this_ax.scatter(
                        probe_locations,
                        x.values,
                        c=the_colour,
                        marker=the_marker,
                        zorder=the_z,
                        alpha=the_alpha,
                        s=the_s,
                        edgecolor='k',
                        linewidth=0.5
                    )
                    ymin = min(x.values.min(), ymin)
                    ymax = max(x.values.max(), ymax)
                    this_ax.set_ylabel(nm)
        self.mdat_min = ymin
        self.mdat_max = ymax

        if self.coord_max is None:
            self.coord_min = probe_locations.min()
            self.coord_max = probe_locations.max()
        else:
            self.coord_min = min(probe_locations.min(), self.coord_min)
            self.coord_max = max(probe_locations.max(), self.coord_max)

    def update_xlims(self, x0, x1):
        self.track_ax.set_xlim([x0, x1])
        self.m_axs[self.row_names[-1]].set_xlim([x0, x1])
        # add arrow to track axis
        self.add_arrow_to_track_ax()

    def apply_edge_buffer(self, edge_buffer=None):
        if self.coord_max is None:
            raise ValueError("Must plot tracks or M values before applying the edge buffer.")
        if edge_buffer is None:
            edge_buffer = self.edge_buffer
        else:
            self.edge_buffer = edge_buffer
        self.update_xlims(self.coord_min - edge_buffer, self.coord_max + edge_buffer)

    def add_arrow_to_track_ax(self):
        # if arrow present, remove it first
        if self.track_arrow is not None:
            for v in self.track_arrow.values():
                v.remove()
            self.track_arrow = None

        # add arrow to track axis
        if self.strand == '+':
            loc_str = 'bottom right'
        else:
            loc_str = 'bottom left'
        # TODO: auto-select hardcoded arrow width parameters?
        self.track_arrow = common.arrowed_spines(
            self.track_ax,
            lw=0.3,
            locations=(loc_str,),
            x_length_fraction=0.025,
            y_width_fraction=0.15
        )

    def fix_axes(self):
        self.apply_edge_buffer()

        repr_dm_ax = self.m_axs[self.row_names[-1]]
        repr_dm_ax.xaxis.set_ticklabels([])
        [ax.set_ylim([self.mdat_min * 1.05, self.mdat_max * 1.05]) for ax in self.m_axs.values()]

        repr_de_ax = self.de_axs[self.row_names[-1]]
        plt.setp([t.yaxis for t in self.de_axs.values()], visible=False)
        repr_de_ax.set_ylim([-1, 1])
        plt.setp(repr_de_ax.xaxis.get_ticklabels(), rotation=90)
        repr_de_ax.set_xlabel('DE logFC')
        [ax.set_xlim([self.de_min * 1.05, self.de_max * 1.05]) for ax in self.de_axs.values()]
        [self.de_axs[k].set_xticklabels([]) for k in self.row_names[:-1]]

        # declutter track axis by removing all borders except the bottom
        for k in ['left', 'right', 'top']:
            self.track_ax.spines[k].set_visible(False)

        self.gs.update(top=0.98, left=0.07, right=0.97, bottom=0.09, hspace=0.15, wspace=0.05)

    def set_title(self, ttl, **kwargs):
        self.m_axs[self.row_names[0]].set_title(ttl, **kwargs)
        self.gs.update(top=0.95)

    def remove_title(self):
        self.m_axs[self.row_names[0]].set_title("")
        self.gs.update(top=0.98)


class MexLocusPlotFloatingDMR(MexLocusPlot):
    """
    Mex locus plot where DMRs are shown on a axis inset into the M value axis (floating at the top) to improve
    clarity.
    """
    def __init__(self, *args, **kwargs):
        self.dm_axs = {}
        super(MexLocusPlotFloatingDMR, self).__init__(*args, **kwargs)

    def update_xlims(self, x0, x1):
        super(MexLocusPlotFloatingDMR, self).update_xlims(x0, x1)
        plt.setp(self.dm_axs.values(), xlim=[x0, x1])

    def fix_axes(self):
        super(MexLocusPlotFloatingDMR, self).fix_axes()
        for ax in self.dm_axs.values():
            ax.patch.set_visible(False)
            plt.setp(ax.spines.values(), visible=False)
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            ax.set_ylim([-.55, .55])

    def setup_figure_axes(self):
        super(MexLocusPlotFloatingDMR, self).setup_figure_axes()
        # add DMR insets
        for nm in self.row_names:
            parent_ax = self.m_axs[nm]
            self.dm_axs[nm] = inset_axes(
                parent_ax,
                width="100%",
                height="15%",
                loc=3,  # lower left
                bbox_to_anchor=(0., .85, 1., 1.),
                bbox_transform=parent_ax.transAxes,
                borderpad=0.
            )

    def plot_dmrs(self, dmr_res, clusters, colours=DIRECTION_COLOUR_GETTER):
        """
        :param dmr_res: Nested dictionary. First level is keyed by cluster ID. Second level is keyed by comparison.
        :param clusters: Dictionary keyed by cluster ID. Each entry is an instance of ProbeCluster
        :return:
        """
        self.check_setup()

        for cid, d1 in dmr_res.items():
            pc = clusters[cid]
            for nm, d2 in d1.items():
                if nm in self.row_names:
                    this_ax = self.dm_axs[nm]
                    the_patch = patches.Rectangle(
                        (pc.coord_range[0], -.5),
                        width=pc.coord_range[1] - pc.coord_range[0],
                        height=1.,
                        facecolor=colours(d2['median_change']),
                        edgecolor='k',
                        linewidth=0.75
                    )
                    this_ax.add_patch(the_patch)


class MexLocusPlotter(object):
    def __init__(self, tax_id=9606, genome_version='GRCh37', gtf_fn=None):
        self.tax_id = tax_id
        self.genome_version = genome_version
        if gtf_fn is None:
            gtf_fn = genomics.get_reference_genome_gtf(tax_id, version=genome_version)
        self.gtf_fn = gtf_fn
        self.db = genomics.GtfAnnotation(gtf_fn)
        self.mdat = None
        self.dmr_res = None
        self.anno = None
        self.de_res = None
        self.dmr_comparison_groups = None
        self.logger = log.get_console_logger(self.__class__.__name__)

        # default plotting parameters
        self.fig_kws = {}
        self.m_plot_kws = {}

        self.de_direction_colour = None
        self.dm_direction_colour = None

        self.set_plot_parameters()


    def check_data_compat(self):

        if self.dmr_comparison_groups is not None:

            if self.dmr_res is not None:
                for grp_name, grp_dict in self.dmr_comparison_groups.items():
                    if grp_name not in self.dmr_res:
                        raise ValueError("Group %s is not in the DMR results" % grp_name)

            if self.de_res is not None:
                for grp_name, grp_dict in self.dmr_comparison_groups.items():
                    if grp_name not in self.de_res:
                        raise ValueError("Group %s is not in the DE results" % grp_name)

            if self.mdat is not None:
                for grp_name, grp_dict in self.dmr_comparison_groups.items():
                    all_samples = list(setops.reduce_union(*grp_dict.values()))
                    if len(self.mdat.columns.intersection(all_samples)) != len(all_samples):
                        raise ValueError("Group %s contains samples that are missing from mdat" % grp_name)

    def set_mvalues(self, dat):
        self.mdat = dat
        self.check_data_compat()

    def set_dmr_res(self, dmr_res, comparison_groups):
        self.dmr_res = dmr_res
        self.dmr_comparison_groups = comparison_groups
        self.anno = dmr_res.anno
        self.check_data_compat()

    def set_de_res(self, de_res):
        self.de_res = de_res
        self.check_data_compat()

    def set_plot_parameters(
            self,
            figsize=(8, 6),
            colours='default',
            markers='default',
            zorder='default',
            alpha='default',
            size='default',
            de_direction_colours=DIRECTION_COLOURS,
            dm_direction_colours=DIRECTION_COLOURS,
            de_vmin=None,
            de_vmax=None,
            dm_vmin=None,
            dm_vmax=None
    ):
        self.m_plot_kws = {
            'colours': colours,
            'markers': markers,
            'zorder': zorder,
            'alpha': alpha,
            'size': size
        }
        self.fig_kws = {
            'figsize': figsize
        }
        self.de_direction_colour = direction_colour_getter(de_direction_colours, vmin=de_vmin, vmax=de_vmax)
        self.dm_direction_colour = direction_colour_getter(dm_direction_colours, vmin=dm_vmin, vmax=dm_vmax)

    def plot_legend(self):
        """
        Generate a figure showing the interpretation of the various colours / markers
        :return:
        """
        pass

    def plot_gene(
        self,
        gene,
        fdr_cutoff=0.01,
        edge_buffer_bp=500,
        dmr_display='floating'
    ):
        class_dict = {
            'floating': MexLocusPlotFloatingDMR,
            'overlay': MexLocusPlot
        }
        if self.dmr_res is None:
            raise AttributeError("Must run set_dmr_res() before plotting.")

        if self.mdat is None:
            raise AttributeError("Must run set_mvalues() befoire plotting.")

        if dmr_display.lower() not in class_dict:
            raise AttributeError("Unsupported dmr_display %s. Options are %s." % (
                dmr_display,
                ', '.join(class_dict.keys())
            ))
        the_class = class_dict[dmr_display.lower()]

        feat_res = self.db.hierarchical_feature_search(
            gene,
            [['transcript'], ['exon', 'three_prime_utr', 'five_prime_utr']]
        )
        if len(feat_res) == 0:
            raise ValueError("No search results for gene %s" % gene)

        # identify clusters (not necessarily DMRs) associated with this gene
        this_clusters = dict([
                                 (k, v) for k, v in self.dmr_res.clusters.items() if len([x for x in v.genes if x[0] == gene])
                                 ])

        # get coordinate range (min, max) - this only requires the parental gene search results
        xmin = 1e12
        xmax = 0
        this_chrom = None

        for ftr in feat_res:
            if this_chrom is None:
                this_chrom = ftr.chrom
            elif this_chrom != ftr.chrom:
                raise AttributeError("Multiple gene search results on DIFFERENT chromosomes are not supported")
            check_intersection = lambda x: \
                (ftr.start <= x.coord_range[0] <= ftr.stop) or \
                (ftr.start <= x.coord_range[1] <= ftr.stop) or \
                (x.coord_range[0] < ftr.start and x.coord_range[1] > ftr.stop)
            this_clusters.update(
                dict([
                         (k, v) for k, v in self.dmr_res.clusters.items() if v.chr == ftr.chrom and check_intersection(v)
                         ])
            )
            xmin = min(xmin, ftr.start)
            xmax = max(xmax, ftr.stop)

        # include specified margins
        xmin -= edge_buffer_bp
        xmax += edge_buffer_bp

        # extract DMRs for plotting

        this_dmrs = {}
        for k, v in this_clusters.items():
            for k2, t in self.dmr_res.results.items():
                if t[k].get('padj', 1.) < fdr_cutoff:
                    this_dmrs.setdefault(k, {})[k2] = t[k]

        # add all probes involved in DMRs (they may be outside of the region)
        include_probes = set()

        for cid, d1 in this_dmrs.items():
            pc = this_clusters[cid]
            include_probes.update(pc.pids)

        # add all probes in the region
        ix = (self.anno.CHR == this_chrom) & (self.anno.MAPINFO >= xmin) & (self.anno.MAPINFO <= xmax)
        include_probes.update(self.anno.index[ix])

        # reduce to those probes in the data
        include_probes = self.mdat.index.intersection(include_probes)

        # sort by coordinate
        probe_locations = self.anno.loc[include_probes, 'MAPINFO'].sort_values()

        # prepare DE results for plotting (if present)
        this_de = {}
        if self.de_res is not None:
            de_res_present = []
            for nm in self.dmr_comparison_groups:
                ix = self.de_res[nm]['Gene Symbol'] == gene
                if ix.sum() == 1:
                    this_de[nm] = self.de_res[nm].loc[ix].squeeze()
                    de_res_present.append(nm)
                elif ix.sum() > 1:
                    self.logger.warn("Gene %s. Found %d matching DE results in comparison %s. Why?", gene, int(ix.sum()),
                                     nm)
                    this_de[nm] = None
                else:
                    self.logger.warn("Gene %s. No matching DE result for comparison %s.", gene, nm)
                    this_de[nm] = None
            if len(de_res_present) > 0:
                self.logger.info(
                    "Gene %s. Found DE results in %d comparisons: %s.",
                    gene,
                    len(de_res_present),
                    ', '.join(de_res_present)
                )
            else:
                self.logger.warn("Gene %s. No DE results found.", gene)


        plot_obj = the_class(
            self.dmr_comparison_groups.keys(),
            fig_kws=self.fig_kws
        )
        plot_obj.plot_tracks(feat_res)
        plot_obj.plot_dmrs(this_dmrs, self.dmr_res.clusters, colours=self.dm_direction_colour)

        plot_obj.plot_m_values(self.mdat, probe_locations, self.dmr_comparison_groups, **self.m_plot_kws)
        plot_obj.plot_de(this_de, fdr_cutoff=fdr_cutoff, colours=self.de_direction_colour)
        plot_obj.fix_axes()

        return plot_obj
