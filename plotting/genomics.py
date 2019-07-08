from matplotlib import pyplot as plt, patches, collections as plt_collections
from utils import genomics, setops, dictionary, log
from plotting import common


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


def plot_gene_transcripts(feature_res, bar_height=0.8, edge_buffer_bp=500, ax=None, patch_kwargs=None):
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

    for gene_ftr, d1 in feature_res.items():
        chroms.add(gene_ftr.chrom)
        if len(chroms) > 1:
            raise AttributeError("The supplied features inhabit different chromosomes (%s)." % ', '.join(chroms))
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

    # go back through and add directions
    # we'll assume all directions are the same across transcripts (FIXME!)
    chrom = None
    for i, (gene_ftr, d1) in enumerate(feature_res.items()):
        midy = (gene_y[i + 1] + gene_y[i]) / 2
        if gene_ftr.strand == '+':
            txt0 = "5'"
            txt1 = "3'"
        else:
            txt1 = "5'"
            txt0 = "3'"
        ax.text(xmin + edge_buffer_bp * 0.5, midy, txt0, horizontalalignment='right', va='center', fontsize=12)
        ax.text(xmax - edge_buffer_bp * 0.5, midy, txt1, horizontalalignment='left', va='center', fontsize=12)

    ax.set_xlabel('Chromosome %s' % the_chrom)
    ax.yaxis.set_visible(False)

    return ax


class MethylationExpressionLocusPlot(object):

    def __init__(self, row_names, fig_kws=None):
        if fig_kws is None:
            fig_kws = {'figsize': (8, 6)}
        self.fig_kws = fig_kws
        self.fig = None
        self.gs = None
        self.track_ax = None
        self.dm_axs = {}
        self.de_axs = {}
        self.coord_min = None
        self.coord_max = None
        self.include_de = False
        self.row_names = row_names
        self.nrows = len(row_names) + 1  # additional row for the track plot
        self.mdat_min = None
        self.mdat_max = None
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
        de_sharex = None

        # plot_gene_transcripts(feat_res, ax=track_ax)
        self.dm_axs = {}
        self.de_axs = {}
        for i, nm in enumerate(self.row_names):
            self.dm_axs[nm] = self.fig.add_subplot(self.gs[i, 0], sharex=dm_sharex)
            dm_sharex = self.dm_axs[nm]
            self.de_axs[nm] = self.fig.add_subplot(self.gs[i, 1], sharex=de_sharex, sharey=de_sharex)
            de_sharex = self.de_axs[nm]

    def check_setup(self):
        if self.fig is None:
            raise AttributeError("Must run setup_figure_axes() before plotting")

    def plot_tracks(self, feat_res):
        self.check_setup()
        plot_gene_transcripts(feat_res, ax=self.track_ax)

    def plot_dmrs(self, dmr_res, clusters):
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
                    # bar showing DMR coverage and delta
                    this_ax.barh(
                        0,
                        pc.coord_range[1] - pc.coord_range[0],
                        left=pc.coord_range[0],
                        height=d2['median_change'],
                        color=DIRECTION_COLOURS['up' if d2['median_change'] > 0 else 'down'],
                        align='edge',
                        zorder=15,
                        edgecolor='k',
                        linewidth=0.75
                    )

    def plot_de(self, de_res, fdr_cutoff=0.01):
        """

        :param de_res:
        :param fdr_cutoff:
        :return:
        """
        self.check_setup()

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
                    color=DIRECTION_COLOURS['up'] if this_res['logFC'] > 0 else DIRECTION_COLOURS['down']
                )

        self.set_de_on(True)

    def plot_m_values(self, mdat, probe_locations, comparisons, group_colours=None):
        """

        :param mdat: pd.DataFrame containing the data to plot. Columns are samples, rows are probes
        :param probe_locations: pd.Series containing the probe IDs to include and their genomic coordinates
        :param comparisons: Dictionary keyed by comparison (equivalent to row_names). Each entry is a dictionary keyed
        by group name (e.g. 'Disease' / 'Healthy') and with values giving the samples in that group. The sample names
        must be in the columns of `mdat`.
        :return:
        """
        if group_colours is None:
            all_groups = sorted(setops.reduce_union(*(t.keys() for t in comparisons.values())))
            n_groups = len(all_groups)
            group_colours = dict(zip(all_groups, common.get_best_cmap(n_groups)))

        # scatter plot individual probes
        ymin = 0
        ymax = 0
        for nm in self.row_names:
            grp_dict = comparisons[nm]
            this_ax = self.dm_axs[nm]
            for grp_nm, grp_samples in grp_dict.items():
                the_colour = group_colours.get(grp_nm, '0.5')
                for col, x in mdat.loc[probe_locations.index, grp_samples].iteritems():
                    this_ax.scatter(
                        probe_locations,
                        x.values,
                        c=the_colour,
                        marker='o',  ## TODO?
                        zorder=20,  ## TODO?
                        alpha=0.6,  ## TODO?
                        s=20,  ## TODO?
                        edgecolor='k',
                        linewidth=0.5
                    )
                    ymin = min(x.values.min(), ymin)
                    ymax = max(x.values.max(), ymax)
                    this_ax.set_ylabel(nm)
        self.mdat_min = ymin
        self.mdat_max = ymax

    def update_xlims(self, x0, x1):
        self.track_ax.set_xlim([x0, x1])
        self.dm_axs[self.row_names[-1]].set_xlim([x0, x1])

    def fix_axes(self):
        repr_dm_ax = self.dm_axs[self.row_names[-1]]
        repr_dm_ax.xaxis.set_ticklabels([])
        [ax.set_ylim([self.mdat_min * 1.05, self.mdat_max * 1.05]) for ax in self.dm_axs.values()]

        repr_de_ax = self.de_axs[self.row_names[-1]]
        plt.setp([t.yaxis for t in self.de_axs.values()], visible=False)
        repr_de_ax.set_ylim([-1, 1])
        plt.setp(repr_de_ax.xaxis.get_ticklabels(), rotation=90)
        repr_de_ax.set_xlabel('DE logFC')

        # declutter track axis by removing all borders except the bottom
        for k in ['left', 'right', 'top']:
            self.track_ax.spines[k].set_visible(False)
        self.track_ax.set_xlim(repr_dm_ax.get_xlim())

        self.gs.update(top=0.98, left=0.07, right=0.97, bottom=0.09, hspace=0.15, wspace=0.05)


class MethylationExpressionLocusPlotter(object):
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

    def plot_gene(self, gene, edge_buffer_bp=500, figsize=(8, 6), group_colours=None, fdr_cutoff=0.01):
        if self.dmr_res is None:
            raise AttributeError("Must run set_dmr_res() before plotting.")

        if self.mdat is None:
            raise AttributeError("Must run set_mvalues() befoire plotting.")

        if group_colours is None:
            all_groups = sorted(setops.reduce_union(*(t.keys() for t in self.dmr_comparison_groups.values())))
            n_groups = len(all_groups)
            group_colours = dict(zip(all_groups, common.get_best_cmap(n_groups)))

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
            for nm in self.dmr_comparison_groups:
                ix = self.de_res[nm]['Gene Symbol'] == gene
                if ix.sum() == 1:
                    this_de[nm] = self.de_res[nm].loc[ix].squeeze()
                elif ix.sum() > 1:
                    self.logger.warn("Gene %s. Found %d matching DE results in comparison %s. Why?", gene, int(ix.sum()),
                                     nm)
                    this_de[nm] = None
                else:
                    self.logger.warn("Gene %s. No matching DE result for comparison %s.", gene, nm)
                    this_de[nm] = None

        plot_obj = MethylationExpressionLocusPlot(
            self.dmr_comparison_groups.keys(),
            fig_kws={'figsize': figsize}
        )
        plot_obj.plot_tracks(feat_res)
        plot_obj.plot_dmrs(this_dmrs, self.dmr_res.clusters)
        plot_obj.plot_m_values(self.mdat, probe_locations, self.dmr_comparison_groups, group_colours=group_colours)
        plot_obj.plot_de(this_de, fdr_cutoff=fdr_cutoff)
        plot_obj.fix_axes()

        return plot_obj
