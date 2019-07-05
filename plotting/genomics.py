from matplotlib import pyplot as plt, patches, collections as plt_collections
from utils import genomics, setops

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
    def __init__(self, tax_id=9606, genome_version='GRCh37', gtf_fn=None):
        self.tax_id = tax_id
        self.genome_version = genome_version
        if gtf_fn is None:
            gtf_fn = genomics.get_reference_genome_fasta(tax_id, version=genome_version)
        self.gtf_fn = gtf_fn
        self.meta = None
        self.mdat = None
        self.dmr_res = None
        self.anno = None
        self.de_res = None
        self.dmr_comparison_groups = None

    def check_data_compat(self):
        if self.meta is not None and self.mdat is not None:
            if self.mdat.columns.sort_values() != self.meta.index.sort_values():
                raise ValueError("meta index and mdat columns must match")

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

    def set_meta(self, meta):
        self.meta = meta
        self.check_data_compat()

    def set_dmr_res(self, dmr_res, comparison_groups):
        self.dmr_res = dmr_res
        self.dmr_comparison_groups = comparison_groups
        self.anno = dmr_res.anno
        self.check_data_compat()

    def set_de_res(self, de_res):
        self.de_res = de_res
        self.check_data_compat()

    def plot_gene(self, gene, rois, figsize=(8, 6)):
        this_clusters = dict([
            (k, v) for k, v in self.dmr_res.clusters.items() if len([x for x in v.genes if x[0] == gene])
        ])
        feat_res = genomics.get_features_from_gtf(gene, tax_id=self.tax_id, version=self.genome_version)
        for i, (ftr, x) in enumerate(feat_res.items()):
            # get clusters and DMRs
            check_intersection = lambda x: \
                (ftr.start <= x.coord_range[0] <= ftr.stop) or \
                (ftr.start <= x.coord_range[1] <= ftr.stop) or \
                (x.coord_range[0] < ftr.start and x.coord_range[1] > ftr.stop)
            this_clusters.update(
                dict([
                    (k, v) for k, v in self.dmr_res.clusters.items() if v.chr == ftr.chrom and check_intersection(v)
                ])
            )
        this_dmrs = {}
        for k, v in this_clusters.items():
            for k2, t in self.dmr_res.results.items():
                if t.get('rej_h0', False):
                    this_dmrs.setdefault(k, {})[k2] = t

        nrows = len(self.dmr_comparison_groups) + 1  # one ax for each comparison plus one for the track

        # first column shows methylation data, second RNA-Seq DE (if supplied)
        fig = plt.figure(figsize=figsize)
        if self.de_res is None:
            gs = plt.GridSpec(nrows=nrows, ncols=2, width_ratios=(10000, 1))
        else:
            gs = plt.GridSpec(nrows=nrows, ncols=2, width_ratios=(10, 1))
        track_ax = fig.add_subplot(gs[-1, 0])
        sharex = None
        de_share = None

        plot_gene_transcripts(feat_res, ax=track_ax)
        dm_axs = {}
        de_axs = {}
        for i, nm in enumerate(self.dmr_comparison_groups):
            dm_axs[nm] = fig.add_subplot(gs[i, 0], sharex=sharex)
            sharex = dm_axs[nm]
            if self.de_res is not None:
                de_axs[nm] = fig.add_subplot(gs[i, 1], sharex=de_share, sharey=de_share)
                de_share = de_axs[nm]

        # keep tabs on all probes involved in DMRs (they may be outside of the region)
        probes_in_dmrs = set()

        for cid, d1 in this_dmrs.items():
            pc = this_clusters[cid]
            probes_in_dmrs.update(pc.pids)
            for pid, d2 in d1.items():
                this_ax = dm_axs[pid]
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

        # scatter showing individual probe values
        ix = (
                     (self.anno.CHR == ftr.chrom) & (self.anno.MAPINFO >= ftr.start) & (self.anno.MAPINFO <= ftr.stop)
             ) | (self.anno.index.isin(probes_in_dmrs))
        this_probe_ids = self.anno.loc[self.anno.index[ix].intersection(self.mdat.index), 'MAPINFO'].sort_values().index
        ymin = 0
        ymax = 0
        for nm, grp in self.dmr_comparison_groups.items():
            this_ax = dm_axs[nm]
            for col, x in self.mdat.loc[this_probe_ids, grp].iteritems():
                the_typ = self.meta.type[col]
                this_ax.scatter(
                    anno.MAPINFO[this_probe_ids],
                    x.values,
                    c=CELL_TYPE_COLOUR_MAP[the_typ],
                    marker=CELL_TYPE_MARKER_MAP[the_typ],
                    zorder=CELL_TYPE_ZORDER_MAP[the_typ],
                    alpha=CELL_TYPE_ALPHA_MAP[the_typ],
                    s=CELL_TYPE_SIZE_MAP[the_typ],
                    edgecolor='k',
                    linewidth=0.5
                )
                ymin = min(x.values.min(), ymin)
                ymax = max(x.values.max(), ymax)
                this_ax.set_ylabel(pid)

            if self.de_res is not None:
                this_de_ax = de_axs[pid]
                # bar chart to represent DE
                ix = de_res_s1[pid]['Gene Symbol'] == g
                if ix.sum() == 1:
                    this_de = de_res_s1[pid].loc[ix].squeeze()
                    this_de_ax.barh(
                        [0],
                        [this_de.logFC],
                        edgecolor='k' if this_de.FDR < de_params['fdr'] else 'none',
                        linewidth=1.,
                        color=consts.METHYLATION_DIRECTION_COLOURS['hyper'] if this_de.logFC > 0 else consts.METHYLATION_DIRECTION_COLOURS['hypo']
                    )
                elif ix.sum() > 1:
                    logger.warn("Gene %s. Found %d matching DE results in patient %s??", g, int(ix.sum()), pid)
                else:
                    logger.warn("Gene %s. No matching DE result for patient %s.", g, pid)


        dm_axs[consts.PIDS[-1]].xaxis.set_ticklabels([])
        [ax.set_ylim([ymin * 1.05, ymax * 1.05]) for ax in dm_axs.values()]

        if self.de_res is not None:
            plt.setp([t.yaxis for t in de_axs.values()], visible=False)
            de_axs[consts.PIDS[-1]].set_ylim([-1, 1])
            plt.setp(de_axs[consts.PIDS[-1]].xaxis.get_ticklabels(), rotation=90)
            de_axs[consts.PIDS[-1]].set_xlabel('DE logFC')


        gs.update(top=0.98, left=0.07, right=0.97, bottom=0.09, hspace=0.15, wspace=0.05)
        fig.savefig(os.path.join(outdir, "%s_mvalues.png" % g), dpi=300)

        if g in rois:
            for i, (x0, x1) in enumerate(rois[g]):
                track_ax.set_xlim([x0 - 200, x1 + 200])
                dm_axs.values()[0].set_xlim([x0 - 200, x1 + 200])
                fig.savefig(os.path.join(outdir, "%s_mvalues_roi_%d.png" % (g, i + 1)), dpi=300)

