from matplotlib import pyplot as plt, patches, collections as plt_collections

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