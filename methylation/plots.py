from matplotlib import pyplot as plt
from matplotlib import patches, collections
from plotting.common import COLOUR_BREWERS
from plotting import venn
import dmr
from utils import dictionary
import os
import seaborn as sns


def illustrate_probes(anno, cls, chr, loc_from, loc_to, cmap=None,
                      alpha_classed=0.7, alpha_unclassed=0.3, ax=None,
                      regions=None):
    """
    Plot showing the probes arranged along a section of the genome, coloured by class
    :param anno: Annotation DataFrame as loaded by load_data.methylation_array.load_illumina_methylationepic_annotation
    :param cls: Series giving the class of each probe, separated by a semicolon when more than one applies.
    Must have an index compatible with anno (i.e. probe names).
    :param chr: Chromosome
    :param loc_from:
    :param loc_to:
    :param cmap: The colours used for the classes. If none, a default option is used.
    :param regions: If provided, this is a dictionary. The keys match the classes used in `cls`. The values are lists
    of lists of probes (each nested list is a region).
    :return:
    """
    if ax is None:
        fig = plt.figure(figsize=(10, 3))
        ax = fig.add_subplot(111)

    if cmap is None:
        # need to find classes ourselves
        t = cls.str.split(';').apply(set)
        all_classes = reduce(lambda x, y: x.union(y), t)
        if '' in all_classes:
            all_classes.remove('')
        try:
            cols = COLOUR_BREWERS[len(all_classes)]
            cmap = dict(zip(all_classes, cols))
        except KeyError:
            raise ValueError(
                "Too many classes (%d) for automatic cmap generation. Please provide cmap." % len(all_classes)
            )

    this_dat = anno.loc[anno.CHR == chr].sort_values(by='MAPINFO', axis=0)
    idx = (loc_from <= this_dat.MAPINFO) & (this_dat.MAPINFO < loc_to)

    strand_y = this_dat.Strand.where(this_dat.Strand == 'F', 0)
    strand_y.loc[strand_y == 'F'] = 1
    strand_y = strand_y.astype(int)

    colour_no_role = [0.5, 0.5, 0.5]

    t = this_dat.loc[idx]
    for grp, c in cmap.items():
        cidx = cls.str.contains(grp)
        x = t.MAPINFO[idx & cidx].values
        y = strand_y[idx & cidx].values
        ax.scatter(x, y, c=c, marker='|', label=grp, alpha=alpha_classed, zorder=3)

    if regions is not None:
        illustrate_regions(anno, regions, cls, chr, loc_from, loc_to, cmap=cmap, alpha_classed=alpha_classed, ax=ax)

    cidx = cls == ''
    x = t.MAPINFO[idx & cidx].values
    y = strand_y[idx & cidx].values
    ax.scatter(x, y, c=colour_no_role, marker='|', alpha=alpha_unclassed, label='none', zorder=2)
    ax.plot([loc_from, loc_to], [0, 0], 'k-', zorder=1, alpha=alpha_unclassed)
    ax.plot([loc_from, loc_to], [1, 1], 'k-', zorder=1, alpha=alpha_unclassed)
    ax.legend(loc='center right', frameon=True, facecolor='w')
    ax.set_ylim([-.3, 1.3])
    ax.yaxis.set_ticks([0, 1])
    ax.yaxis.set_ticklabels(['R', 'F'])
    ax.set_xlabel('Chromosomal coordinate')
    ax.set_title('Chromosome {}: {:,} - {:,}'.format(chr, loc_from, loc_to))

    return ax


def patch_from_cluster(anno, probes, y, height, **kwargs):
    xs = anno.loc[probes, 'MAPINFO']
    min_x = xs.min()
    max_x = xs.max()
    xy = (min_x, y)
    width = max_x - min_x
    return patches.Rectangle(xy, width, height, **kwargs)


def illustrate_regions(anno, regions, cls_list, cmap=None,
                       alpha=0.7, ax=None, ylim=None, linestyle='-'):
    """
    :param regions: From `get_clusters_by_location`,this is a nested dict.
    """
    if ax is None:
        ax = plt.gca()

    if cmap is None:
        # need to find classes ourselves
        try:
            cols = COLOUR_BREWERS[len(cls_list)]
            cmap = dict(zip(cls_list, cols))
        except KeyError:
            raise ValueError(
                "Too many classes (%d) for automatic cmap generation. Please provide cmap." % len(cls_list)
            )

    if ylim is None:
        sq_buffer = 0.2
        ymin = -sq_buffer
        ymax = 1 + sq_buffer
    else:
        ymin, ymax = ylim

    height = ymax - ymin

    for grp, c in cmap.items():
        ptchs = []
        if grp in regions:
            for p_arr in regions[grp].values():
                ptchs.append(patch_from_cluster(anno, p_arr, ymin, height))
            coll = collections.PatchCollection(
                ptchs,
                alpha=alpha,
                linewidths=1.5,
                linestyles=linestyle,
                facecolors='none',
                edgecolors=c
            )
            ax.add_collection(coll)
            print "Added collection of %d patches with colour %s" % (len(ptchs), c)


def venn_dmr_counts(
        dmr_results,
        pids=None,
        outdir=None,
        figname='dmr_venn',
        comparisons=('matched', 'gibco'),
        set_labels=('Isogenic', 'Reference')
):
    """
    For each supplied patient (top level of keys in the supplied results), generate 3 Venn diagrams
    - all DMR count
    - hypermethylated DMR counts
    - hypomethylated DMR counts
    In each case, splitting according to `comparisons`
    :param dmr_results: Dictionary of results, e.g. from the results (results_significant) attribute of DmrResults.
    :param outdir:
    :param comparisons: Iterable giving the comparisons to use.
    :param set_labels: Iterable giving the set labels to use for the plot. Coresponds to comparisons.
    :return:
    """
    if len(set_labels) != len(comparisons):
        raise AttributeError("set_labels must be the same length as comparisons")
    if pids is None:
        pids = dmr_results.keys()

    # isogenic vs reference in each patient, all, hyper and hypomethylated
    fig, axs = plt.subplots(3, len(pids), figsize=(11, 7), num=figname)
    for i, pid in enumerate(pids):
        inputs = [
            dmr_results[pid][k] for k in comparisons
        ]

        inputs_all = [set(x.keys()) for x in inputs]

        inputs_up = [set(
            [k for k, v in x.items() if v['median_change'] > 0]
        ) for x in inputs]
        inputs_down = [set(
            [k for k, v in x.items() if v['median_change'] < 0]
        ) for x in inputs]

        sl = set_labels if i == 0 else None

        venn.venn_diagram(set_labels=sl, ax=axs[0, i], *inputs_all)
        axs[0, i].set_title("GBM%s" % pid)

        venn.venn_diagram(set_labels=sl, ax=axs[1, i], *inputs_up)
        venn.venn_diagram(set_labels=sl, ax=axs[2, i], *inputs_down)

    fig.tight_layout()
    if outdir is not None:
        fig.savefig(os.path.join(outdir, "%s.png" % figname), dpi=200)
        fig.savefig(os.path.join(outdir, "%s_venn.pdf" % figname))


def dmr_overlap(
        dmr_results,
        pids=None,
        comparisons=('matched', 'gibco'),
        comparison_titles=('Isogenic', 'Reference'),
        outdir=None,
        figname='dmr_individual_overlap'
):
    """
    Plot a Venn diagram showing the overlap of DMRs in the the individuals. One plot per comparison (default)
    behaviour is to include two comparisons: isogenic and reference.
    This can be limited to a single DMR probe class, or (default) include all.
    :param dmr_results:
    :param pids:
    :param comparisons:
    :param comparison_titles:
    :param outdir:
    :param figname:
    :return:
    """
    if len(comparisons) != len(comparison_titles):
        raise AttributeError("Length of comparisons and comparison_titles must be equal.")
    if pids is None:
        pids = dmr_results.keys()
    fig, axs = plt.subplots(1, len(comparisons), num=figname)

    for i, cmp in enumerate(comparisons):
        inputs = [
            dmr_results[pid][cmp] for pid in pids
        ]
        # n_level = 3
        # if probe_class is not None:
        #     inputs = [dmr.dict_by_sublevel(t, 2, probe_class) for t in inputs]
        #     n_level = 2
        #     hasher = lambda x: tuple(x)
        # else:
        #     hasher = lambda x: (x[0], x[2])
        #
        # inputs = [
        #     set([hasher(t) for t, _ in dmr.dict_iterator(x, n_level=n_level)]) for x in inputs
        # ]
        venn.venn_diagram(set_labels=pids, ax=axs[i], *inputs)
        axs[i].set_title(comparison_titles[i])

    fig.tight_layout()
    if outdir is not None:
        fig.savefig(os.path.join(outdir, "%s.png" % figname), dpi=200)
        fig.savefig(os.path.join(outdir, "%s.pdf" % figname))


def dmr_cluster_count_by_class(
        blocks,
        show_labels=True,
        set_labels=None,
        outdir=None,
        figname="dmr_cluster_count_by_class",
        ax=None
):
    """
    Venn diagram showing the distribution of clusters amongst the classes for a single result.
    :param blocks: dictionary, keyed by class, values are iterables of IDs
    :param outdir: If supplied, write plot to a file.
    :param figname:
    :return:
    """
    created = False
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        created = True

    # get classes
    # all_classes = dmr_results.classes

    sl, blocks = zip(*blocks.items())

    # blocks = []
    # sl = []
    # for cls in all_classes:
    #     d = dmr_results[cls]
    #     blocks.append(d.keys())
    #     sl.append(cls)

    if set_labels is None:
        set_labels = sl

    if not show_labels:
        set_labels = None

    venn.venn_diagram(*blocks, ax=ax, set_labels=set_labels)

    if created:
        fig.tight_layout()

    if outdir is not None:
        ax.figure.savefig(os.path.join(outdir, "%s.png" % figname), dpi=200)
        ax.figure.savefig(os.path.join(outdir, "%s.pdf" % figname))
    return ax


def dmr_cluster_count_array(
        dmr_results,
        comparisons=('matched', 'gibco'),
        comparison_labels=True,
        outdir=None,
        figname='dmr_sign_cluster_count_array',
        figsize=None
):
    """
    :param dmr_results: Full dictionary output by compute_dmr
    :param comparisons: If supplied, limit comparisons to this/these.
    """
    n_pid = len(dmr_results)
    n_cmp = len(comparisons)

    if comparison_labels is not None and comparison_labels != True and len(comparison_labels) != n_cmp:
        raise AttributeError("Length of comparisons and comparison_labels must be equal.")
    elif comparison_labels is True:
        comparison_labels = comparisons

    if figsize is None:
        figsize = (n_pid ** .9 * 2.4, n_cmp ** .9 * 2.7)
    fig, axs = plt.subplots(n_cmp, n_pid, figsize = figsize)
    show_labels = True
    for j, pid in enumerate(sorted(dmr_results.keys())):
        for i, cmp in enumerate(comparisons):
            ax = axs[i, j]
            dmr_cluster_count_by_class(
                dmr_results[pid][cmp],
                ax=ax,
                show_labels=show_labels,
            )
            if comparison_labels and j == 0:
                ax.text(0.1, 0.5, cmp, verticalalignment='center', horizontalalignment='right', transform=ax.transAxes)
            show_labels = False

    fig.tight_layout(rect=(0.03, 0.03, 1, 1), pad=0.)

    # add text labels
    for j, pid in enumerate(sorted(dmr_results.keys())):
        fig.text((j + 1.) / (n_pid + 1.), 0.02, pid)

    if outdir is not None:
        fig.savefig(os.path.join(outdir, "%s.png" % figname), dpi=200)
        fig.savefig(os.path.join(outdir, "%s.pdf" % figname))