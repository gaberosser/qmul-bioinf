from matplotlib import pyplot as plt
from plotting.utils import COLOUR_BREWERS
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


def illustrate_regions(anno, regions, cls_list, chr, loc_from, loc_to, cmap=None,
                       alpha_classed=0.7, ax=None):
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

    this_dat = anno.loc[anno.CHR == chr].sort_values(by='MAPINFO', axis=0)
    idx = (loc_from <= this_dat.MAPINFO) & (this_dat.MAPINFO < loc_to)
    t = this_dat.loc[idx]

    sq_buffer = 0.2
    for grp, c in cmap.items():
        if grp in regions:
            for p_arr in regions[grp]:
                # get coords
                xs = t.MAPINFO.loc[p_arr]
                min_x = xs.min()
                max_x = xs.max()
                sq_x = [min_x, min_x, max_x, max_x, min_x]
                sq_y = [-sq_buffer, 1 + sq_buffer, 1 + sq_buffer, -sq_buffer, -sq_buffer]
                ax.plot(sq_x, sq_y, c=c, ls='-', lw=1.5, label=None, alpha=alpha_classed)
