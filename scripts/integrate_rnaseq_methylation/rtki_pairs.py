from methylation import dmr, process, plots
from load_data import methylation_array
from settings import DATA_DIR
from utils.output import unique_output_dir
import operator
import os
import numpy as np
import pandas as pd
import multiprocessing as mp
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3, venn2


def construct_contingency(x, y):
    return np.array([
        [((x < 0) & (y < 0)).sum(), ((x > 0) & (y < 0)).sum()],
        [((x < 0) & (y > 0)).sum(), ((x > 0) & (y > 0)).sum()],
    ])


def plot_n_region_heatmap(dat, n_arr, d_arr, ax=None, **kwargs):
    # heatmap plots for individual classes and combined
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    sns.heatmap(
        dat,
        cmap='RdBu_r',
        xticklabels=d_arr,
        yticklabels=n_arr,
        ax=ax,
        **kwargs
    )
    ax.set_xlabel("Maximum distance between probes")
    ax.set_ylabel("Minimum number of probes")
    ax.figure.tight_layout()

    return ax


if __name__ == '__main__':
    outdir = unique_output_dir("rtk1_de_dmr", reuse_empty=True)
    d_max = 200
    n_min = 4
    dm_min = 1.4  # minimum median delta M required to declare a cluster relevant
    alpha = 0.05

    patient_ids = ['018', '019', '031']
    ref_samples = ['insc', 'ensc_h9', 'ensc_fe']
    # n_jobs = mp.cpu_count()
    n_jobs = 4

    ## load all DE gene lists

    indir_de = os.path.join(DATA_DIR, 'rnaseq_de', 'rtk1')

    de = {}

    for p in patient_ids + ['all']:
        for ref in ref_samples:
            fn = os.path.join(indir_de, 'gbm-%s-%s.csv' % (ref, p))
            de[(ref, p)] = pd.read_csv(fn, header=0, index_col=0)

    ## compute DMR

    anno = methylation_array.load_illumina_methylationepic_annotation()
    b = methylation_array.gbm_nsc_methylationepic('swan')
    m = process.m_from_beta(b)

    # reduce anno down to probes in the data
    anno = anno.loc[anno.index.intersection(b.index)]

    # add merged class column to annotation
    dmr.add_merged_probe_classes(anno)

    # split gene symbols and store as a set
    anno.loc[:, 'UCSC_RefGene_Name'] = \
        anno.UCSC_RefGene_Name.str.split(';').apply(lambda x: set(x) if isinstance(x, list) else None)

    if False:
        # carry out a parameter sweep to see how the n_min and d_max parameters affect the number of valid regions

        d_max_arr = np.arange(100, 1100, 100)
        n_min_arr = np.arange(4, 11)
        cluster_sweep, nreg, npro = dmr.dmr_region_parameter_sweep(anno, d_max_arr=d_max_arr, n_min_arr=n_min_arr, n_jobs=n_jobs)

        # combine results for the distinct classes
        nreg_all = reduce(operator.add, nreg.values())
        npro_all = reduce(operator.add, npro.values())
        n_pro_tot = float(anno.shape[0])

        # heatmap plots
        # number of regions

        ax = plot_n_region_heatmap(nreg_all, n_min_arr, d_max_arr)
        cax = [a for a in ax.figure.get_axes() if a is not ax][0]
        cax.set_ylabel('Number of regions considered', rotation=270, labelpad=14)

        ax.figure.savefig(os.path.join(outdir, 'parameter_sweep_regions_total.png'), dpi=200)
        ax.figure.savefig(os.path.join(outdir, 'parameter_sweep_regions_total.pdf'))

        fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True)
        for i, k in enumerate(nreg.keys()):
            plot_n_region_heatmap(nreg[k], n_min_arr, d_max_arr, ax=axs[i], cbar=False)
            axs[i].set_title(k)
            plt.setp(axs[i].xaxis.get_ticklabels(), rotation=90)
            if i == 0:
                plt.setp(axs[i].yaxis.get_ticklabels(), rotation=0)
            else:
                axs[i].yaxis.label.set_visible(False)
        fig.tight_layout()

        fig.savefig(os.path.join(outdir, 'parameter_sweep_regions_by_class.png'), dpi=200)
        fig.savefig(os.path.join(outdir, 'parameter_sweep_regions_by_class.pdf'))

        ax = plot_n_region_heatmap(npro_all / n_pro_tot * 100, n_min_arr, d_max_arr, vmin=0, vmax=100)
        cax = [a for a in ax.figure.get_axes() if a is not ax][0]
        cax.set_ylabel('% probes retained', rotation=270, labelpad=14)

        ax.figure.savefig(os.path.join(outdir, 'parameter_sweep_probes_total.png'), dpi=200)
        ax.figure.savefig(os.path.join(outdir, 'parameter_sweep_probes_total.pdf'))

        fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True)
        for i, k in enumerate(npro.keys()):
            nt = anno.merged_class.str.contains(k).sum()
            plot_n_region_heatmap(npro[k] / nt * 100, n_min_arr, d_max_arr, ax=axs[i], cbar=False, vmin=0, vmax=100)
            axs[i].set_title(k)
            plt.setp(axs[i].xaxis.get_ticklabels(), rotation=90)
            if i == 0:
                plt.setp(axs[i].yaxis.get_ticklabels(), rotation=0)
            else:
                axs[i].yaxis.label.set_visible(False)
        fig.tight_layout()

        fig.savefig(os.path.join(outdir, 'parameter_sweep_probes_by_class.png'), dpi=200)
        fig.savefig(os.path.join(outdir, 'parameter_sweep_probes_by_class.pdf'))


    # carry out full relevance / significance analysis for a number of parameter values
    n_jobs = 2
    d_max_arr = [200, 500, 800]
    n_min_arr = [4, 6, 8]

    all_results = {}
    all_results_relevant = {}
    all_results_significant = {}

    for d in d_max_arr:
        for n in n_min_arr:
            print "d = %d; n = %d" % (d, n)

            clusters = dmr.identify_clusters(anno, n_min=n, d_max=d, n_jobs=n_jobs)
            test_results = {}
            test_results_relevant = {}
            test_results_significant = {}
            for sid in ['018', '019', '031']:
                samples = ('GBM%s' % sid, 'Dura%s' % sid)
                test_results[sid] = dmr.test_clusters(clusters, m, samples=samples, min_median_change=dm_min, n_jobs=n_jobs)
                test_results_relevant[sid] = dmr.mht_correction(test_results[sid], alpha=alpha)
                test_results_significant[sid] = dmr.filter_dictionary(
                    test_results_relevant[sid],
                    filt=lambda x: x['rej_h0'],
                    n_level=3
                )

            # add list of annotated genes to each relevant cluster
            for sid in ['018', '019', '031']:
                for (chr, cls, cid), attrs in dmr.dict_iterator(test_results_relevant[sid], n_level=3):
                    pids = clusters[chr][cls][cid]
                    genes = anno.loc[attrs['probes']].UCSC_RefGene_Name.dropna()
                    geneset = reduce(lambda x, y: x.union(y), genes, set())
                    attrs['genes'] = geneset

            all_results[(d, n)] = test_results
            all_results_relevant[(d, n)] = test_results_relevant
            all_results_significant[(d, n)] = test_results_significant

    raise StopIteration

    # 1: What is the joint distribution of methylation / mRNA fold change?
    # Get methylation level and DE fold change for linked genes (pairwise only)

    meth_de_joint = {}
    de_cols = ['genes', 'logFC', 'ensembl', 'direction', 'FDR', 'logCPM']
    meth_cols = ['me_genes', 'chr', 'me_cid', 'me_mediandelta', 'me_medianfc', 'me_fdr']

    for sid in patient_ids:
        print sid
        this_de = de[('insc', sid)]
        meth_de_joint[sid] = {}
        for (chr, cls, cid), attrs in dmr.dict_iterator(test_results_significant[sid], n_level=3):
            meth_de_joint[sid].setdefault(cls, pd.DataFrame(columns=de_cols + meth_cols))
            if len(attrs['genes']) == 0:
                # print "No annotated genes: (%s, %s, %d)" % (chr, cls, cid)
                continue
            try:
                # matching entry in DE
                de_match = this_de.loc[this_de.loc[:, 'genes'].isin(attrs['genes'])]
                if de_match.shape[0] == 0:
                    continue
                me_data = np.tile(
                    [
                        chr, cid, attrs['median_change'], attrs['median_fc'], attrs['padj']
                    ],
                    (de_match.shape[0], 1))
                me_data = np.concatenate(
                    (
                        np.reshape(de_match.genes.values, (de_match.shape[0], 1)),
                        me_data
                    ),
                    axis=1
                )
                me_match = pd.DataFrame(data=me_data, columns=meth_cols, index=de_match.index)
                this_match = pd.concat((de_match, me_match), axis=1)
            except Exception:
                print "Failed to add data: (%s, %s, %d)" % (chr, cls, cid)
                continue
            meth_de_joint[sid][cls] = pd.concat(
                (meth_de_joint[sid][cls], this_match), axis=0, ignore_index=True
            )
        meth_de_joint[sid]['all'] = pd.concat(meth_de_joint[sid].values(), axis=0, ignore_index=True)

    # Generate table giving the number of overlaps in each patient and cluster class
    # this includes the number of absolute overlaps AND the number of unique overlaps
    the_cols = ['DE genes', 'DMR', 'DMR genes', 'overlaps', 'unique overlaps']
    the_cols += reduce(lambda x, y: x + y, [['%s' % t, '%s_unique' % t] for t in dmr.CLASSES], [])
    de_dmr_matches = pd.DataFrame(
        columns=the_cols,
        index=pd.Index(patient_ids, name='patient'),
    )

    for sid in patient_ids:
        this_all = meth_de_joint[sid]['all']
        n_de = de[('insc', sid)].shape[0]
        n_dmr = len(list(dmr.dict_iterator(test_results_significant[sid], n_level=3)))
        n_dmr_genes = len(
            reduce(
                lambda x, y: x.union(y), [t[1]['genes'] for t in dmr.dict_iterator(test_results_significant[sid], n_level=3)], set()
            )
        )
        n_overlaps = meth_de_joint[sid]['all'].shape[0]
        n_overlaps_unique = meth_de_joint[sid]['all'].me_genes.unique().shape[0]
        this_datum = [n_de, n_dmr, n_dmr_genes, n_overlaps, n_overlaps_unique]
        for cls in dmr.CLASSES:
            this_datum.append(meth_de_joint[sid][cls].shape[0])
            this_datum.append(meth_de_joint[sid][cls].me_genes.unique().shape[0])
        de_dmr_matches.loc[sid] = this_datum

    # Scatter plots
    for sid in ['018', '019', '031']:
        fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        for i, cls in enumerate(['all', 'tss', 'gene', 'island']):
            ax = axs.flat[i]

            # get values for ALL DMR clusters of this class
            x = meth_de_joint[sid][cls].loc[:, 'logFC'].astype(float)
            y = meth_de_joint[sid][cls].loc[:, 'me_mediandelta'].astype(float)

            # contingency table for Fisher's exact test
            conting = construct_contingency(x, y)
            logodds, fisherp = stats.fisher_exact(conting)

            # get values for DMR clusters that are ONLY in this class (no overlaps)
            if cls == 'all':
                print "%s - %s p = %.3e" % (
                    sid, cls, fisherp
                )
                ax.scatter(x, y, c='k')
                ax.axhline(0, c=0.4 * np.ones(3))
                ax.axvline(0, c=0.4 * np.ones(3))
                ttl = "%s; p={0}" % cls
                if fisherp < 0.001:
                    ttl = ttl.format('%.2e') % fisherp
                else:
                    ttl = ttl.format('%.3f') % fisherp

            else:
                cid_other = set(
                    np.unique(
                        np.concatenate([meth_de_joint[sid][t].me_cid.values for t in dmr.CLASSES.difference({cls,})])
                    )
                )
                xu = x.loc[~meth_de_joint[sid][cls].loc[:, 'me_cid'].isin(cid_other)].astype(float)
                yu = y.loc[~meth_de_joint[sid][cls].loc[:, 'me_cid'].isin(cid_other)].astype(float)
                contingu = construct_contingency(xu, yu)
                logoddsu, fisherpu = stats.fisher_exact(contingu)

                print "%s - %s p = %.3e (incl overlaps), p = %.3e (excl overlaps)" % (
                    sid, cls, fisherp, fisherpu
                )

                ax.scatter(x, y, c='k')
                ax.scatter(xu, yu, c='b')
                ax.axhline(0, c=0.4 * np.ones(3))
                ax.axvline(0, c=0.4 * np.ones(3))

                ttl = "%s; p={0}; unique p={0}" % cls

                if fisherp < 0.001:
                    ttl = ttl.format('%.2e') % (fisherp, fisherpu)
                else:
                    ttl = ttl.format('%.3f') % (fisherp, fisherpu)

            ax.set_title(ttl)

        fig.text(0.5, 0.04, 'RNASeq DE logFC', ha='center')
        fig.text(0.04, 0.5, 'EPIC DMR median delta M', va='center', rotation='vertical')
        fig.savefig(os.path.join(outdir, "de_vs_dmr_%s.png" % sid), dpi=200)
        fig.savefig(os.path.join(outdir, "de_vs_dmr_%s.pdf" % sid))

    # 2: To what extent do the same genes appear in all RTK 1 samples?
    # Venn diagram
    all_genes = reduce(
        lambda x, y: x.union(y),
        (set(meth_de_joint[sid]['all'].genes.unique()) for sid in ['018', '019', '031']),
        set()
    )

    fig, axs = plt.subplots(ncols=3, figsize=(8, 3.2))
    f = open(os.path.join(outdir, "core_genes_de_dmr.txt"), 'wb')

    sids = ('018', '019', '031')
    venn_counts = {}
    venn_sets = {}
    for i, cls in enumerate(dmr.CLASSES):
        this_genecount = {}
        this_geneset = {}

        # all
        for j in range(1, 8):
            bn = "{0:03b}".format(j)
            this_intersection = set(all_genes)
            for k in range(3):
                if bn[k] == '1':
                    this_intersection = this_intersection.intersection(
                        meth_de_joint[sids[k]][cls].genes.unique()
                    )
            this_genecount[bn] = len(this_intersection)
            this_geneset[bn] = list(this_intersection)
        venn_counts[cls] = this_genecount
        venn_sets[cls] = this_geneset
        venn = venn3(subsets=venn_counts[cls], set_labels=sids, ax=axs[i])
        axs[i].set_title(cls)
        print "%s core genes: %s" % (cls, ', '.join(this_geneset['111']))
        f.write("%s core genes: %s\n" % (cls, ', '.join(this_geneset['111'])))
    core_all = set(venn_sets['tss']['111']).intersection(venn_sets['gene']['111']).intersection(venn_sets['island']['111'])
    print "Core genes shared across all classes: %s" % ', '.join(list(core_all))
    f.write("Core genes shared across all classes: %s\n" % ', '.join(list(core_all)))
    f.close()

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_and_de_overlap.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "dmr_and_de_overlap.pdf"))


    # TODO: consider running a large parameter sweep (possibly on Apocrita)
    # matrix of (n, d) parameters
    # for each parameter set and for each sample, store clusters, test_results
    # analyse the number of clusters, relevant and significant. Why is 031 so weird?!