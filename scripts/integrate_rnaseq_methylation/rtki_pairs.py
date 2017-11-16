from methylation import dmr, process, plots
from load_data import methylation_array
from settings import GIT_LFS_DATA_DIR
from utils.output import unique_output_dir
import operator
import os
import numpy as np
import pandas as pd
from scipy import stats
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
    d_max = 400
    n_min = 6
    dm_min = 1.4  # minimum median delta M required to declare a cluster relevant
    alpha = 0.05

    patient_ids = ['018', '019', '031']
    # n_jobs = mp.cpu_count()
    n_jobs = 4

    ## load all DE gene lists

    indir_de = os.path.join(GIT_LFS_DATA_DIR, 'rnaseq_de', 'rtk1', 'insc_h9nsc')
    de = {}

    for p in patient_ids + ['all']:
        fn = os.path.join(indir_de, 'gbm-insc-ensc-%s.csv' % p)
        this_de_insc_only = pd.read_csv(fn, header=0, index_col=None)
        in_insc = ~this_de_insc_only.iloc[:, 1].isnull()
        in_ensc = ~this_de_insc_only.iloc[:, 7].isnull()

        # DE genes in iNSC comparison only
        de[(p, 'insc_only')] = this_de_insc_only.loc[in_insc & ~in_ensc].iloc[:, :6]
        # DE genes in both comparisons
        # here we use the logFC etc from the iNSC comparison, since this is what we're interested in
        de[(p, 'insc_and_h9')] = this_de_insc_only.loc[in_insc & in_ensc].iloc[:, :6]
        # DE genes in H9 comparison only
        de[(p, 'h9_only')] = this_de_insc_only.loc[~in_insc & in_ensc].iloc[:, 6:]
        de[(p, 'h9_only')].columns = this_de_insc_only.columns[:6]

    ## compute DMR

    anno = methylation_array.load_illumina_methylationepic_annotation()
    b, me_meta = methylation_array.gbm_rtk1_and_paired_nsc(norm_method='swan')
    m = process.m_from_beta(b)

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(b.index)
    anno = anno.loc[common_probes]
    b = b.loc[common_probes]
    m = m.loc[common_probes]

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
    # TODO: this is very SLOW. I have saved the results using dill.

    n_jobs = 12
    # d_max_arr = [200, 500, 800]
    d_max_arr = [500]
    # n_min_arr = [4, 6, 8]
    n_min_arr = [4]

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
                ## FIXME: the function  mht_correction now modifies in place; this will error
                test_results_relevant[sid] = dmr.mht_correction(test_results[sid], alpha=alpha)
                test_results_significant[sid] = dmr.filter_dictionary(
                    test_results_relevant[sid],
                    filt=lambda x: x['rej_h0'],
                    n_level=3
                )

            # add list of annotated genes to all clusters
            for sid in ['018', '019', '031']:
                for (chr, cls, cid), attrs in dmr.dict_iterator(test_results[sid], n_level=3):
                    pids = attrs['probes']
                    genes = anno.loc[attrs['probes']].UCSC_RefGene_Name.dropna()
                    geneset = reduce(lambda x, y: x.union(y), genes, set())
                    attrs['genes'] = geneset

            all_results[(d, n)] = test_results
            all_results_relevant[(d, n)] = test_results_relevant
            all_results_significant[(d, n)] = test_results_significant

    # generate table with the number of proposed clusters, relevant clusters and significant clusters for each
    # sample at each of the parameter values tested
    cols = (
        'd_max', 'n_min', 'sample', 'clusters_proposed', 'clusters_relevant', 'clusters_significant',
        'genes_proposed', 'genes_relevant', 'genes_significant'
    )
    table_cluster_numbers = pd.DataFrame(columns=cols)

    def count_genes(res, sid):
        the_genes = reduce(
            lambda x, y: x.union(y),
            [t[1]['genes'] for t in dmr.dict_iterator(res[(d, n)][sid], n_level=3)],
            set()
        )
        return len(the_genes)

    for d in d_max_arr:
        for n in n_min_arr:
            ncl = len(list(
                dmr.dict_iterator(all_results[(d, n)][patient_ids[0]], n_level=3)
            ))
            ng = count_genes(all_results, patient_ids[0])

            for sid in ['018', '019', '031']:
                ncl_re = len(list(
                    dmr.dict_iterator(all_results_relevant[(d, n)][sid], n_level=3)
                ))
                ncl_si = len(list(
                    dmr.dict_iterator(all_results_significant[(d, n)][sid], n_level=3)
                ))
                ng_re = count_genes(all_results_relevant, sid)
                ng_si = count_genes(all_results_significant, sid)
                this_row = pd.Series({
                    'd_max': d,
                    'n_min': n,
                    'sample': sid,
                    'clusters_proposed': ncl,
                    'clusters_relevant': ncl_re,
                    'clusters_significant': ncl_si,
                    'genes_proposed': ng,
                    'genes_relevant': ng_re,
                    'genes_significant': ng_si
                })
                table_cluster_numbers = table_cluster_numbers.append(this_row, ignore_index=True)

    # investigate issue with 031 and (200, 4): very low number of significant clusters
    # start by looking at the distribution of adjusted pvalues in each patient
    if False:

        pvals = {}
        padj = {}
        dupe_pvals = {}
        dupe_padj = {}
        for n in n_min_arr:
            for d in d_max_arr:
                pvals[(d, n)] = {}
                padj[(d, n)] = {}
                dupe_pvals[(d, n)] = {}
                dupe_padj[(d, n)] = {}
                for sid in patient_ids:
                    already_seen = set()
                    pvals[(d, n)][sid] = []
                    padj[(d, n)][sid] = []
                    dupe_pvals[(d, n)][sid] = []
                    dupe_padj[(d, n)][sid] = []
                    g = dmr.dict_iterator(all_results_relevant[(d, n)][sid], n_level=3)
                    for k, attrs in g:
                        dupe_pvals[(d, n)][sid].append(attrs['pval'])
                        dupe_padj[(d, n)][sid].append(attrs['padj'])
                        probes = tuple(attrs['probes'])
                        if probes not in already_seen:
                            pvals[(d, n)][sid].append(attrs['pval'])
                            padj[(d, n)][sid].append(attrs['padj'])
                            already_seen.add(probes)

        # histogram of (adj) pval

        cutoff = np.log10(0.05)
        x0 = np.linspace(-4, cutoff, 40)
        x1 = np.arange(cutoff, 0., x0[1] - x0[0])
        bins = np.concatenate((x0[:-1], x1))

        # we need to use (200, 4) here to demonstrate the 'issue' being investigated
        (d, n) = (500, 4)

        fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)
        for j, pid in enumerate(patient_ids):
            axs[j].hist(np.log10(dupe_padj[(d, n)][pid]), bins, label='Incl duplicates')
            axs[j].hist(np.log10(padj[(d, n)][pid]), bins, label='Excl duplicates')
            axs[j].axvline(cutoff, ls='--', c='k', label='FDR cutoff')
            axs[j].set_title(pid)

        axs[0].legend(loc='upper left')
        axs[-1].set_xlabel('log10 adjusted pvalue')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'padj_distribution_%d_%d.png' % (d, n)), dpi=200)
        fig.savefig(os.path.join(outdir, 'padj_distribution_%d_%d.pdf' % (d, n)))

        fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)
        for j, pid in enumerate(patient_ids):
            axs[j].hist(np.log10(dupe_pvals[(d, n)][pid]), bins, label='Incl duplicates')
            axs[j].hist(np.log10(pvals[(d, n)][pid]), bins, label='Excl duplicates')
            axs[j].axvline(cutoff, ls='--', c='k', label='FDR cutoff')
            axs[j].set_title(pid)

        axs[0].legend(loc='upper left')
        axs[-1].set_xlabel('log10 unadjusted pvalue')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'pval_distribution_%d_%d.png' % (d, n)), dpi=200)
        fig.savefig(os.path.join(outdir, 'pval_distribution_%d_%d.pdf' % (d, n)))

        padj_max = 0
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for j, pid in enumerate(patient_ids):
            this_pval = np.array(dupe_pvals[(d, n)][pid])
            this_padj = np.array(dupe_padj[(d, n)][pid])
            sort_idx = np.argsort(this_pval)
            this_pval = this_pval[sort_idx]
            this_padj = this_padj[sort_idx]
            first_idx = np.where(this_pval > alpha)[0][0]
            padj_max = max(padj_max, this_padj[first_idx])
            ax.scatter(this_pval, this_padj, label=pid)
            ax.axhline(alpha, ls='--', c='k', label='FDR cutoff' if j == 0 else None)

        ax.legend(loc='upper left')
        ax.set_xlim([0, alpha])
        ax.set_ylim([0, padj_max])
        ax.set_xlabel("Unadjusted pvalue")
        ax.set_ylabel("Adjusted pvalue")

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'pval_vs_padj_%d_%d.png' % (d, n)), dpi=200)
        fig.savefig(os.path.join(outdir, 'pval_vs_padj_%d_%d.pdf' % (d, n)))

    """
    The histograms show a clear 'spiky' pattern, indicative of the pvalues coming from a finite set of values.
    This makes sense: the Mann-Whitney U value depends on the ordering and the number of probes in a cluster. There
    are a limited number of possibilities.
    """

    # Match DMR and DE.
    # We do this based on gene symbol, as it is (unfortunately) the only information contained in the EPIC annotation.
    # We are interested in two comparisons:
    # a) DE (in GBM vs iNSC but not in GBM vs H9)
    # b) DE (in GBM vs iNSC and in GBM vs H9)

    # 1: What is the joint distribution of methylation / mRNA fold change?
    # Get methylation level and DE fold change for linked genes (pairwise only)

    def compute_joint_de_dmr(this_test_results, this_de):
        res = {}

        for sid in patient_ids:
            print sid
            res[sid] = {}

            de_cols = ['genes', 'logFC', 'ensembl', 'direction', 'FDR', 'logCPM']
            meth_cols = ['me_genes', 'chr', 'me_cid', 'me_mediandelta', 'me_medianfc', 'me_fdr']

            for (chr, cls, cid), attrs in dmr.dict_iterator(this_test_results[sid], n_level=3):
                res[sid].setdefault(cls, pd.DataFrame(columns=de_cols + meth_cols))

                if len(attrs['genes']) == 0:
                    continue

                try:
                    # matching entry in DE
                    de_match = this_de[sid].loc[this_de[sid].loc[:, 'genes'].isin(attrs['genes'])]

                    if de_match.shape[0] > 0:
                        # form the DMR data block
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
                        res[sid][cls] = pd.concat(
                            (res[sid][cls], this_match), axis=0, ignore_index=True
                        )
                except Exception as exc:
                    print "Failed to add data (iNSC only): (%s, %s, %d)" % (chr, cls, cid)
                    print repr(exc)
                    continue

            res[sid]['all'] = pd.concat(res[sid].values(), axis=0, ignore_index=True)

        return res

    this_de_insc_only = dict([(sid, de[(sid, 'insc_only')]) for sid in patient_ids])
    this_de_insc_h9 = dict([(sid, de[(sid, 'insc_and_h9')]) for sid in patient_ids])
    meth_de_joint_insc_only = compute_joint_de_dmr(test_results_significant, this_de_insc_only)
    meth_de_joint_insc_h9 = compute_joint_de_dmr(test_results_significant, this_de_insc_h9)

    # Generate table giving the number of overlaps in each patient and cluster class
    # this includes the number of absolute overlaps AND the number of unique overlaps
    the_cols = ['DE genes', 'DMR', 'DMR genes', 'overlaps', 'unique overlaps']
    the_cols += reduce(lambda x, y: x + y, [['%s' % t, '%s_unique' % t] for t in dmr.CLASSES], [])
    de_dmr_matches_insc_only = pd.DataFrame(
        columns=the_cols,
        index=pd.Index(patient_ids, name='patient'),
    )
    de_dmr_matches_insc_h9 = pd.DataFrame.copy(de_dmr_matches_insc_only)

    def n_overlap_datum(sid, meth_de, this_de):
        n_de = this_de.shape[0]
        n_dmr = len(list(dmr.dict_iterator(test_results_significant[sid], n_level=3)))
        n_dmr_genes = len(
            reduce(
                lambda x, y: x.union(y), [t[1]['genes'] for t in dmr.dict_iterator(test_results_significant[sid], n_level=3)], set()
            )
        )
        n_overlaps = meth_de[sid]['all'].shape[0]
        n_overlaps_unique = meth_de[sid]['all'].me_genes.unique().shape[0]
        this_datum = [n_de, n_dmr, n_dmr_genes, n_overlaps, n_overlaps_unique]
        for cls in dmr.CLASSES:
            this_datum.append(meth_de[sid][cls].shape[0])
            this_datum.append(meth_de[sid][cls].me_genes.unique().shape[0])
        return this_datum

    for sid in patient_ids:
        de_dmr_matches_insc_only.loc[sid] = n_overlap_datum(sid, meth_de_joint_insc_only, de[(sid, 'insc_only')])
        de_dmr_matches_insc_h9.loc[sid] = n_overlap_datum(sid, meth_de_joint_insc_h9, de[(sid, 'insc_and_h9')])


    def scatter_plot_dmr_de(meth_de, fig_filestem):
        for sid in ['018', '019', '031']:
            fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
            for i, cls in enumerate(['all', 'tss', 'gene', 'island']):
                ax = axs.flat[i]

                # get values for ALL DMR clusters of this class
                x = meth_de[sid][cls].loc[:, 'logFC'].astype(float)
                y = meth_de[sid][cls].loc[:, 'me_mediandelta'].astype(float)

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
                            np.concatenate([meth_de[sid][t].me_cid.values for t in dmr.CLASSES.difference({cls, })])
                        )
                    )
                    xu = x.loc[~meth_de[sid][cls].loc[:, 'me_cid'].isin(cid_other)].astype(float)
                    yu = y.loc[~meth_de[sid][cls].loc[:, 'me_cid'].isin(cid_other)].astype(float)
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
            fig.savefig("%s_%s.png" % (fig_filestem, sid), dpi=200)
            fig.savefig("%s_%s.png" % (fig_filestem, sid), dpi=200)

    print "*** Genes that match and are DE in GBM vs iNSC only ***"
    scatter_plot_dmr_de(meth_de_joint_insc_only, os.path.join(outdir, "de_vs_dmr_insc_only"))
    print "*** Genes that match and are DE in GBM vs iNSC AND GBM vs H9 ***"
    scatter_plot_dmr_de(meth_de_joint_insc_h9, os.path.join(outdir, "de_vs_dmr_insc_h9"))

    # 2: To what extent do the same genes appear in all RTK 1 samples?
    def venn_diagram_and_core_genes(meth_de, text_file, fig_file):
        all_genes = reduce(
            lambda x, y: x.union(y),
            (set(meth_de[sid]['all'].genes.unique()) for sid in patient_ids),
            set()
        )

        fig, axs = plt.subplots(ncols=3, figsize=(8, 3.2))
        f = open(text_file, 'wb')

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
                            meth_de[patient_ids[k]][cls].genes.unique()
                        )
                this_genecount[bn] = len(this_intersection)
                this_geneset[bn] = list(this_intersection)
            venn_counts[cls] = this_genecount
            venn_sets[cls] = this_geneset
            venn = venn3(subsets=venn_counts[cls], set_labels=patient_ids, ax=axs[i])
            axs[i].set_title(cls)
            print "%s core genes: %s" % (cls, ', '.join(this_geneset['111']))
            f.write("%s core genes: %s\n" % (cls, ', '.join(this_geneset['111'])))
        core_all = set(venn_sets['tss']['111']).intersection(venn_sets['gene']['111']).intersection(
            venn_sets['island']['111'])
        print "Core genes shared across all classes: %s" % ', '.join(list(core_all))
        f.write("Core genes shared across all classes: %s\n" % ', '.join(list(core_all)))
        f.close()

        fig.tight_layout()
        fig.savefig("%s.png" % fig_file, dpi=200)
        fig.savefig("%s.pdf" % fig_file)

        return venn_sets


    print "*** GBM vs iNSC and NOT GBM vs H9 Venn overlaps ***"
    venn_insc_only = venn_diagram_and_core_genes(
        meth_de_joint_insc_only,
        os.path.join(outdir, "core_genes_de_dmr_insc_only.txt"),
        os.path.join(outdir, "dmr_and_de_overlap_insc_only")
    )
    print "*** GBM vs iNSC and GBM vs H9 Venn overlaps ***"
    venn_insc_h9 = venn_diagram_and_core_genes(
        meth_de_joint_insc_h9,
        os.path.join(outdir, "core_genes_de_dmr_insc_h9.txt"),
        os.path.join(outdir, "dmr_and_de_overlap_insc_h9")
    )

    # analysing the length distribution of the clusters
    # this varies depending on the parameters to define clusters
    cluster_size_dist = {}
    cluster_size_max = 40

    for n in n_min_arr:
        for d in d_max_arr:
            g = dmr.dict_iterator(all_results[(d, n)], n_level=4)
            probesets = set([tuple(attrs['probes']) for _, attrs in g])
            cluster_size_dist[(d, n)] = np.histogram([len(t) for t in probesets], bins=range(cluster_size_max))

    fig, axs = plt.subplots(nrows=2, sharex=True)
    [axs[0].plot(t[1][:-1], t[0], label="d_max = %d" % k[0]) for k, t in cluster_size_dist.items() if k[1] == 4]
    axs[0].legend(loc='upper right')
    [axs[1].plot(t[1][:-1], t[0], label="n_min = %d" % k[1]) for k, t in cluster_size_dist.items() if k[0] == 200]
    axs[1].legend(loc='upper right')
    axs[1].set_xlabel('Number in cluster')
    fig.savefig(os.path.join(outdir, 'cluster_size_dist.png'), dpi=200)
    fig.savefig(os.path.join(outdir, 'cluster_size_dist.pdf'))

    # permutation studies for the null hypothesis that probes are uniform randomly drawn from the full set of observed
    # data
    m1 = m.loc[:, 'Dura018']
    m2 = m.loc[:, 'GBM018']
    this_cluster = all_results_relevant[(d_max, n_min)]['018']['1']['gene'].values()[0]
    probes = this_cluster['probes']
    n_perm = 10000

    simulated_meds_1 = {}
    simulated_meds_2 = {}

    for n_probe in [4, 6, 8, 10]:
        print "Running permutation study for a cluster with %d probes..." % n_probe
        print "Model 1: randomly drawn probes"
        simulated_meds_1[n_probe] = dmr.cluster_permutation_test(m1, m2, n_probe=n_probe, n_perm=10000)
        print "Model 2: randomly drawn consecutive probes (not spanning chromosomes)"
        simulated_meds_2[n_probe] = dmr.cluster_confined_permutation_test(m1, m2, anno, n_probe=n_probe, n_perm=10000)

    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True)
    for j, n_probe in enumerate([4, 6, 8, 10]):
        ax = axs.flat[j]
        ax.hist(simulated_meds_1[n_probe], 100, normed=True, color='k', label='Random probe selection')
        ax.hist(simulated_meds_2[n_probe], 100, normed=True, color='r', alpha=0.5, label='Consecutive probe selection')
        ax.set_xlabel('Median difference in M (GBM - iNSC)')
        ax.set_ylabel('Density')
        ax.set_title("%d probes in cluster" % n_probe)
    axs.flat[0].legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "median_difference_permutation_test.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "median_difference_permutation_test.pdf"))

    # TODO: compute some actual observed values for different cluster sizes and look at the pvalue similarity
    # we can use the lookup table simulated_meds_2 to get the permutation background

    # if obs_2 < 0:
    #     pval_2 = 2 * np.where(meds_2 < obs_2)[0][0] / float(n_perm)
    # else:
    #     pval_2 = 2 * (1 - np.where(meds_2 > obs_2)[0][0] / float(n_perm))
    #
    #
    # print "Observed median delta: %.2f.\nMann-Whitney p value: %.3e.\nPermutation value: %.3f." % (
    #     obs_1,
    #     this_cluster['pval'],
    #     pval_2
    # )
