from methylation import dmr, process
from load_data import methylation_array
from settings import DATA_DIR
import references
from utils.output import unique_output_dir
import os
import numpy as np
from scipy import stats
import pandas as pd
import multiprocessing as mp
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3, venn2


if __name__ == '__main__':
    outdir = unique_output_dir("rtk1_de_dmr", reuse_empty=True)
    d_max = 200
    n_min = 4
    dm_min = 1.4  # minimum median delta M required to declare a cluster relevant
    alpha = 0.05

    patient_ids = ['018', '019', '031']
    ref_samples = ['insc', 'ensc_h9', 'ensc_fe']
    n_jobs = mp.cpu_count()

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

    clusters = dmr.identify_clusters(anno, n_min=n_min, d_max=d_max, n_jobs=n_jobs)
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
            genes = anno.loc[clusters[chr][cls][cid]].UCSC_RefGene_Name.dropna()
            geneset = reduce(lambda x, y: x.union(y), genes, set())
            attrs['genes'] = geneset

    # 1: What is the joint distribution of methylation / mRNA fold change?
    # Get methylation level and DE fold change for linked genes (pairwise only)

    meth_de_joint = {}
    de_cols = ['genes', 'logFC', 'ensembl', 'direction', 'FDR', 'logCPM']
    meth_cols = ['me_genes', 'me_mediandelta', 'me_medianfc', 'me_fdr']

    for sid in ['018', '019', '031']:
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
                    # print "No matching DE: (%s, %s, %d)" % (chr, cls, cid)
                    continue
                else:
                    print "Found %d matching DE: (%s, %s, %d)" % (de_match.shape[0], chr, cls, cid)
                me_data = np.tile([attrs['median_change'], attrs['median_fc'], attrs['padj']], (de_match.shape[0], 1))
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
                (meth_de_joint[sid][cls], this_match), axis=0
            )
        meth_de_joint[sid]['all'] = pd.concat(meth_de_joint[sid].values(), axis=0)

    # Scatter plots
    for sid in ['018', '019', '031']:
        fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        for i, cls in enumerate(['all', 'tss', 'gene', 'island']):
            ax = axs.flat[i]
            x = meth_de_joint[sid][cls].loc[:, 'logFC']
            y = meth_de_joint[sid][cls].loc[:, 'me_mediandelta']
            # contingency table for Fisher's exact test
            conting = np.array([
                [((x < 0) & (y < 0)).sum(), ((x > 0) & (y < 0)).sum()],
                [((x < 0) & (y > 0)).sum(), ((x > 0) & (y > 0)).sum()],
            ])
            logodds, fisherp = stats.fisher_exact(conting)
            ax.scatter(x, y)
            ax.axhline(0, c=0.4 * np.ones(3))
            ax.axvline(0, c=0.4 * np.ones(3))
            if fisherp < 0.001:
                ax.set_title("%s (Fisher's p=%.2e)" % (cls, fisherp))
            else:
                ax.set_title("%s (Fisher's p=%.3f)" % (cls, fisherp))
        fig.text(0.5, 0.04, 'RNASeq DE logFC', ha='center')
        fig.text(0.04, 0.4, 'EPIC DMR median delta M', va='center', rotation='vertical')
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
            b = "{0:03b}".format(j)
            this_intersection = set(all_genes)
            for k in range(3):
                if b[k] == '1':
                    this_intersection = this_intersection.intersection(
                        meth_de_joint[sids[k]][cls].genes.unique()
                    )
            this_genecount[b] = len(this_intersection)
            this_geneset[b] = list(this_intersection)
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