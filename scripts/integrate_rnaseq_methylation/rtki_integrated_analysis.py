from methylation import dmr, process, plots
from load_data import methylation_array
from settings import DATA_DIR
from utils.output import unique_output_dir
import collections
import os
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing as mp
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import venn


def construct_contingency(x, y):
    return np.array([
        [((x < 0) & (y < 0)).sum(), ((x > 0) & (y < 0)).sum()],
        [((x < 0) & (y > 0)).sum(), ((x > 0) & (y > 0)).sum()],
    ])


def compute_joint_de_dmr(this_test_results, this_de):
    res = {}

    for sid in this_test_results:
        print sid
        res[sid] = {}

        de_cols = ['genes', 'logFC', 'ensembl', 'direction', 'FDR', 'logCPM']
        meth_cols = ['me_genes', 'chr', 'me_cid', 'me_mediandelta', 'me_median1', 'me_median2', 'me_fdr']
        meth_attrs = ['median_change', 'median1', 'median2']

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
                        [chr, cid] + [attrs[k] for k in meth_attrs] + [attrs['padj']],
                        (de_match.shape[0], 1)
                    )
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
                print "Failed to add data: (%s, %s, %d)" % (chr, cls, cid)
                print repr(exc)
                continue

        res[sid]['all'] = pd.concat(res[sid].values(), axis=0, ignore_index=True)

    return res


if __name__ == '__main__':
    outdir = unique_output_dir("rtk1_de_dmr", reuse_empty=True)
    d_max = 400
    n_min = 6
    dm_min = 1.4  # minimum median delta M required to declare a cluster relevant
    alpha = 0.05
    dmr_method = 'mwu'

    patient_pairs = collections.OrderedDict([
        ('018', (
            ('GBM018_P10', 'GBM018_P12'), ('DURA018_NSC_N4_P4', 'DURA018_NSC_N2_P6')
        )),
        ('019', ('GBM019_P4', 'DURA019_NSC_N8C_P2')),
        ('030', ('GBM030_P5', 'DURA030_NSC_N16B6_P1')),
        ('031', ('GBM031_P4', 'DURA031_NSC_N44B_P2')),
    ])
    n_jobs = mp.cpu_count()
    # n_jobs = 12

    ## load all DE gene lists
    ncol_per_de_block = 6

    indir_de = os.path.join(DATA_DIR, 'rnaseq_de', 'rtk1', 'insc_gibco')
    de = {}

    # blank string corresponds to full set
    for p in patient_pairs.keys() + ['']:
        # label used in DE structure
        if p == '':
            lbl = 'all'
        else:
            lbl = p

        # fn = os.path.join(indir_de, 'gbm-insc-ensc-%s.csv' % p)
        fn = os.path.join(indir_de, 'GBM{0}.vs.iNSC{0}-GBM{0}.vs.refNSC.csv'.format(p))
        this_de_insc_only = pd.read_csv(fn, header=0, index_col=None)
        in_insc = ~this_de_insc_only.iloc[:, 1].isnull()
        in_ensc = ~this_de_insc_only.iloc[:, ncol_per_de_block + 1].isnull()

        # DE genes in iNSC comparison only
        de[(lbl, 'insc_only')] = this_de_insc_only.loc[in_insc & ~in_ensc].iloc[:, :ncol_per_de_block]
        # DE genes in both comparisons
        # here we use the logFC etc from the iNSC comparison, since this is what we're interested in
        de[(lbl, 'insc_and_ref')] = this_de_insc_only.loc[in_insc & in_ensc].iloc[:, :ncol_per_de_block]
        # DE genes in H9 comparison only
        de[(lbl, 'ref_only')] = this_de_insc_only.loc[~in_insc & in_ensc].iloc[:, ncol_per_de_block:]
        # replace column labels to avoid the .1 suffix introduced by R
        de[(lbl, 'ref_only')].columns = this_de_insc_only.columns[:ncol_per_de_block]

    ## compute DMR

    anno = methylation_array.load_illumina_methylationepic_annotation()
    b, me_meta = methylation_array.gbm_rtk1_and_paired_nsc(norm_method='swan')
    b.dropna(inplace=True)
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

    clusters = dmr.identify_clusters(anno, n_min=n_min, d_max=d_max, n_jobs=n_jobs)
    test_results = {}
    test_results_relevant = {}
    test_results_significant = {}
    for sid, samples in patient_pairs.items():
        test_results[sid] = dmr.test_clusters(
            clusters,
            m,
            samples=samples,
            min_median_change=dm_min,
            n_jobs=n_jobs,
            method=dmr_method
        )
        test_results_relevant[sid] = dmr.mht_correction(test_results[sid], alpha=alpha)
        test_results_significant[sid] = dmr.filter_dictionary(
            test_results_relevant[sid],
            filt=lambda x: x['rej_h0'],
            n_level=3
        )

    # add list of annotated genes to all clusters
    for sid in test_results:
        for (chr, cls, cid), attrs in dmr.dict_iterator(test_results[sid], n_level=3):
            pids = attrs['probes']
            genes = anno.loc[attrs['probes']].UCSC_RefGene_Name.dropna()
            geneset = reduce(lambda x, y: x.union(y), genes, set())
            attrs['genes'] = geneset

    # create a table of the numbers of DMRs
    cols = (
        'sample', 'clusters_proposed', 'clusters_relevant', 'clusters_significant',
        'genes_proposed', 'genes_relevant', 'genes_significant'
    )
    table_cluster_numbers = pd.DataFrame(columns=cols)

    def count_genes(res, sid):
        the_genes = reduce(
            lambda x, y: x.union(y),
            [t[1]['genes'] for t in dmr.dict_iterator(res[sid], n_level=3)],
            set()
        )
        return len(the_genes)

    ncl = len(list(
        dmr.dict_iterator(test_results[patient_pairs.keys()[0]], n_level=3)
    ))
    ng = count_genes(test_results, patient_pairs.keys()[0])

    for sid in test_results:
        ncl_re = len(list(
            dmr.dict_iterator(test_results_relevant[sid], n_level=3)
        ))
        ncl_si = len(list(
            dmr.dict_iterator(test_results_significant[sid], n_level=3)
        ))
        ng_re = count_genes(test_results_relevant, sid)
        ng_si = count_genes(test_results_significant, sid)
        this_row = pd.Series({
            'sample': sid,
            'clusters_proposed': ncl,
            'clusters_relevant': ncl_re,
            'clusters_significant': ncl_si,
            'genes_proposed': ng,
            'genes_relevant': ng_re,
            'genes_significant': ng_si
        })
        table_cluster_numbers = table_cluster_numbers.append(this_row, ignore_index=True)

    # 1: What is the joint distribution of methylation / mRNA fold change?
    # Get methylation level and DE fold change for linked genes (pairwise only)

    this_de_insc_only = dict([(sid, de[(sid, 'insc_only')]) for sid in patient_pairs])
    this_de_insc_ref = dict([(sid, de[(sid, 'insc_and_ref')]) for sid in patient_pairs])
    meth_de_joint_insc_only = compute_joint_de_dmr(test_results_significant, this_de_insc_only)
    meth_de_joint_insc_ref = compute_joint_de_dmr(test_results_significant, this_de_insc_ref)

    # Generate table giving the number of overlaps in each patient and cluster class
    # this includes the number of absolute overlaps AND the number of unique overlaps
    the_cols = ['DE genes', 'DMR', 'DMR genes', 'overlaps', 'unique overlaps']
    the_cols += reduce(lambda x, y: x + y, [['%s' % t, '%s_unique' % t] for t in dmr.CLASSES], [])
    de_dmr_matches_insc_only = pd.DataFrame(
        columns=the_cols,
        index=pd.Index(patient_pairs, name='patient'),
    )
    de_dmr_matches_insc_ref = pd.DataFrame.copy(de_dmr_matches_insc_only)

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

    for sid in test_results:
        de_dmr_matches_insc_only.loc[sid] = n_overlap_datum(sid, meth_de_joint_insc_only, de[(sid, 'insc_only')])
        de_dmr_matches_insc_ref.loc[sid] = n_overlap_datum(sid, meth_de_joint_insc_ref, de[(sid, 'insc_and_ref')])


    def scatter_plot_dmr_de(meth_de, fig_filestem):
        for sid in meth_de:
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
    print "*** Genes that match and are DE in GBM vs iNSC AND GBM vs REF ***"
    scatter_plot_dmr_de(meth_de_joint_insc_ref, os.path.join(outdir, "de_vs_dmr_insc_ref"))

    # 2: Venn diagrams: to what extent do the same genes appear in all RTK 1 samples?
    def venn_diagram_and_core_genes(meth_de, text_file, fig_file):
        n_sample = len(meth_de)
        b_str = "{0:0%db}" % n_sample
        all_str = ''.join(['1'] * n_sample)
        all_genes = reduce(
            lambda x, y: x.union(y),
            (set(meth_de[sid]['all'].genes.unique()) for sid in patient_pairs),
            set()
        )

        fig, axs = plt.subplots(ncols=3, figsize=(8, 3.2))
        f = open(text_file, 'wb')

        venn_counts = {}
        venn_sets = {}
        # one Venn diagram per class
        # Venn components are patients
        # we therefore build a structure containing patient-wise gene sets for each class
        for i, cls in enumerate(dmr.CLASSES):
            # reordered dict by sublevel
            tmp = dmr.dict_by_sublevel(meth_de, 2, cls)
            gene_sets = [set(tmp[sid].genes) for sid in patient_pairs]
            venn_res, venn_sets[cls], venn_counts[cls] = venn.venn_diagram(
                *gene_sets,
                set_labels=patient_pairs.keys(),
                ax=axs[i]
            )
            axs[i].set_title(cls)
            print "%s core genes: %s" % (cls, ', '.join(venn_sets[cls][all_str]))
            f.write("%s core genes: %s\n" % (cls, ', '.join(venn_sets[cls][all_str])))
        core_all = set(venn_sets['tss'][all_str]).intersection(venn_sets['gene'][all_str]).intersection(
            venn_sets['island'][all_str])
        print "Core genes shared across all classes: %s" % ', '.join(list(core_all))
        f.write("Core genes shared across all classes: %s\n" % ', '.join(list(core_all)))
        f.close()

        fig.tight_layout()
        fig.savefig("%s.png" % fig_file, dpi=200)
        fig.savefig("%s.pdf" % fig_file)

        return venn_sets


    print "*** GBM vs iNSC and NOT GBM vs REF Venn overlaps ***"
    venn_insc_only = venn_diagram_and_core_genes(
        meth_de_joint_insc_only,
        os.path.join(outdir, "core_genes_de_dmr_insc_only.txt"),
        os.path.join(outdir, "dmr_and_de_overlap_insc_only")
    )
    print "*** GBM vs iNSC and GBM vs REF Venn overlaps ***"
    venn_insc_ref = venn_diagram_and_core_genes(
        meth_de_joint_insc_ref,
        os.path.join(outdir, "core_genes_de_dmr_insc_ref.txt"),
        os.path.join(outdir, "dmr_and_de_overlap_insc_ref")
    )