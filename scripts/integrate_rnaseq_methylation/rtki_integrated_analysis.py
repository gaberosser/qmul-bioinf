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
from utils import setops
import json


class TestResultEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, set):
            return tuple(o)
        elif isinstance(o, bool) or isinstance(o, np.bool_):
            return int(o)
        return super(TestResultEncoder, self).default(o)



def construct_contingency(x, y):
    return np.array([
        [((x < 0) & (y < 0)).sum(), ((x > 0) & (y < 0)).sum()],
        [((x < 0) & (y > 0)).sum(), ((x > 0) & (y > 0)).sum()],
    ])


def compute_joint_de_dmr(dmr_results, de_results):
    res = {}

    for sid in dmr_results:
        print sid
        res[sid] = {}

        de_cols = ['genes', 'logFC', 'ensembl', 'direction', 'FDR', 'logCPM']
        meth_cols = ['me_genes', 'chr', 'me_cid', 'me_mediandelta', 'me_median1', 'me_median2', 'me_fdr']
        meth_attrs = ['median_change', 'median1', 'median2']

        for (chr, cls, cid), attrs in dmr.dict_iterator(dmr_results[sid], n_level=3):
            res[sid].setdefault(cls, pd.DataFrame(columns=de_cols + meth_cols))

            if len(attrs['genes']) == 0:
                continue

            try:
                # matching entry in DE (by gene name)
                de_match = de_results[sid].loc[de_results[sid].loc[:, 'genes'].isin(attrs['genes'])]

                if de_match.shape[0] > 0:
                    # form the DMR data block by repeating the same row
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

        # combine all methylation cluster classes
        res[sid]['all'] = pd.concat(res[sid].values(), axis=0, ignore_index=True)

    return res


if __name__ == '__main__':

    PARAMS = {
        'core_min_sample_overlap': 3,  # 3 / 4 samples must match
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'fdr': 0.05,
        # 'dmr_test_method': 'mwu_permute',  # 'mwu', 'mwu_permute'
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        # 'test_kwargs': {'n_max': 1999},
        'test_kwargs': {},
    }

    outdir = unique_output_dir(
        "rtk1_de_dmr.%s.gibco_reference" % PARAMS['dmr_test_method'],
        reuse_empty=True
    )
    indir_de = os.path.join(DATA_DIR, 'rnaseq_de', 'rtk1', 'insc_gibco')  # can also use H9

    # write params to a file
    with open(os.path.join(outdir, 'parameters.json'), 'wb') as f:
        json.dump(PARAMS, f)

    patient_pairs = collections.OrderedDict([
        ('018', (
            ('GBM018_P10', 'GBM018_P12'), ('DURA018_NSC_N4_P4', 'DURA018_NSC_N2_P6')
        )),
        ('019', ('GBM019_P4', 'DURA019_NSC_N8C_P2')),
        ('030', ('GBM030_P5', 'DURA030_NSC_N16B6_P1')),
        ('031', ('GBM031_P4', 'DURA031_NSC_N44B_P2')),
        ('all', (
            ('GBM018_P10', 'GBM018_P12', 'GBM019_P4', 'GBM030_P5', 'GBM031_P4'),
            ('DURA018_NSC_N4_P4', 'DURA018_NSC_N2_P6', 'DURA019_NSC_N8C_P2', 'DURA030_NSC_N16B6_P1', 'DURA031_NSC_N44B_P2')
        )
        ),
    ])
    ref_name = 'GIBCONSC_P4'
    comparisons = ['gbm_insc', 'gbm_ref']
    n_jobs = mp.cpu_count()

    ## load all DE gene lists
    ncol_per_de_block = 6

    de = {}

    # blank string corresponds to full set
    for lbl in patient_pairs.keys():
        # label used in DE structure
        if lbl == 'all':
            # blank string for file naming purposes
            p = ''
        else:
            p = lbl

        # fn = os.path.join(indir_de, 'gbm-insc-ensc-%s.csv' % p)
        fn = os.path.join(indir_de, 'GBM{0}.vs.iNSC{0}-GBM{0}.vs.refNSC.csv'.format(p))
        this_de_insc_only = pd.read_csv(fn, header=0, index_col=None)
        in_insc = ~this_de_insc_only.iloc[:, 1].isnull()
        in_ensc = ~this_de_insc_only.iloc[:, ncol_per_de_block + 1].isnull()

        # FIXME: we could use a nested dict here, like we do for DMR?

        # DE genes in iNSC comparison only
        de[(lbl, 'gbm_insc')] = this_de_insc_only.loc[in_insc & ~in_ensc].iloc[:, :ncol_per_de_block]
        # DE genes in both comparisons
        # here we use the logFC etc from the iNSC comparison, since this is what we're interested in
        de[(lbl, 'gbm_insc_and_ref')] = this_de_insc_only.loc[in_insc & in_ensc].iloc[:, :ncol_per_de_block]

        # DE genes in H9 comparison only - not required
        # de[(lbl, 'ref_only')] = this_de_insc_only.loc[~in_insc & in_ensc].iloc[:, ncol_per_de_block:]
        # replace column labels to avoid the .1 suffix introduced by R
        # de[(lbl, 'ref_only')].columns = this_de_insc_only.columns[:ncol_per_de_block]

    ## compute DMR

    anno = methylation_array.load_illumina_methylationepic_annotation()
    b, me_meta = methylation_array.gbm_rtk1_and_paired_nsc(norm_method='swan', ref='gibco')
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

    clusters = dmr.identify_clusters(anno, n_min=PARAMS['n_min'], d_max=PARAMS['d_max'], n_jobs=n_jobs)

    test_results = {}
    test_results_relevant = {}
    test_results_significant = {}
    for sid, samples in patient_pairs.items():

        test_results.setdefault(sid, {})
        test_results_relevant.setdefault(sid, {})
        test_results_significant.setdefault(sid, {})

        test_results[sid]['gbm_insc'] = dmr.test_clusters(
            clusters,
            m,
            samples=samples,
            min_median_change=PARAMS['delta_m_min'],
            n_jobs=n_jobs,
            method=PARAMS['dmr_test_method'],
            test_kwargs=PARAMS['test_kwargs']
        )

        # for insc and ref: use a replacement `samples`
        samples_ref = (samples[0], (ref_name,))
        test_results[sid]['gbm_ref'] = dmr.test_clusters(
            clusters,
            m,
            samples=samples_ref,
            min_median_change=PARAMS['delta_m_min'],
            n_jobs=n_jobs,
            method=PARAMS['dmr_test_method'],
            test_kwargs=PARAMS['test_kwargs']
        )

        for typ in comparisons:
            test_results_relevant[sid][typ] = dmr.mht_correction(
                test_results[sid][typ],
                alpha=PARAMS['fdr']
            )
            test_results_significant[sid][typ] = dmr.filter_dictionary(
                test_results_relevant[sid][typ],
                filt=lambda x: x['rej_h0'],
                n_level=3
            )

    # add list of annotated genes to all clusters
    for sid in test_results:
        for (typ, chr, cls, cid), attrs in dmr.dict_iterator(test_results[sid], n_level=4):
            pids = attrs['probes']
            genes = anno.loc[attrs['probes']].UCSC_RefGene_Name.dropna()
            geneset = reduce(lambda x, y: x.union(y), genes, set())
            attrs['genes'] = geneset

    # write full set of results to disk (json format)
    # NB this converts gene sets to lists and cluster ID to a string
    with open(os.path.join(outdir, "dmr_results.json"), 'wb') as f:
        json.dump(test_results, f, cls=TestResultEncoder)

    # create a table of the numbers of DMRs
    cols = (
        'sample', 'clusters_proposed', 'clusters_relevant', 'clusters_significant',
        'genes_proposed', 'genes_relevant', 'genes_significant'
    )
    # table_cluster_numbers = pd.DataFrame(columns=cols)
    table_cluster_numbers = {}

    def count_genes(res):
        the_genes = reduce(
            lambda x, y: x.union(y),
            [t[1]['genes'] for t in dmr.dict_iterator(res, n_level=3)],
            set()
        )
        return len(the_genes)

    # the number of proposed clusters and proposed genes is constant across all results, so just use the first
    k0 = patient_pairs.keys()[0]
    k1 = test_results.values()[0].keys()[0]
    ncl = len(list(
        dmr.dict_iterator(test_results[k0][k1], n_level=3)
    ))
    ng = count_genes(test_results[k0][k1])

    for sid in test_results:
        for typ in comparisons:
            table_cluster_numbers.setdefault(typ, pd.DataFrame(columns=cols))
            ncl_re = len(list(
                dmr.dict_iterator(test_results_relevant[sid][typ], n_level=3)
            ))
            ncl_si = len(list(
                dmr.dict_iterator(test_results_significant[sid][typ], n_level=3)
            ))
            ng_re = count_genes(test_results_relevant[sid][typ])
            ng_si = count_genes(test_results_significant[sid][typ])
            this_row = pd.Series({
                'sample': sid,
                'clusters_proposed': ncl,
                'clusters_relevant': ncl_re,
                'clusters_significant': ncl_si,
                'genes_proposed': ng,
                'genes_relevant': ng_re,
                'genes_significant': ng_si
            })
            table_cluster_numbers[typ] = table_cluster_numbers[typ].append(this_row, ignore_index=True)

    for typ in comparisons:
        table_cluster_numbers[typ].to_csv(os.path.join(outdir, "cluster_numbers.%s.csv" % typ))

    # split into (insc not ref) and (insc and ref) groups, like DE

    this_dmr_insc_only = {}
    this_dmr_insc_ref = {}

    # for convenience:
    insc_results = dmr.dict_by_sublevel(test_results_significant, 2, 'gbm_insc')
    ref_results = dmr.dict_by_sublevel(test_results_significant, 2, 'gbm_ref')

    for sid in test_results_significant:

        this_dmr_insc_only.setdefault(sid, {})
        this_dmr_insc_ref.setdefault(sid, {})

        for chr in insc_results[sid]:
            this_dmr_insc_only[sid].setdefault(chr, {})
            this_dmr_insc_ref[sid].setdefault(chr, {})
            for cls in dmr.CLASSES:
                this_dmr_insc_only[sid][chr].setdefault(cls, {})
                this_dmr_insc_ref[sid][chr].setdefault(cls, {})

                in_insc = set(insc_results[sid][chr][cls].keys())
                in_ref = set(ref_results[sid][chr][cls].keys())
                in_insc_only = in_insc.difference(in_ref)
                in_insc_and_ref = in_insc.intersection(in_ref)

                # use the iNSC results in both cases, as that is what we are interested in
                this_dmr_insc_only[sid][chr][cls] = dict([
                    (cid, insc_results[sid][chr][cls][cid]) for cid in in_insc_only
                ])
                this_dmr_insc_ref[sid][chr][cls] = dict([
                    (cid, insc_results[sid][chr][cls][cid]) for cid in in_insc_and_ref
                ])

    # counts of clusters and genes in DMR
    cols += (
        'clusters_significant_gbm_insc',
        'genes_significant_gbm_insc',
        'clusters_significant_gbm_ref',
        'genes_significant_gbm_ref',
    )
    table_clusters_ref_insc = pd.DataFrame(columns=cols)

    for sid in insc_results:
        old_row = table_cluster_numbers['gbm_insc'].loc[
            table_cluster_numbers['gbm_insc'].loc[:, 'sample'] == sid
        ]
        this_row = pd.Series(old_row.iloc[0], copy=True)
        this_row.loc['clusters_significant_gbm_insc'] = len(list(
            dmr.dict_iterator(this_dmr_insc_only[sid], n_level=3)
        ))
        this_row.loc['clusters_significant_gbm_ref'] = len(list(
            dmr.dict_iterator(this_dmr_insc_ref[sid], n_level=3)
        ))
        this_row.loc['genes_significant_gbm_insc'] = count_genes(this_dmr_insc_only[sid])
        this_row.loc['genes_significant_gbm_ref'] = count_genes(this_dmr_insc_ref[sid])
    table_clusters_ref_insc.to_csv(os.path.join(outdir, "cluster_numbers_gbm_ref.csv"))

    """
    Before we start working with Venn plots or anything, I want to check the overlap between DMR on individuals and
    DMR on the pooled samples.
    We use the chr, class and cluster ID as unique IDs here.
    """
    print "(iNSC vs GBM) AND NOT (iNSC vs reference)"
    for cls in dmr.CLASSES:
        a = dmr.dict_by_sublevel(
            this_dmr_insc_only['all'],
            2,
            cls
        )
        ga = dmr.dict_iterator(a, n_level=2)
        a_keys = set([tuple(t[0]) for t in ga])

        b_keys = set()
        for sid in [t for t in patient_pairs if t != 'all']:
            gb = dmr.dict_iterator(
                dmr.dict_by_sublevel(this_dmr_insc_only[sid], 2, cls),
                n_level=2
            )
            b_keys.update(set([tuple(t[0]) for t in gb]))

        print "***"
        print "Individual comparisons: %d DMRs in class %s" % (len(a_keys), cls)
        print "Lumped comparisons: %d DMRs in class %s" % (len(b_keys), cls)
        print "Intersection of those: %d DMRs in class %s" % (len(b_keys.intersection(a_keys)), cls)
    print "***"

    print "(iNSC vs GBM) AND (iNSC vs reference)"
    for cls in dmr.CLASSES:
        a = dmr.dict_by_sublevel(
            this_dmr_insc_ref['all'],
            2,
            cls
        )
        ga = dmr.dict_iterator(a, n_level=2)
        a_keys = set([tuple(t[0]) for t in ga])

        b_keys = set()
        for sid in [t for t in patient_pairs if t != 'all']:
            gb = dmr.dict_iterator(
                dmr.dict_by_sublevel(this_dmr_insc_ref[sid], 2, cls),
                n_level=2
            )
            b_keys.update(set([tuple(t[0]) for t in gb]))

        print "***"
        print "Individual comparisons: %d DMRs in class %s" % (len(a_keys), cls)
        print "Lumped comparisons: %d DMRs in class %s" % (len(b_keys), cls)
        print "Intersection of those: %d DMRs in class %s" % (len(b_keys.intersection(a_keys)), cls)
    print "***"

    # 1: What is the joint distribution of methylation / mRNA fold change?
    # Get methylation level and DE fold change for linked genes (pairwise only)

    this_de_insc_only = dict([(sid, de[(sid, 'gbm_insc')]) for sid in patient_pairs])
    this_de_insc_ref = dict([(sid, de[(sid, 'gbm_insc_and_ref')]) for sid in patient_pairs])

    meth_de_joint_insc_only = compute_joint_de_dmr(this_dmr_insc_only, this_de_insc_only)
    meth_de_joint_insc_ref = compute_joint_de_dmr(this_dmr_insc_ref, this_de_insc_ref)

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
        n_dmr = len(list(dmr.dict_iterator(insc_results[sid], n_level=3)))
        n_dmr_genes = len(
            reduce(
                lambda x, y: x.union(y), [t[1]['genes'] for t in dmr.dict_iterator(insc_results[sid], n_level=3)], set()
            )
        )
        n_overlaps = meth_de['all'].shape[0]
        n_overlaps_unique = meth_de['all'].me_genes.unique().shape[0]
        this_datum = [n_de, n_dmr, n_dmr_genes, n_overlaps, n_overlaps_unique]
        for cls in dmr.CLASSES:
            this_datum.append(meth_de[cls].shape[0])
            this_datum.append(meth_de[cls].me_genes.unique().shape[0])
        return this_datum

    for sid in test_results:
        de_dmr_matches_insc_only.loc[sid] = n_overlap_datum(sid, meth_de_joint_insc_only[sid], de[(sid, 'gbm_insc')])
        de_dmr_matches_insc_ref.loc[sid] = n_overlap_datum(sid, meth_de_joint_insc_ref[sid], de[(sid, 'gbm_insc_and_ref')])

    def scatter_plot_dmr_de(meth_de, fig_filestem, fig_titlestem=''):
        for sid in meth_de:
            fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, num="%s %s" % (fig_titlestem, sid))
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
    scatter_plot_dmr_de(meth_de_joint_insc_only, os.path.join(outdir, "de_vs_dmr_insc_only"), "GBM v iNSC AND NOT GBM vs REF")
    print "*** Genes that match and are DE in GBM vs iNSC AND GBM vs REF ***"
    scatter_plot_dmr_de(meth_de_joint_insc_ref, os.path.join(outdir, "de_vs_dmr_insc_ref"), "GBM v iNSC AND GBM vs REF")

    # 2: Venn diagrams: to what extent do the same genes appear in all RTK 1 samples?
    def venn_diagram_and_core_genes(meth_de, text_file, fig_file, min_overlap=4, fig_title=None):
        n_sample = len(meth_de)
        bns = list(setops.binary_combinations_sum_gte(n_sample, min_overlap))

        fig_kws = {'figsize': (8, 3.2)}
        if fig_title is not None:
            fig_kws['num'] = fig_title
        fig, axs = plt.subplots(ncols=3, **fig_kws)

        f = open(text_file, 'wb')
        venn_counts = {}
        venn_sets = {}
        core_sets = {}
        # one Venn diagram per class
        # Venn components are patients
        # we therefore build a structure containing patient-wise gene sets for each class
        for i, cls in enumerate(dmr.CLASSES):
            # reordered dict by sublevel
            tmp = dmr.dict_by_sublevel(meth_de, 2, cls)
            # build labels and gene lists to control the order (otherwise we're relying on dict iteration order)
            set_labels = []
            gene_sets = []
            for k, t in tmp.items():
                gene_sets.append(set(t.genes))
                set_labels.append(k)
            venn_res, venn_sets[cls], venn_counts[cls] = venn.venn_diagram(
                *gene_sets,
                set_labels=set_labels,
                ax=axs[i]
            )
            axs[i].set_title(cls)

            core_genes = set()
            for bn in bns:
                core_genes = core_genes.union(venn_sets[cls][bn])
            core_sets[cls] = core_genes
            core_genes = sorted(core_genes)
            print "%s core genes [%d]: %s" % (cls, len(core_genes), ', '.join(core_genes))
            f.write("%s core genes [%d]: %s\n" % (cls, len(core_genes), ', '.join(core_genes)))

        # intersection of all
        core_all = sorted(reduce(set.intersection, core_sets.values()))
        print "Core genes shared across all classes [%d]: %s" % (len(core_all), ', '.join(core_all))
        f.write("Core genes shared across all classes [%d]: %s\n" % (len(core_all), ', '.join(core_all)))

        # union of all
        core_union = sorted(reduce(set.union, core_sets.values()))
        print "Core genes in >=1 class [%d]: %s" % (len(core_union), ', '.join(core_union))
        f.write("Core genes in >= 1 class [%d]: %s\n" % (len(core_union), ', '.join(core_union)))
        f.close()

        fig.tight_layout()
        fig.savefig("%s.png" % fig_file, dpi=200)
        fig.savefig("%s.pdf" % fig_file)

        return venn_sets

    # since a 5-way Venn is not supported (and is horrible to look at even if you can draw it...), remove the
    # 'all' results here and report them separately

    print "*** GBM vs iNSC and NOT GBM vs REF Venn overlaps - individual patient ***"
    this_de_dmr = dict(meth_de_joint_insc_only)
    this_de_dmr.pop('all')
    venn_insc_only = venn_diagram_and_core_genes(
        this_de_dmr,
        os.path.join(outdir, "core_genes_de_dmr_insc_only.by_patient.txt"),
        os.path.join(outdir, "dmr_and_de_overlap_insc_only"),
        min_overlap=PARAMS['core_min_sample_overlap'],
        fig_title='Individual comparison, iNSC and not ref'
    )
    print "*** GBM vs iNSC and GBM vs REF Venn overlaps - individual patient ***"
    this_de_dmr = dict(meth_de_joint_insc_ref)
    this_de_dmr.pop('all')
    venn_insc_ref = venn_diagram_and_core_genes(
        this_de_dmr,
        os.path.join(outdir, "core_genes_de_dmr_insc_ref.by_patient.txt"),
        os.path.join(outdir, "dmr_and_de_overlap_insc_ref"),
        min_overlap=PARAMS['core_min_sample_overlap'],
        fig_title='Individual comparison, iNSC and ref'
    )
    print "*** GBM vs iNSC and NOT GBM vs REF Venn overlaps - whole cohort ***"
    core_sets = {}
    with open(os.path.join(outdir, "core_genes_de_dmr_insc_only.by_cohort.txt"), 'wb') as f:
        for i, cls in enumerate(dmr.CLASSES):
            core_sets[cls] = sorted(set(meth_de_joint_insc_only['all'][cls].genes))
            print "%s core genes [%d]: %s" % (cls, len(core_sets[cls]), ', '.join(core_sets[cls]))
            f.write("%s core genes [%d]: %s\n" % (cls, len(core_sets[cls]), ', '.join(core_sets[cls])))

        # intersection of all
        core_all = sorted(reduce(lambda x, y: set(x).intersection(y), core_sets.values()))
        print "Core genes shared across all classes [%d]: %s" % (len(core_all), ', '.join(core_all))
        f.write("Core genes shared across all classes [%d]: %s\n" % (len(core_all), ', '.join(core_all)))

        # union of all
        core_union = sorted(reduce(lambda x, y: set(x).union(y), core_sets.values()))
        print "Core genes in >=1 class [%d]: %s" % (len(core_union), ', '.join(core_union))
        f.write("Core genes in >= 1 class [%d]: %s\n" % (len(core_union), ', '.join(core_union)))

    print "*** GBM vs iNSC and GBM vs REF Venn overlaps - whole cohort ***"
    core_sets = {}
    with open(os.path.join(outdir, "core_genes_de_dmr_insc_ref.by_cohort.txt"), 'wb') as f:
        for i, cls in enumerate(dmr.CLASSES):
            core_sets[cls] = sorted(set(meth_de_joint_insc_ref['all'][cls].genes))
            print "%s core genes [%d]: %s" % (cls, len(core_sets[cls]), ', '.join(core_sets[cls]))
            f.write("%s core genes [%d]: %s\n" % (cls, len(core_sets[cls]), ', '.join(core_sets[cls])))

        # intersection of all
        core_all = sorted(reduce(lambda x, y: set(x).intersection(y), core_sets.values()))
        print "Core genes shared across all classes [%d]: %s" % (len(core_all), ', '.join(core_all))
        f.write("Core genes shared across all classes [%d]: %s\n" % (len(core_all), ', '.join(core_all)))

        # union of all
        core_union = sorted(reduce(lambda x, y: set(x).union(y), core_sets.values()))
        print "Core genes in >=1 class [%d]: %s" % (len(core_union), ', '.join(core_union))
        f.write("Core genes in >= 1 class [%d]: %s\n" % (len(core_union), ', '.join(core_union)))

    # represent these results on a plot
    core_sets_insc_only = {}
    core_sets_insc_ref = {}

    bns = list(setops.binary_combinations_sum_gte(len(patient_pairs) - 1, len(patient_pairs) - 2))
    for i, cls in enumerate(dmr.CLASSES):
        core_genes = set()
        for bn in bns:
            core_genes = core_genes.union(venn_insc_only[cls][bn])
        core_sets_insc_only[cls] = core_genes

        core_genes = set()
        for bn in bns:
            core_genes = core_genes.union(venn_insc_ref[cls][bn])
        core_sets_insc_ref[cls] = core_genes


    # check the direction of change in each sample - these must be consistent for the gene to be considered core
    de_values = {}
    dmr_values = {}
    for i, cls in enumerate(dmr.CLASSES):
        # meth_de = dmr.dict_by_sublevel(meth_de_joint_insc_ref, 2, cls)
        meth_de = dmr.dict_by_sublevel(meth_de_joint_insc_only, 2, cls)

        per_sample_de = dict([(k, set(x.genes)) for k, x in meth_de.items()])
        all_genes = reduce(set.union, per_sample_de.values())
        de_values[cls] = pd.DataFrame(index=sorted(all_genes), columns=patient_pairs.keys())

        # Step 1: check each gene for DE compatibility
        for g in all_genes:
            this_row = []
            for sid, this_dat in meth_de.items():
                this_genes = this_dat.genes
                if g in this_genes.values:
                    de_values[cls].loc[g, sid] = this_dat.loc[this_genes == g, 'logFC'].values[0]

        # find those genes shared by >1 sample
        role_call = (~de_values[cls].isnull()).sum(axis=1)
        shared_de = de_values[cls].loc[
            role_call > 1
        ].fillna(0)
        obs = np.sign(shared_de).astype(int).sum(axis=1).abs()
        expctd = shared_de.astype(bool).astype(int).sum(axis=1)
        if not (obs == expctd).all():
            to_remove = obs.loc[obs != expctd].index
            print "Not all DE genes had the same direction change. Removing %s." % ', '.join(to_remove)
            ## TODO: remove
            raise NotImplementedError

        # Step 2: check each matching DMR for compatibility
        per_sample_dmr = dict(
            [(k, set([tuple(t) for t in x.loc[:, ['chr', 'me_cid']].values])) for k, x in meth_de.items()]
        )
        all_dmr = reduce(set.union, per_sample_dmr.values())

        dmr_values[cls] = pd.DataFrame(index=sorted(all_dmr), columns=patient_pairs.keys())
        for chr, cid in all_dmr:
            this_row = []
            for sid, this_dat in meth_de.items():
                idx = ((this_dat.loc[:, 'chr'] == chr) & (this_dat.loc[:, 'me_cid'] == cid))
                if idx.any():
                    dmr_values[cls].loc[(chr, cid), sid] = this_dat.loc[idx, 'me_mediandelta'].values[0]
        role_call = (~dmr_values[cls].isnull()).sum(axis=1)
        shared_dmr = dmr_values[cls].loc[
            role_call > 1
        ].fillna(0)
        obs = np.sign(shared_dmr).astype(int).sum(axis=1).abs()
        expctd = shared_dmr.astype(bool).astype(int).sum(axis=1)
        if not (obs == expctd).all():
            to_remove = obs.loc[obs != expctd].index
            print "Not all DE genes had the same direction change. Removing %s." % ', '.join(to_remove)
            ## TODO: remove
            raise NotImplementedError


    # OK!




