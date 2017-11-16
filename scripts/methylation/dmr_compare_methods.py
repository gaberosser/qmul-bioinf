from methylation import dmr, process
from load_data import methylation_array
import collections
import multiprocessing as mp


if __name__ == "__main__":

    CORE_PARAMS = {
        'core_min_sample_overlap': 3,  # 3 / 4 samples must match
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'fdr': 0.05,
        'reference': 'gibco',
        'ref_name': 'GIBCONSC_P4',
    }

    PARAMS_ARR = {
        'mwu': dict(CORE_PARAMS),
        'permute': dict(CORE_PARAMS),
    }

    # MWU
    PARAMS_ARR['mwu'].update({
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
    })

    # MWU permute
    PARAMS_ARR['permute'].update({
        'dmr_test_method': 'mwu_permute',  # 'mwu', 'mwu_permute'
        'test_kwargs': {'n_max': 1999},
    })

    n_jobs = mp.cpu_count()

    patient_pairs = collections.OrderedDict([
        ('018', (
            ('GBM018_P10', 'GBM018_P12'), ('DURA018_NSC_N4_P4', 'DURA018_NSC_N2_P6')
        )),
        ('019', ('GBM019_P4', 'DURA019_NSC_N8C_P2')),
        ('030', ('GBM030_P5', 'DURA030_NSC_N16B6_P1')),
        ('031', ('GBM031_P4', 'DURA031_NSC_N44B_P2'))
    ])
    patient_pairs['all'] = tuple([
        tuple(reduce(
            lambda x, y: list([x] if isinstance(x, str) else x) + list([y] if isinstance(y, str) else y),
            t,
            []
        )) for t in zip(*patient_pairs.values())
    ])

    ## compute DMR

    anno = methylation_array.load_illumina_methylationepic_annotation()
    b, me_meta = methylation_array.gbm_rtk1_and_paired_nsc(norm_method='swan', ref=CORE_PARAMS['reference'])
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

    clusters = dmr.identify_clusters(anno, n_min=CORE_PARAMS['n_min'], d_max=CORE_PARAMS['d_max'], n_jobs=n_jobs)

    test_results = {}
    test_results_relevant = {}
    test_results_significant = {}

    for k, params in PARAMS_ARR.items():

        print "Run %s." % k

        this_test_results = {}
        this_test_results_relevant = {}
        this_test_results_significant = {}

        for sid, samples in patient_pairs.items():

            this_res = dmr.test_clusters(
                clusters,
                m,
                samples=samples,
                min_median_change=CORE_PARAMS['delta_m_min'],
                n_jobs=n_jobs,
                method=params['dmr_test_method'],
                test_kwargs=params['test_kwargs']
            )

            this_test_results.setdefault(sid, {})
            this_test_results_relevant.setdefault(sid, {})
            this_test_results_significant.setdefault(sid, {})

            ## FIXME: the function  mht_correction now modifies in place; this will error

            this_test_results[sid]['insc_only'] = this_res
            this_test_results_relevant[sid]['insc_only'] = dmr.mht_correction(
                this_test_results[sid]['insc_only'],
                alpha=CORE_PARAMS['fdr']
            )
            this_test_results_significant[sid]['insc_only'] = dmr.filter_dictionary(
                this_test_results_relevant[sid]['insc_only'],
                filt=lambda x: x['rej_h0'],
                n_level=3
            )

            # for insc and ref: use a replacement `samples`
            samples_ref = (samples[0], (CORE_PARAMS['ref_name'],))
            this_res = dmr.test_clusters(
                clusters,
                m,
                samples=samples_ref,
                min_median_change=CORE_PARAMS['delta_m_min'],
                n_jobs=n_jobs,
                method=params['dmr_test_method'],
                test_kwargs=params['test_kwargs']
            )
            this_test_results[sid]['insc_ref'] = this_res
            this_test_results_relevant[sid]['insc_ref'] = dmr.mht_correction(
                this_test_results[sid]['insc_ref'],
                alpha=CORE_PARAMS['fdr']
            )
            this_test_results_significant[sid]['insc_ref'] = dmr.filter_dictionary(
                this_test_results_relevant[sid]['insc_ref'],
                filt=lambda x: x['rej_h0'],
                n_level=3
            )

        # add list of annotated genes to all clusters
        # typ is insc_ref (if present) and insc_only
        for sid in this_test_results:
            for (typ, chr, cls, cid), attrs in dmr.dict_iterator(this_test_results[sid], n_level=4):
                pids = attrs['probes']
                genes = anno.loc[attrs['probes']].UCSC_RefGene_Name.dropna()
                geneset = reduce(lambda x, y: x.union(y), genes, set())
                attrs['genes'] = geneset

        test_results[k] = this_test_results
        test_results_relevant[k] = this_test_results_relevant
        test_results_significant[k] = this_test_results_significant

