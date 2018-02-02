import multiprocessing as mp
import os
import pandas as pd
from utils import output, setops
from methylation import dmr, process, loader as methylation_loader
from rnaseq import loader as rnaseq_loader
from analysis import cross_comparison


def load_methylation(pids, norm_method='swan'):
    """
    Load and prepare the Illumina methylation data
    """
    obj = methylation_loader.load_by_patient(pids, norm_method=norm_method)
    me_data = obj.data
    me_data.dropna(inplace=True)
    me_data = process.m_from_beta(me_data)
    anno = methylation_loader.load_illumina_methylationepic_annotation()

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(me_data.index)
    anno = anno.loc[common_probes]
    dmr.add_merged_probe_classes(anno)
    me_data = me_data.loc[common_probes]
    obj.data = me_data

    return obj, anno


def load_rnaseq(pids, ref_names):
    # Load RNA-Seq from STAR
    obj = rnaseq_loader.load_by_patient(pids)

    # load additional references
    ref_obj = rnaseq_loader.load_references(ref_names)
    obj = rnaseq_loader.loader.MultipleBatchLoader([obj, ref_obj])

    return obj


def dmr_results_hash(pids, dmr_params):
    hash_elements = tuple(sorted(pids)) + (
        dmr_params['d_max'],
        dmr_params['n_min'],
        dmr_params['delta_m_min'],
        dmr_params['dmr_test_method'],
    ) + tuple(dmr_params['test_kwargs'].items())
    return hash(hash_elements)


if __name__ == "__main__":

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'GLM'
    }

    dmr_params = {
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'fdr': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
        'n_jobs': mp.cpu_count(),
    }

    pids = ['019', '030', '031', '017', '050', '054']

    subgroups = {
        'RTK I': ['019', '030', '031'],
        'RTK II': ['017', '050', '054'],
    }

    external_ref_names_de = ['GSE61794', 'GSE38993']

    external_refs_de = [
        ('GIBCO', 'NSC'),
        ('H9', 'NSC'),
        ('H1', 'NSC'),
    ]

    external_refs_dm = [
        'GIBCO'
    ]

    # plotting parameters
    cmap = 'RdYlGn_r'

    # file location parameters
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

    # load methylation
    me_obj, anno = load_methylation(pids)
    me_data = me_obj.data

    # Compute DMR
    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    the_hash = dmr_results_hash(pids, dmr_params)
    filename = 'dmr_results.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    loaded = False
    if os.path.isfile(fn):
        dmr_res = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        dmr_res = dmr.compute_cross_dmr(me_data, me_obj.meta, anno, pids, dmr_params)
        # Save DMR results to disk
        dmr_res.to_pickle(fn, include_annotation=False)
        print "Saved DMR results to %s" % fn

    # generate cross-comparison correction
    dmr_classes = ['tss', 'island']

    # extract significant results for tss and island
    dmr_res_sign_all_cls = dmr_res.results_significant_by_class()
    dmr_sign = {}
    for k1 in dmr_res_sign_all_cls:
        for k2 in dmr_res_sign_all_cls[k1]:
            dmr_sign[(k1, k2)] = {}
            for dc in dmr_classes:
                dmr_sign[(k1, k2)].update(dmr_res_sign_all_cls[k1][k2][dc])
            dmr_sign[(k1, k2)] = dmr_sign[(k1, k2)].keys()

    # Intersection across the references gives us a final list that need correcting
    dm_specific_to_all_refs = cross_comparison.compute_cross_comparison_correction(
        dmr_sign,
        pids,
        external_refs=external_refs_dm,
        set_type='pair_only'
    )

    rnaseq_obj = load_rnaseq(pids, external_ref_names_de)

    # discard unneeded samples
    rnaseq_obj.meta = rnaseq_obj.meta.loc[rnaseq_obj.meta.index.str.contains('NSC')]
    rnaseq_obj.data = rnaseq_obj.data.loc[:, rnaseq_obj.meta.index]
