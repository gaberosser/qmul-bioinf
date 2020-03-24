import os
import pandas as pd

from utils import output, setops, log, genomics
from settings import INTERMEDIATE_DIR
from scripts.hgic_final import two_strategies_combine_de_dmr as tscdd
from scripts.hgic_final import consts
from methylation import dmr

logger = log.get_console_logger()

if __name__ == "__main__":
    outdir = output.unique_output_dir(reuse_empty=True)

    # Set this to True to include reference methylation data
    # This will limit the number of available probes (to 450K)
    include_external_dm_refs = False

    de_params = consts.DE_PARAMS
    dmr_params = consts.DMR_PARAMS

    norm_method_s1 = 'swan'

    pids = consts.PIDS

    if include_external_dm_refs:
        external_ref_names_dm = ['GSE38216']
        external_ref_samples_dm = ['H9 NPC 1', 'H9 NPC 2']
    else:
        external_ref_names_dm = None
        external_ref_samples_dm = None

    if include_external_dm_refs:
        external_refs_dm = [
            ('GIBCO', 'NSC'),
            ('H9', 'NSC'),
        ]
    else:
        external_refs_dm = [
            ('GIBCO', 'NSC'),
        ]

    external_refs_dm_labels = [t[0] for t in external_refs_dm]

    # plotting parameters
    cmap = 'RdYlGn_r'

    # file location parameters
    DMR_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'dmr')

    if not os.path.isdir(DMR_LOAD_DIR):
        raise AttributeError("To run this script, we require pre-computed DMR results residing in %s" % DMR_LOAD_DIR)

    #########################################################################
    ### STRATEGY 1: No references, just compare GBM-iNSC for each patient ###
    #########################################################################

    # some quantities relating to set membership

    # load methylation
    # For S1, we just need the paired comparison. This is important - any other samples lead to a change in the final
    # probe list (any NA rows are dropped from ALL samples).
    me_obj, anno = tscdd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
    me_data = me_obj.data
    me_meta = me_obj.meta
    me_meta.insert(0, 'patient_id', me_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    the_hash = tscdd.dmr_results_hash(me_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise AttributeError("Unable to load pre-computed DMR results, expected at %s" % fn)

    # extract results
    dmr_res_full_s1 = dmr_res_s1.results
    dmr_res_sign_s1 = dmr_res_s1.results_significant

    ##################
    ### STRATEGY 2 ###
    ##################

    norm_method_s2 = 'pbc'

    # load methylation with external references
    me_obj_with_ref, anno = tscdd.load_methylation(
        pids,
        ref_names=external_ref_names_dm,
        ref_name_filter=external_ref_samples_dm,
        norm_method=norm_method_s2,
        patient_samples=consts.ALL_METHYL_SAMPLES  # technically contains H9 too, but these will be ignored
    )
    me_data_with_ref = me_obj_with_ref.data

    # Load DMR cross-comparison

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise this is written after computing the results

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s2

    the_hash = tscdd.dmr_results_hash(me_obj_with_ref.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_cross_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S2 DMR results from %s", fn)
        dmr_res_s2 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise AttributeError("Unable to load pre-computed DMR results, expected at %s" % fn)

    dmr_res_full_s2 = {}
    dmr_res_sign_s2 = {}

    for k1 in dmr_res_s2.keys():
        dmr_res_full_s2.setdefault('syngeneic', {})[k1] = dmr_res_s2[k1][k1]
        dmr_res_sign_s2.setdefault('syngeneic', {})[k1] = dmr_res_s2[k1][k1]
        for k2 in external_refs_dm_labels:
            dmr_res_full_s2.setdefault(k2, {})[k1] = dmr_res_s2[k1][k2]
            dmr_res_sign_s2.setdefault(k2, {})[k1] = dmr_res_s2[k1][k2]


    # export results with genomic region
    # it turns out that these are actually THE SAME for S1 and S2 (in the case tested here - but possibly NOT for other
    # cases!)
    # Blow for blow, DMR results differ due to the different norming routine: SWAN for S1 vs PBC for S2.
    clusters = dmr_res_s1.clusters
    coords = anno.MAPINFO
    cluster_mapping = {}
    for cid in clusters:
        cl = clusters[cid]
        pids = cl.pids
        ch = anno.loc[cl.pids, 'CHR'].unique()
        assert len(ch) == 1, "Cluster {} maps to multiple chromosomes.".format(cid)
        ch = ch[0]
        pl = coords.reindex(cl.pids)
        x0 = pl.min()
        x1 = pl.max()
        cluster_mapping[cid] = {'chromosome': ch, 'start': x0, 'end': x1}

    cluster_mapping = pd.DataFrame(cluster_mapping).transpose()

    cluster_mapping.to_csv(os.path.join(outdir, 'cluster_mapping_grch37.csv'))

    # bed format
    region_data = {}
    for cid, row in cluster_mapping.iterrows():
        region_data[str(cid)] = [row['chromosome'], row['start'], row['end'], '.']
    genomics.write_bed_file(region_data, os.path.join(outdir, 'cluster_mapping_grch37.bed'))
