import multiprocessing as mp
import os
import pickle
import pandas as pd

from rnaseq import loader as rnaseq_loader, general
from methylation import loader as meth_loader, process, dmr
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd
from scripts.hgic_final import two_strategies_combine_de_dmr as tscdd
from scripts.hgic_final import consts
from hgic_consts import NH_ID_TO_PATIENT_ID_MAP
from utils import setops, output, log
from settings import INTERMEDIATE_DIR

"""
Run DE and DMR for the comparison GIC - iAPC
Apply a similar pipeline to that in tscdd to generate a longlist and shortlist of DE/DMR pairs.
"""
outdir = output.unique_output_dir()
logger = log.get_console_logger()

# file location parameters
DMR_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'dmr')
DE_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'de')

de_params = consts.DE_PARAMS
dmr_params = consts.DMR_PARAMS
dmr_params['n_jobs'] = mp.cpu_count()

norm_method_s1 = 'swan'
min_cpm = 1.

pids = [
    '019',
    '031',
    '050',
    '052'
]

# load RNA-Seq data
rna_obj = rnaseq_loader.load_by_patient(
    pids,
    type='cell_culture',
    source='star',
    include_control=False
)


# filter CC to include only HIC and iAPC
rna_obj.filter_samples(rna_obj.meta.index.isin(consts.S1_RNASEQ_SAMPLES_IAPC + consts.S1_RNASEQ_SAMPLES_GIC))
# FIXME: this is a bug in the loader?
rna_obj.batch_id = rna_obj.meta.batch


# load or run DE
the_hash = tsgd.de_results_hash(rna_obj.meta.index, de_params)
filename = 'de_results_paired_comparison.%d.pkl' % the_hash
fn = os.path.join(DE_LOAD_DIR, filename)

if os.path.isfile(fn):
    logger.info("Reading S1 DE results from %s", fn)
    with open(fn, 'rb') as f:
        de_res_full_s1 = pickle.load(f)
else:
    groups_s1 = pd.Series(index=rna_obj.meta.index)
    comparisons_s1 = {}
    for pid in pids:
        groups_s1[rna_obj.meta.index.str.contains(pid) & (rna_obj.meta['type'] == 'GBM')] = "GBM%s" % pid
        groups_s1[groups_s1.index.str.contains(pid) & (rna_obj.meta['type'] == 'iAPC')] = "iAPC%s" % pid

        comparisons_s1[("GBM%s" % pid, "iAPC%s" % pid)] = "GBM%s - iAPC%s" % (pid, pid)

    logger.info("Computing DE results from scratch")
    de_res_full_s1 = tsgd.de_grouped_dispersion(
        rna_obj.data,
        groups_s1,
        comparisons_s1,
        min_cpm=min_cpm,
        return_full=True,
        **de_params
    )
    # rename the keys to simplify
    de_res_full_s1 = dict([(pid, de_res_full_s1[("GBM%s" % pid, "iAPC%s" % pid)]) for pid in pids])

    with open(fn, 'wb') as f:
        pickle.dump(de_res_full_s1, f)

    logger.info("Saved S1 DE results to %s", fn)

# extract only significant DE genes
de_res_s1 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s1.items()])

# generate wide-form lists and save to Excel file
de_by_member = [de_res_s1[pid].index for pid in pids]
venn_set, venn_ct = setops.venn_from_arrays(*de_by_member)

# add null set manually from full DE results
de_genes_all = setops.reduce_union(*venn_set.values())
k_null = ''.join(['0'] * len(pids))
venn_set[k_null] = list(de_res_full_s1[pids[0]].index.difference(de_genes_all))
venn_ct[k_null] = len(venn_set[k_null])

de_data = setops.venn_set_to_wide_dataframe(
    de_res_s1,
    venn_set,
    pids,
    full_data=de_res_full_s1,
    cols_to_include=['logFC', 'FDR'],
    consistency_check_col='logFC',
    consistency_check_method='sign'
)
# add gene symbols back in
general.add_gene_symbols_to_ensembl_data(de_data)
de_data.to_excel(os.path.join(outdir, 'full_de.xlsx'))


# load methylation data
me_obj, anno = tsgd.load_methylation(
    pids,
    type='cell_culture',
    norm_method=norm_method_s1
)
# filter CC to include only iNSC
me_obj.filter_samples(me_obj.meta.index.isin(consts.S1_METHYL_SAMPLES_IAPC + consts.S1_METHYL_SAMPLES_GIC))
# FIXME: this is a bug in the loader?
me_obj.batch_id = me_obj.meta.batch

# generate hash for loading before any relabelling
dmr_hash_dict = dict(dmr_params)
dmr_hash_dict['norm_method'] = norm_method_s1
the_hash = tsgd.dmr_results_hash(me_obj.meta.index, dmr_hash_dict)
filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
fn = os.path.join(DMR_LOAD_DIR, filename)

if os.path.isfile(fn):
    logger.info("Reading S1 DMR results from %s", fn)
    dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
else:
    logger.info("Computing S1 DMR results from scratch.")
    dmr_res_s1 = tsgd.paired_dmr(me_obj.data, me_obj.meta, anno, pids, dmr_params, type1='GBM', type2='iAPC')
    # Save DMR results to disk
    dmr_res_s1.to_pickle(fn, include_annotation=False)
    print "Saved DMR results to %s" % fn

# very adhoc lookup for targets of interest
GOI = [
    'PTGER4',
    'NTRK2',
    'ALDH3B1'
]
dmr_res_s1_sign = dmr_res_s1.results_significant

for g in GOI:
    # any DE results?
    ix = de_data['Gene Symbol'] == g
    assert ix.sum() == 1, "Expected a single DE row to match gene %s" % g
    the_row = de_data[ix].squeeze()
    p_ix = the_row[pids] == 'Y'
    logger.info("Gene %s is DE in %d patients (GIC vs iAPC)", g, p_ix.sum())
    if p_ix.sum() > 0:
        logger.info(', '.join([str(t) for t in p_ix.index[p_ix]]))
        logger.info(', '.join(["%.3f" % the_row["%s_logFC" % t] for t in p_ix.index[p_ix]]))

    # any DM results?
    the_cids = set([k for k in dmr_res_s1.clusters if g in [t[0] for t in dmr_res_s1.clusters[k].genes]])
    logger.info("Gene %s corresponds to %d methylation probe clusters.", g, len(the_cids))
    if len(the_cids) > 0:
        logger.info(', '.join([str(t) for t in the_cids]))
        # get significant results
        the_dm_rel = {k: the_cids.intersection(dmr_res_s1_sign[k]) for k in pids}
        for c in the_cids:
            logger.info("Cluster %d", c)
            the_dm_rel = [k for k in pids if c in dmr_res_s1_sign[k]]
            logger.info("Significant DM in %d patients", len(the_dm_rel))
            if len(the_dm_rel) > 0:
                logger.info("Patients: %s", ', '.join([str(t) for t in the_dm_rel]))
                logger.info("Results: %s", ', '.join([str(dmr_res_s1_sign[t][c]['median_change']) for t in the_dm_rel]))