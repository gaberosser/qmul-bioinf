import multiprocessing as mp
import os
import pandas as pd
import numpy as np
import pickle
from utils import output, setops, excel, log
from methylation import dmr, process, loader as methylation_loader, annotation_gene_to_ensembl
from rnaseq import loader as rnaseq_loader, differential_expression, general, filter
from integrator import rnaseq_methylationarray
from analysis import cross_comparison
from load_data import loader
from plotting import venn
from matplotlib import pyplot as plt, text, patches
import seaborn as sns


logger = log.get_console_logger()


def load_methylation(pids, ref_names=None, norm_method='swan', ref_name_filter=None):
    """
    Load and prepare the Illumina methylation data
    """
    # patient data
    obj = methylation_loader.load_by_patient(pids, norm_method=norm_method)
    anno = methylation_loader.load_illumina_methylationepic_annotation()

    # reference data
    if ref_names is not None:
        ref_obj = methylation_loader.load_reference(ref_names, norm_method=norm_method)
        if ref_name_filter is not None:
            ref_obj.filter_by_sample_name(ref_name_filter, exact=True)
        obj = loader.MultipleBatchLoader([obj, ref_obj])

    me_data = obj.data.dropna()
    me_data = process.m_from_beta(me_data)

    # reduce anno and data down to common probes
    common_probes = anno.index.intersection(me_data.index)

    anno = anno.loc[common_probes]
    # dmr.add_merged_probe_classes(anno)
    me_data = me_data.loc[common_probes]
    obj.data = me_data

    return obj, anno


def load_rnaseq(pids, ref_names, ref_name_filter='NSC', discard_filter='IPSC', strandedness=None):
    """
    :param strandedness: Iterable of same length as ref_names giving the strandedness of each ref
    """
    if strandedness is None:
        strandedness = ['u'] * len(ref_names)
    else:
        if len(strandedness) != len(ref_names):
            raise ValueError("Supplied strandedness must be a list of the same length as the ref_names.")

    # Load RNA-Seq from STAR
    obj = rnaseq_loader.load_by_patient(pids)

    # load additional references
    ref_objs = []
    for rn, strnd in zip(ref_names, strandedness):
        ref_obj = rnaseq_loader.load_references(rn, strandedness=strnd)
        if ref_name_filter is not None:
            # only keep relevant references
            ref_obj.meta = ref_obj.meta.loc[ref_obj.meta.index.str.contains(ref_name_filter)]
            ref_obj.data = ref_obj.data.loc[:, ref_obj.meta.index]
        ref_objs.append(ref_obj)
    obj = loader.MultipleBatchLoader([obj] + ref_objs)

    if discard_filter is not None:
        if not hasattr(discard_filter, '__iter__'):
            discard_filter = [discard_filter]
        for d in discard_filter:
            obj.meta = obj.meta.loc[~obj.meta.index.str.contains(d)]
            obj.data = obj.data.loc[:, obj.meta.index]

    obj.batch_id = obj.batch_id.loc[obj.meta.index]

    return obj


def paired_dmr(me_data, me_meta, anno, pids, dmr_params):
    """
    Compute DMRs for paired GBM-iNSC comparisons (defined in that order) for all patients
    :param me_data:
    :param me_meta:
    :param anno:
    :param pids:
    :param dmr_params:
    :return:
    """
    dmr_res_obj = dmr.DmrResults(anno=anno)
    dmr_res_obj.identify_clusters(**dmr_params)
    dmr_res = {}

    for pid in pids:
        the_idx1 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == 'GBM')
        the_idx2 = me_meta.index.str.contains(pid) & (me_meta.loc[:, 'type'] == 'iNSC')
        # control comparison order
        the_samples = [
            me_meta.index[the_idx1],
            me_meta.index[the_idx2],
        ]

        # the_idx = the_idx1 | the_idx2
        # the_groups = me_meta.loc[the_idx, 'type'].values
        # the_samples = me_meta.index[the_idx].groupby(the_groups).values()

        the_obj = dmr_res_obj.copy()
        the_obj.test_clusters(me_data,
                              samples=the_samples,
                              n_jobs=dmr_params['n_jobs'],
                              min_median_change=dmr_params['delta_m_min'],
                              method=dmr_params['dmr_test_method'],
                              alpha=dmr_params['alpha'],
                              **dmr_params['test_kwargs']
                              )
        dmr_res[pid] = the_obj

    return dmr.DmrResultCollection(**dmr_res)


def dmr_results_hash(pids, dmr_params):
    hash_elements = tuple(sorted(pids)) + (
        dmr_params['d_max'],
        dmr_params['n_min'],
        dmr_params['delta_m_min'],
        dmr_params['dmr_test_method'],
        dmr_params['alpha'],
        dmr_params['norm_method']
    ) + tuple(dmr_params['test_kwargs'].items())
    return hash(hash_elements)


def de_results_hash(pids, de_params):
    hash_elements = tuple(sorted(pids)) + tuple([de_params[k] for k in sorted(de_params.keys())])
    return hash(hash_elements)


def expanded_core_sets(venn_set, subgroup_ind):
    """
    Compute the sets that belong to the 'expanded core', which comprises any combination of patients that bridges
    multiple subgroups.
    :param venn_set: As returned by setops.venn_from_arrays
    :param subgroup_ind: Dict, keys are the names of the subgroups, values are boolean arrays with the membership for
    each patient. It's important that the ordering here is the same as that used downstream.
    :return:
    """
    ecs = []
    for k in venn_set:
        this_k = np.array([t for t in k]).astype(bool)
        nmatch = 0
        for grp, grp_idx in subgroup_ind.items():
            if this_k[grp_idx].any():
                nmatch += 1
        if nmatch > 1:
            # add to the expanded core set
            ecs.append(k)
    return ecs


def de_dmr_hash(dat, cols=('gene', 'dmr_cid')):
    """
    Generate a list of hash values from the supplied DataFrame
    :param dat:
    :param cols:
    :return:
    """
    return list(dat.loc[:, cols].itertuples(index=False, name=None))


def annotate_de_dmr_wide_form(data):
    data.insert(0, 'dmr_cluster_id', [t[1] for t in data.index])
    data.insert(0, 'gene', [t[0] for t in data.index])


def expand_dmr_results_table_by_gene_and_cluster(df, drop_geneless=False):
    df.insert(0, 'cluster_id', df.index)

    if not drop_geneless:
        # we need to protect the geneless results by adding dummy values
        g_count = df.genes.apply(len)
        df.loc[g_count == 0, 'genes'] = [(('!', '!'),)] * (g_count == 0).sum()

    # open out the the (gene, relation) lists contained in the `genes` column
    g = pd.DataFrame(
        df.genes.apply(pd.Series).stack().reset_index(level=-1, drop=True),
        columns=['gene_relation']
    )
    g.insert(0, 'relation', g.gene_relation.apply(lambda x: x[1]))
    g.insert(0, 'gene', g.gene_relation.apply(lambda x: x[0]))
    g = g.drop('gene_relation', axis=1).reset_index()

    # pivot table around the relation
    tbl = g.pivot_table(index=['cluster_id', 'gene'], columns='relation', fill_value=0, aggfunc='size').reset_index()

    # extract other attributes from original dataframe
    attrs = df.drop(['cluster_id', 'genes'], axis=1).loc[tbl.cluster_id]

    # combine them
    res = pd.concat((tbl.set_index('cluster_id'), attrs), axis=1).reset_index()

    # if maintaining geneless results, deprotect these now
    if not drop_geneless:
        new_gene_col = res.loc[:, 'gene'].replace('!', np.nan)  # don't use None as the replacement variable!
        res.loc[:, 'gene'] = new_gene_col
        # drop the pivot column, but only if it actually exists
        if '!' in res.columns:
            res.drop('!', axis=1, inplace=True)
        # quick sanity check
        cols = sorted(set(g.relation.unique()).difference(['!']))
        ix = res.gene.isnull()
        if not (res.loc[ix, cols].sum(axis=1) == 0).all():
            logger.error("Some 'geneless' clusters still have gene relations defined.")
            raise ValueError("Some 'geneless' clusters still have gene relations defined.")

    res.index = zip(res.cluster_id, res.gene)
    return res


def patient_and_subgroup_specific_set(pids, subgroups):
    """
    For each PID, get the list of sets that contains that patient and is specific to the right subgroup
    :param pids:
    :param subgroups:
    :return: Dict, one entry per PID. Values are the lists of corresponding sets.
    """
    subgroup_ind = dict([
        (k, pd.Index(pids).isin(v)) for k, v in subgroups.items()
    ])
    candidates = list(setops.binary_combinations(len(pids), include_zero=False))
    # for each PID, what subgroup are they in?
    subgroup_membership = {}
    res = {}
    for i, pid in enumerate(pids):
        res[pid] = []
        for sg, arr in subgroup_ind.items():
            if arr[i]:
                if pid in subgroup_membership:
                    raise ValueError("Patient %s matches two or more subgroups: %s %s." % (
                        str(pid),
                        subgroup_membership[pid],
                        sg
                    ))
                subgroup_membership[pid] = sg
        if pid not in subgroup_membership:
            raise ValueError("Patient %s has no subgroup." % str(pid))

        not_sg_idx = np.where(~subgroup_ind[subgroup_membership[pid]])[0]
        for c in candidates:
            c = np.array([t for t in c])
            if any(c[not_sg_idx] == '1'):
                # Not subgroup specific
                continue
            if c[i] == '1':
                res[pid].append(''.join(c))
    return res


def partial_subgroup_specific(pids, subgroup_ind):
    """
    For each subgroup, get the list of sets that corresponds to >= 2 members of that subgroup and no members of any
    other.
    :param pids: List of PIDs
    :param subgroup_ind: Dict, keyed by subgroup name. Values are np.array of booleans (i.e. boolean indexers). Order
    of indexers must match pids.
    :return: Dict, keyed by subgroup name. Each value is a list of sets.
    """
    ss_part_sets = {}
    candidates = list(setops.binary_combinations_sum_gte(len(pids), 2))
    for grp in subgroup_ind:
        ss_part_sets[grp] = [t for t in candidates if np.array([x for x in t]).astype(int)[~subgroup_ind[grp]].sum() == 0]
    return ss_part_sets


def de_grouped_dispersion(dat, groups, comparisons, min_cpm=1., **de_params):
    cpm = dat.divide(dat.sum(), axis=1) * 1e6
    keep = (cpm > min_cpm).sum(axis=1) > 0
    dat = dat.loc[keep]

    res = differential_expression.run_multiple_de(
        the_data=dat,
        the_groups=groups,
        comparisons=comparisons,
        **de_params
    )

    for the_comparison, this_res in res.items():
        try:
            the_genes = this_res.index
            the_groups = groups[groups.isin(the_comparison)]
            the_cpm = cpm.loc[the_genes, the_groups.index]
            keep = filter.filter_cpm_by_group(
                the_cpm,
                groups,
                min_cpm=min_cpm
            )
            res[the_comparison] = this_res.loc[keep]
        except Exception as exc:
            print repr(exc)

    return res


def de_export_to_ipa(de_wideform, pids):
    cols = []
    for p in pids:
        cols.append("%s logFC" % p)
        cols.append("%s FDR" % p)
    ipa_export = pd.DataFrame(index=de_wideform.index, columns=cols)
    for pid in pids:
        ix = de_wideform.loc[:, pid] == 'Y'
        ipa_export.loc[ix, "%s logFC" % pid] = de_wideform.loc[ix, "%s_logFC" % pid]
        ipa_export.loc[ix, "%s FDR" % pid] = de_wideform.loc[ix, "%s_FDR" % pid]
    return ipa_export


def export_dmr_data_for_ipa(results_tbl):
    # output data for IPA
    all_genes = set()
    for v in results_tbl.values():
        all_genes.update(v.gene.dropna())
    # convert to Ensembl and drop any that cannot be found
    genes_for_ipa = annotation_gene_to_ensembl.gene_to_ens(all_genes).dropna().sort_values()
    # data import is made easier if the columns contain letters as well as numbers
    ipa_export = pd.DataFrame(index=genes_for_ipa.index, columns=["Patient %s" % p for p in results_tbl.keys()])

    for pid, this in results_tbl.items():
        ix = this.gene.isin(genes_for_ipa.index)
        this = this.loc[ix]
        this = this.loc[~this.padj.isnull()]
        this = this.groupby('gene').padj.min()
        ipa_export.loc[this.index, "Patient %s" % pid] = this
    # easier to leave non-significant values as NaN (missing)
    # alternative: set p value to 1.0
    # ipa_export = ipa_export.fillna(1.)

    # replace the index with the Ensembl IDs
    ipa_export = ipa_export.set_index(genes_for_ipa)
    return ipa_export


if __name__ == "__main__":
    """
    The contents of this script (i.e. not including the functions) have been transferred to
    hgic_final.two_strategies_grouped_dispersion.

    Functions are left because other scripts rely on them. Will eventually move these too.
    """
    pass