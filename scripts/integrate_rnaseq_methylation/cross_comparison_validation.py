import multiprocessing as mp
import os
import re
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import references
from rnaseq import filter, differential_expression
from settings import LOCAL_DATA_DIR
from utils import output, setops
from load_data import rnaseq_data


def add_gene_symbols(df):
    """
    Add gene symbols to the DataFrame df which is indexed by Ensembl IDs
    """
    gs = references.ensembl_to_gene_symbol(df.index)
    # resolve any duplicates arbitrarily (these should be rare)
    gs = gs.loc[~gs.index.duplicated()]
    df.insert(0, 'Gene Symbol', gs)


def add_fc_direction(df):
    direction = pd.Series(index=df.index, name='Direction')
    direction.loc[df.logFC < 0] = 'down'
    direction.loc[df.logFC > 0] = 'up'
    df.insert(df.shape[1], 'Direction', direction)


def run_one_de(the_data, the_groups, the_comparison, lfc=1, fdr=0.01, method='QLGLM'):
    if method == 'QLGLM':
        the_contrast = "%s - %s" % (the_comparison[0], the_comparison[1])
        res = differential_expression.edger_glmqlfit(the_data, the_groups, the_contrast, lfc=lfc, fdr=fdr)
    elif method == 'GLM':
        the_contrast = "%s - %s" % (the_comparison[0], the_comparison[1])
        res = differential_expression.edger_glmfit(the_data, the_groups, the_contrast, lfc=lfc, fdr=fdr)
    elif method == 'exact':
        res = differential_expression.edger_exacttest(the_data, the_groups, pair=the_comparison[::-1], lfc=lfc, fdr=fdr)

    add_gene_symbols(res)
    add_fc_direction(res)

    return res



def compute_cross_de(rnaseq_obj, pids, external_references=(('GIBCO', 'NSC'),), lfc=1, fdr=0.01, method='QLGLM',
                     njob=None):
    """
    Compute DE between every patient GBM sample and every _other_ healthy patient sample, in addition to paired DE.
    We can also include one or more external references (e.g. Gibco, the default).
    :param rnaseq_obj:
    :param pids:
    :param external_references:
    :param lfc:
    :param fdr:
    :param method:
    :param njob:
    :return:
    """
    if method not in {'QLGLM', 'GLM', 'exact'}:
        raise NotImplementedError("Unsupported method.")
    # if njob is None:
    #     njob = mp.cpu_count()
    # if njob > 1:
    #     pool = mp.Pool(njob)

    de = {}

    for pid in pids:

        # # usual paired comparison
        # the_idx = rnaseq_obj.meta.index.str.contains(pid)
        # the_data = rnaseq_obj.data.loc[:, the_idx]
        # the_data = filter.filter_by_cpm(the_data, min_n_samples=1)
        #
        # the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
        # the_comparison = ['GBM', 'iNSC']
        # de_matched[pid] = run_one_de(the_data, the_groups, the_comparison, lfc=lfc, fdr=fdr, method=method)

        # cross comparison
        for pid2 in pids:
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, 'type'] == 'GBM')) | \
                      (rnaseq_obj.meta.index.str.contains(pid2) & (rnaseq_obj.meta.loc[:, 'type'] == 'iNSC'))
            the_data = rnaseq_obj.data.loc[:, the_idx]
            the_data = filter.filter_by_cpm(the_data, min_n_samples=1)
            the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
            the_comparison = ['GBM', 'iNSC']

            de[(pid, pid2)] = run_one_de(the_data, the_groups, the_comparison, lfc=lfc, fdr=fdr, method=method)

        # external reference comparison
        for er, er_type in external_references:
            the_idx = (rnaseq_obj.meta.index.str.contains(pid) & (rnaseq_obj.meta.loc[:, 'type'] == 'GBM')) | \
                      (rnaseq_obj.meta.index.str.contains(er) & (rnaseq_obj.meta.loc[:, 'type'] == er_type))
            the_data = rnaseq_obj.data.loc[:, the_idx]
            the_data = filter.filter_by_cpm(the_data, min_n_samples=1)
            the_groups = rnaseq_obj.meta.loc[the_idx, 'type'].values
            the_comparison = ['GBM', er_type]
            de[(pid, er)] = run_one_de(the_data, the_groups, the_comparison, lfc=lfc, fdr=fdr, method=method)


    return de


if __name__ == "__main__":
    # if this is specified, we load the DMR results from a JSON rather than recomputing them to save time
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'integrate_rnaseq_methylation')

    outdir = output.unique_output_dir("cross_validate_de_dmr", reuse_empty=True)
    ref_name = 'GIBCONSC_P4'
    # all n=2 samples and RTK II samples
    pids = ['017', '019', '030', '031', '050', '054']
    cmap = 'RdYlGn_r'

    de_params = {
        'lfc': 1,
        'fdr': 0.01,
        'method': 'GLM'
    }

    # dmr_params = {
    #     'core_min_sample_overlap': 3,  # 3 / 4 samples must match
    #     'd_max': 400,
    #     'n_min': 6,
    #     'delta_m_min': 1.4,
    #     'fdr': 0.01,
    #     'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
    #     'test_kwargs': {},
    #     'n_jobs': 8,
    # }

    # Load RNA-Seq from STAR
    rnaseq_obj = rnaseq_data.load_by_patient(pids, annotate_by='Ensembl Gene ID')
    # discard unmapped, etc
    rnaseq_obj.data = rnaseq_obj.data.loc[rnaseq_obj.data.index.str.contains('ENSG')]
    rnaseq_obj.meta = rnaseq_obj.meta.loc[~rnaseq_obj.meta.index.str.contains('IPSC')]
    rnaseq_obj.data = rnaseq_obj.data.loc[:, rnaseq_obj.meta.index]

    # load RNA-Seq from Salmon (for normalised comparison)
    salmon_dat = rnaseq_data.load_salmon_by_patient_id(pids)
    idx = salmon_dat.index.str.replace(r'.[0-9]+$', '')
    salmon_dat.index = idx
    fn = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.p10.release90',
                      'gene_to_transcript.txt')
    gene_transcript = pd.read_csv(fn, header=0, sep='\t').set_index('Transcript stable ID')

    # aggregate to gene level
    genes = gene_transcript.loc[salmon_dat.index, 'Gene stable ID']
    salmon_dat = salmon_dat.groupby(genes).sum()

    # discard unmapped, etc
    salmon_dat = salmon_dat.loc[:, ~salmon_dat.columns.str.contains('IPSC')]


    de_res = compute_cross_de(rnaseq_obj, pids, **de_params)

    # counts of DE genes
    de_counts = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    for pid in pids:
        for pid2 in pids + ['GIBCO']:
            de_counts.loc[pid, pid2] = de_res[(pid, pid2)].shape[0]

    # now we need to compare the paired results with every other result (Gibco and other iNSC)
    pair_only = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    ref_only = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    pair_and_ref_concordant = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    pair_and_ref_discordant = pd.DataFrame(index=pids, columns=pids + ['GIBCO'])
    # loop over GBM samples
    for pid in pids:
        # syngeneic comparison
        the_pair = de_res[(pid, pid)]

        # loop over (i)NSC samples
        # when this is the same as the syngeneic comparison, there will (obviously) be no 'pair only' or 'ref only'
        # genes!
        for pid2 in pids + ['GIBCO']:
            the_ref = de_res[(pid, pid2)]
            the_sets, _ = setops.venn_from_arrays(the_pair.index, the_ref.index)
            pair_only.loc[pid, pid2] = the_sets['10']
            ref_only.loc[pid, pid2] = the_sets['01']
            # for overlapping genes: separate based on direction (matching or non matching)
            the_conc_idx = (the_pair.loc[the_sets['11']].Direction == the_ref.loc[the_sets['11']].Direction)
            pair_and_ref_concordant.loc[pid, pid2] = the_pair.loc[the_sets['11']].loc[the_conc_idx].index
            pair_and_ref_discordant.loc[pid, pid2] = the_pair.loc[the_sets['11']].loc[~the_conc_idx].index

    # can get counts like this
    po_counts = pair_only.applymap(len)
    ro_counts = ref_only.applymap(len)

    # look at intersection of Gibco and all others for a given GBM
    po_int_gibco = pd.DataFrame(index=pair_only.index, columns=pair_only.columns)
    for pid in pids:
        the_ref = pair_only.loc[pid].iloc[-1]
        po_int_gibco.loc[pid] = pair_only.loc[pid].apply(lambda x: sorted(set(x).intersection(the_ref)))
    po_pct_overlap_with_gibco = po_int_gibco.applymap(len) / pair_only.applymap(len) * 100.

    # now look at it the other way: what is present in X vs Y_i that isn't in X vs any other Y?
    po_diff = pd.DataFrame(index=pair_only.index, columns=pair_only.columns)
    for pid in pids:
        for pid2 in pair_only.columns:
            the_ref = pair_only.loc[pid, pid2]
            all_else = pair_only.loc[pid, pair_only.columns != pid2]
            union_all_else = reduce(set.union, all_else, set())
            po_diff.loc[pid, pid2] = sorted(set(the_ref).difference(union_all_else))

    # get counts like this
    po_diff.applymap(len)

    # this computes the set of genes that ALWAYS appears in the pair_only comparison for each possible ref
    po_each = [
        sorted(
            reduce(lambda x, y: set(x).intersection(y), po_diff.loc[~po_diff.index.str.contains(pid), pid])
        ) for pid in pids + ['GIBCO']
    ]
    po_each = pd.Series(po_each, index=pids + ['GIBCO'])

    # get the genes that consistently differ in the pair comparison only and NOT in Gibco (across all patients)
    # these will have an expression pattern in Gibco similar to GBM, so that they do NOT appear
    # po_gibco_diff = sorted(reduce(lambda x, y: set(y).intersection(x), po_diff.loc[:, 'GIBCO']))
    po_gibco_diff = po_each.loc['GIBCO']
    po_gibco_diff_gs = references.ensembl_to_gene_symbol(po_gibco_diff)
    po_gibco_diff_gs = po_gibco_diff_gs.where(~po_gibco_diff_gs.isnull(), po_gibco_diff)

    po_dat = rnaseq_obj.data.loc[po_gibco_diff]
    po_dat.index = po_gibco_diff_gs
    po_dat = np.log2(po_dat + 1)

    # po_dat = salmon_dat.loc[po_gibco_diff]
    # po_dat.index = po_gibco_diff_gs
    # # dropna() here loses one gene - LINC01090 / ENSG00000231689
    # # all others are present
    # po_dat = np.log2(po_dat.dropna() + 0.01)

    # rearrange columns
    cols = (
        po_dat.columns[po_dat.columns.str.contains('GBM')].tolist() +
        ['GIBCO_NSC_P4'] +
        po_dat.columns[po_dat.columns.str.contains('DURA')].tolist()
    )
    po_dat = po_dat.loc[:, cols]
    # insert spacing columns
    idx = np.where(po_dat.columns.str.contains('GIBCO'))[0][0]
    po_dat.insert(idx, '', np.nan)
    po_dat.insert(idx + 2, ' ', np.nan)

    fig = plt.figure(figsize=(7, 10))
    ax = fig.add_subplot(111)
    ax = sns.heatmap(po_dat, cmap=cmap, ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "consistently_in_pair_only.png"), dpi=200)

    # get the genes that consistently appear in the Gibco reference comparison only and NOT in any other reference
    # these will have a different expression pattern in Gibco to GBM (while the level in iNSC will not differ from GBM)
    ro_diff = pd.DataFrame(index=ref_only.index, columns=ref_only.columns)
    for pid in pids:
        for pid2 in ref_only.columns:
            the_ref = ref_only.loc[pid, pid2]
            all_else = ref_only.loc[pid, ref_only.columns != pid2]
            union_all_else = reduce(set.union, all_else, set())
            ro_diff.loc[pid, pid2] = sorted(set(the_ref).difference(union_all_else))

    # get counts like this
    ro_diff.applymap(len)

    # this computes the set of genes that ALWAYS appears in the ref_only comparison for each possible ref
    ro_each = [
        sorted(
            reduce(lambda x, y: set(x).intersection(y), ro_diff.loc[~ro_diff.index.str.contains(pid), pid])
        ) for pid in pids + ['GIBCO']
    ]
    ro_each = pd.Series(ro_each, index=pids + ['GIBCO'])

    # ro_gibco_diff = sorted(reduce(lambda x, y: set(y).intersection(x), ro_diff.loc[:, 'GIBCO']))
    ro_gibco_diff = ro_each.loc['GIBCO']
    ro_gibco_diff_gs = references.ensembl_to_gene_symbol(ro_gibco_diff)
    # the lincRNA symbols are missing, so keep ENSG for those
    ro_gibco_diff_gs = ro_gibco_diff_gs.where(~ro_gibco_diff_gs.isnull(), other=ro_gibco_diff)

    ro_dat = rnaseq_obj.data.loc[ro_gibco_diff]
    ro_dat.index = ro_gibco_diff_gs
    ro_dat = np.log2(ro_dat + 1)

    # ro_dat = salmon_dat.loc[ro_gibco_diff]
    # ro_dat.index = ro_gibco_diff_gs
    # # dropna() here
    # # all others are present
    # ro_dat = np.log2(ro_dat.dropna() + 0.01)

    # rearrange columns
    cols = (
        ro_dat.columns[ro_dat.columns.str.contains('GBM')].tolist() +
        ['GIBCO_NSC_P4'] +
        ro_dat.columns[ro_dat.columns.str.contains('DURA')].tolist()
    )
    ro_dat = ro_dat.loc[:, cols]

    # insert spacing columns
    idx = np.where(ro_dat.columns.str.contains('GIBCO'))[0][0]
    ro_dat.insert(idx, '', np.nan)
    ro_dat.insert(idx + 2, ' ', np.nan)

    fig = plt.figure(figsize=(7, 10))
    ax = fig.add_subplot(111)
    ax = sns.heatmap(ro_dat, cmap=cmap, ax=ax)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "consistently_in_ref_only.png"), dpi=200)
