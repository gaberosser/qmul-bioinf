"""
This script is designed to report the (aggregated) beta values for each patient in GBM and NSC (and FFPE?)

The unit of analysis is a gene, and we use concordant DMRs to link probes to one or more genes. Only concordant DMRs
with significant differential methylation between  iNSC and GBM are retained. Of these, we aggregate over all samples
(typically just 2 replicates) and all probes (however many are in the combined DMRs) to get a single beta value for
that gene and cell type.

For simplicity, we'll start the analysis with a previously generated list.
"""

import pandas as pd
from methylation import loader
import numpy as np
import os
import collections
from utils import output, setops
from methylation import dmr, process
from scripts.hgic_de_dmr.two_strategies_grouped_dispersion import dmr_results_hash, paired_dmr, load_methylation


if __name__ == "__main__":
    outdir = output.unique_output_dir("hgic_dmr_concordant_beta")

    # this file is generated by two_strategies_grouped_dispersion
    fn = '~/Dropbox/research/qmul/results/hgic_two_strategies/2018-04-10/s1/full_de_dmr_concordant.xlsx'
    dat = pd.read_excel(fn)
    pids = dat.columns[dat.columns.str.contains(r'^[0-9][0-9][0-9]$')].tolist()

    # also load or compute the DMR results and probe data
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')
    dmr_params = {
        'd_max': 400,
        'n_min': 6,
        'delta_m_min': 1.4,
        'alpha': 0.01,
        'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
        'test_kwargs': {},
    }
    norm_method_s1 = 'swan'

    ############
    # 1: FFPE  #
    ############

    # in this case, we want the median beta value over all probes that are associated with a given gene
    # we'll exclude those associated with gene body only

    ffpe_obj = loader.load_by_patient(pids, type='ffpe', norm_method=norm_method_s1)


    ############
    # 2: DMRs  #
    ############

    # load raw methylation data
    me_obj, anno = load_methylation(pids, norm_method=norm_method_s1)
    me_data = process.beta_from_m(me_obj.data)
    me_meta = me_obj.meta

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    the_hash = dmr_results_hash(me_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        dmr_res_s1 = paired_dmr(me_data, me_meta, anno, pids, dmr_params)
        # Save DMR results to disk
        dmr_res_s1.to_pickle(fn, include_annotation=False)
        print "Saved DMR results to %s" % fn

    median_beta = collections.OrderedDict()

    for pid in pids:
        this_res = {}
        this_res_median = {}
        this_res_cluster_count = collections.Counter()
        this_res_probe_count = collections.Counter()

        the_dmr_res = dmr_res_s1[pid]
        idx = dat[pid] == 'Y'
        the_genes = dat.loc[idx, 'gene'].values
        the_clusters = dat.loc[idx, 'dmr_cluster_id'].values
        for i, c in enumerate(the_clusters):
            g = the_genes[i]
            a = the_dmr_res.clusters[c]

            this_res_cluster_count[g] += 1
            this_res_probe_count[g] += len(a.pids)

            for t in ['GBM', 'DURA']:
                this_res.setdefault(t, {})
                the_dat = me_data.loc[a.pids, me_data.columns.str.contains('%s%s' % (t, pid))]
                if g in this_res[t]:
                    this_res[t][g] = pd.concat((this_res[t][g], the_dat), axis=0)
                else:
                    this_res[t][g] = the_dat

        for t in ['GBM', 'DURA']:
            for g in this_res[t]:
                this_res_median.setdefault(t, {})[g] = np.median(this_res[t][g].values)

        median_beta['GBM%s' % pid] = pd.Series(this_res_median['GBM'])
        median_beta['NSC%s' % pid] = pd.Series(this_res_median['DURA'])

    all_genes = sorted(setops.reduce_union(*[t.index for t in median_beta.values()]))
    all_median = pd.DataFrame(index=all_genes, columns=median_beta.keys())

    i = 2
    for pid in pids:
        this_genes = median_beta['GBM%s' % pid].index
        all_median.loc[this_genes, 'GBM%s' % pid] = median_beta['GBM%s' % pid]
        all_median.loc[this_genes, 'NSC%s' % pid] = median_beta['NSC%s' % pid]
        d = pd.Series(index=all_genes)
        d.loc[this_genes[
            (median_beta['GBM%s' % pid] > median_beta['NSC%s' % pid]).values
        ]] = 'Hyper'
        d.loc[this_genes[
            (median_beta['GBM%s' % pid] < median_beta['NSC%s' % pid]).values
        ]] = 'Hypo'
        all_median.insert(i, '%s_direction' % pid, d)
        i += 3

    all_median.to_excel(os.path.join(outdir, 'median_beta_by_gene.xlsx'))