import multiprocessing as mp
import os
import pandas as pd
import numpy as np
import pickle
import itertools
from utils import output, setops, excel, log, go_analysis
from methylation import dmr, process, loader as methylation_loader, annotation_gene_to_ensembl
from rnaseq import loader as rnaseq_loader, differential_expression, general, filter
from integrator import rnaseq_methylationarray
from analysis import cross_comparison
from load_data import loader
from plotting import venn
from matplotlib import pyplot as plt, text, patches
import seaborn as sns
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd
from scripts.hgic_final import consts


if __name__ == '__main__':
    """
    Run a cross-comparison on our syngeneic patient lines.
    This allows us to use (non-matching) healthy lines as if they were an external reference.
    """

    # define parameters
    min_cpm = 1.
    pids = consts.PIDS
    de_params = consts.DE_PARAMS
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    # file location
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')
    DE_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'de')

    # boilerplate
    outdir = output.unique_output_dir()
    logger = log.get_console_logger()

    ## DE

    # load data (ours only, no references)
    rnaseq_obj = rnaseq_loader.load_by_patient(pids, include_control=False)
    rnaseq_obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES)

    dat_s1 = rnaseq_obj.data
    meta_s1 = rnaseq_obj.meta

    the_hash = tsgd.de_results_hash(meta_s1.index.tolist(), de_params)
    filename = 'de_results_s1_cross_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 cross-comparison DE results from %s", fn)
        with open(fn, 'rb') as f:
            de_res = pickle.load(f)
    else:
        groups = pd.Series(index=meta_s1.index)
        comparisons = {}
        for pid in pids:
            groups[groups.index.str.contains('GBM') & groups.index.str.contains(pid)] = "GBM%s" % pid
            groups[groups.index.str.contains('NSC') & groups.index.str.contains(pid)] = "iNSC%s" % pid
        for gic_pid, insc_pid in itertools.product(pids, pids):
            comparisons[("GBM%s" % gic_pid, "iNSC%s" % insc_pid)] = "GBM%s - iNSC%s" % (gic_pid, insc_pid)
        de_res = tsgd.de_grouped_dispersion(
            dat_s1,
            groups,
            comparisons,
            min_cpm=min_cpm,
            return_full=True,
            **de_params
        )
        with open(fn, 'wb') as f:
            pickle.dump(de_res, f)

        logger.info("Saved S1 DE results to %s", fn)

    # extract only significant DE genes
    de_res_sign = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res.items()])

    de_count = pd.DataFrame(index=pd.Index(pids, name='GIC'), columns=pd.Index(pids, name='iNSC'))
    for gic_pid, insc_pid in itertools.product(pids, pids):
        de_count.loc[gic_pid, insc_pid] = len(de_res_sign[("GBM%s" % gic_pid, "iNSC%s" % insc_pid)])

    # for each GIC line: get DE genes in syngeneic comparison but NOT in any cross-comparison
    n_syn_only = pd.DataFrame(index=pd.Index(pids, name='GIC'), columns=pd.Index(pids, name='iNSC'))
    syn_only = {}
    for gic_pid, insc_pid in itertools.product(pids, pids):
        # this_syn isn't necessarily syngeneic, but we're acting as if it were here
        this_syn = de_res_sign[("GBM%s" % gic_pid, "iNSC%s" % insc_pid)]
        others = [de_res_sign[("GBM%s" % gic_pid, "iNSC%s" % p)] for p in pd.Index(pids).drop(insc_pid)]
        # 'syn only' index
        this_so_ix = this_syn.index.difference(setops.reduce_union(*[t.index for t in others]))
        syn_only[("GBM%s" % gic_pid, "iNSC%s" % insc_pid)] = this_syn.loc[this_so_ix]
        n_syn_only.loc[gic_pid, insc_pid] = this_so_ix.size

    true_syn_only = dict([
        (p, syn_only[("GBM%s" % p, "iNSC%s" % p)]) for p in pids
    ])

    # export to list
    excel.pandas_to_excel(true_syn_only, os.path.join(outdir, "de_only_in_syngeneic.xlsx"))

    # upset plots x 10, one for each 'row' (fix GIC, change comparator)
    ss = setops.specific_sets(pids)

    for gic_pid in pids:
        for_plot = [de_res_sign[('GBM%s' % gic_pid, 'iNSC%s' % pid)].index for pid in pids]
        set_colours = [
            ('True syngeneic only', {'sets': [ss[gic_pid]], 'colour': '#fc3535'}),
            ("Spurious 'syngeneic only'", {'sets': [ss[pid] for pid in pd.Index(pids).drop(gic_pid)], 'colour': '#3c9100'})

        ]
        plt_dict = venn.upset_set_size_plot(
            for_plot,
            pids,
            set_colours=set_colours,
            n_plot=20,
            include_singletons=True
        )
        plt_dict['figure'].savefig(os.path.join(outdir, "upset_vary_comparator_%s.png" % gic_pid), dpi=200)

    # GO analysis

    # only include the following namespaces
    namespaces = {'BP', 'MF'}
    alpha = 0.05

    bg_genes = setops.reduce_union(*[t.index for t in de_res_sign.values()])
    go_obj = go_analysis.GOAnalysis()
    gc = go_analysis.ens_to_entrez(bg_genes)
    go_obj.set_gene_conversion(gc)
    go_obj.set_bg_genes(bg_genes)

    # for each GIC / comparator, run GO on the 'syngeneic only' genes
    # this results in VERY few terms
    go_syn_only = {}
    for gic_pid, insc_pid in itertools.product(pids, pids):
        k = ("GBM%s" % gic_pid, "iNSC%s" % insc_pid)
        go_syn_only[k] = go_obj.run_one(syn_only[k].index)

    for k, df in go_syn_only.items():
        df.to_csv(os.path.join(outdir, "all_go_terms_%s_%s.csv" % k))

    # filter - but this basically leaves nothing!
    go_syn_only_filtered = {}
    for k, df in go_syn_only.items():
        this_ix = (df.p_bonferroni <= alpha) & (df.NS.isin(namespaces)) & (df.enrichment == 'e')
        df_filt = df.loc[this_ix]
        # include bottom-most nodes only
        ix = []
        for go_id in df_filt.index:
            ix.append(len(go_obj.obo[go_id].get_all_children()) == 0)
        go_syn_only_filtered[k] = df_filt.loc[ix]
