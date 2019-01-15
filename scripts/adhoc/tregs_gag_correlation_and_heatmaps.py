from scripts.hgic_final import consts, two_strategies_grouped_dispersion as tsgd
from rnaseq import loader
import numpy as np
from scipy import stats
import pandas as pd
from stats import nht
from utils import output, setops
from settings import HGIC_LOCAL_DIR
import os
import pickle
import references
from matplotlib import pyplot as plt
import seaborn as sns


if __name__ == "__main__":
    outdir = output.unique_output_dir()
    pids = consts.PIDS
    DE_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'de')

    eps = .1  # offset for log transform

    target_gene = 'CD274'
    target_ens = references.gene_symbol_to_ensembl(target_gene)

    # load Salmon data

    obj_cc = loader.load_by_patient(pids, source='salmon')
    ix = obj_cc.meta.index.isin(consts.S1_RNASEQ_SAMPLES)
    obj_cc.filter_samples(ix)
    obj_ff = loader.load_by_patient(pids, source='salmon', type='ffpe')

    obj_cc.meta.insert(
        0,
        'pid',
        obj_cc.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>')
    )

    # pull out logged TPM
    log2_tpm_cc = np.log2(obj_cc.data + eps)
    log2_tpm_ff = np.log2(obj_ff.data + eps)

    # load DE results
    the_hash = tsgd.de_results_hash(obj_cc.meta.index.tolist(), consts.DE_PARAMS)
    filename = 'de_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        with open(fn, 'rb') as f:
            de_res_full_s1 = pickle.load(f)
    else:
        raise IOError("No pre-computed results file found: %s" % fn)

    # load genes in GAG-related IPA pathways
    pathways = [
        'Chondroitin Sulfate Biosynthesis',
        'Dermatan Sulfate Biosynthesis',
        'Dermatan Sulfate Biosynthesis (Late Stages)',
        'Chondroitin Sulfate Biosynthesis (Late Stages)',
        # 'Caveolar-mediated Endocytosis Signaling',
    ]

    # pathways by DE
    fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current',
        'core_pipeline',
        'rnaseq',
        's0_individual_patients_direct_comparison',
        'ipa',
        'pathways',
        'full_de_all.xls'
    )
    # result is a dict keyed by PID
    res = pd.read_excel(fn, sheet_name=None, index_col=0)

    # get the list of genes for each pathway
    pathway_genes = {}
    for p in pathways:
        pathway_genes[p] = set()
        for this in res.values():
            if p in this.index:
                pathway_genes[p].update(this.loc[p, 'genes'].split(','))

    # collapse to a single list
    combined_pathway_genes = sorted(setops.reduce_union(*pathway_genes.values()))
    # add target gene for heatmap
    lookup_gene = combined_pathway_genes + [target_gene]
    lookup_ens = references.gene_symbol_to_ensembl(lookup_gene)

    # lookup and plot in patients (GIC)
    log_gic_cc = log2_tpm_cc[consts.S1_RNASEQ_SAMPLES_GIC]
    log_insc_cc = log2_tpm_cc[consts.S1_RNASEQ_SAMPLES_INSC]
    log_dat_cc = pd.concat((log_gic_cc, log_insc_cc), axis=1).loc[lookup_ens]
    log_dat_cc.index = lookup_gene

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.heatmap(
        log_dat_cc,
        cmap='Reds',
        vmin=0.,
        ax=ax,
        cbar_kws={'label': r'$\log_2(\mathrm{TPM} + 0.1)$'}
    )
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    fig.tight_layout()

    fig.savefig(os.path.join(outdir, "log2_tpm_cc_gag_gene_list_heatmap.png"), dpi=200)

    # same but for logFC
    de_logfc = pd.concat(
        [de_res_full_s1[pid].reindex(lookup_ens)['logFC'] for pid in pids],
        axis=1
    )
    de_logfc.columns = pids
    de_logfc.index = lookup_gene

    de_fdr = pd.concat(
        [de_res_full_s1[pid].reindex(lookup_ens)['FDR'] for pid in pids],
        axis=1
    )
    de_fdr.columns = pids
    de_fdr.index = lookup_gene

    yy = de_fdr.astype(float).fillna(1.).values
    xx = de_logfc.astype(float).mask(yy > 0.01)
    mm = np.ma.masked_where(xx.isnull(), xx)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.heatmap(
        xx,
        cmap='RdBu',
        vmin=-6.,
        vmax=6.,
        ax=ax,
        edgecolor='none',
        linewidth=0.,
        cbar_kws={'label': r'$\log_2(\mathrm{FC})$'}
    )
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "log2_fc_gag_gene_list_heatmap.png"), dpi=200)

    # run correlation with TPM values
    x = log_gic_cc.loc[target_ens]
    y = log_gic_cc.loc[lookup_ens]
    y.index = lookup_gene
    y.drop(target_gene, axis=0, inplace=True)

    corr_pdl1_pvals = {}
    corr_pdl1_r = {}

    for gene, vals in y.iterrows():
        corr_pdl1_r[gene], corr_pdl1_pvals[gene] = nht.spearman_exact(x, vals.loc[x.index], nperm=10000)
    corr_pdl1_r = pd.Series(corr_pdl1_r)
    corr_pdl1_pvals = pd.Series(corr_pdl1_pvals)

    # plot the correlation
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sm = plt.cm.ScalarMappable(plt.Normalize(vmin=-1, vmax=1), cmap=plt.cm.get_cmap('RdBu_r'))
    cols = [sm.to_rgba(u) for u in corr_pdl1_r]
    lws = [2. if u < 0.05 else .5 for u in corr_pdl1_pvals]
    ax.bar(range(y.shape[0]), corr_pdl1_r, color=cols)
    ax.bar(range(y.shape[0]), corr_pdl1_r, edgecolor='k', facecolor='none', linewidth=lws)
    ax.set_xticks(range(y.shape[0]))
    ax.set_xticklabels(y.index, rotation=90)
    ax.set_ylim([-.6, .6])
    ax.set_ylabel('Spearman rank correlation coefficient')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "spearman_correlation_with_pdl1.png"), dpi=200)

    # plot scatterplots for the two significant ones
    to_plot = ['CHSY3', 'HS6ST3']
    for g in to_plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(x, y.loc[g, x.index], s=30, edgecolor='k', linewidths=1., facecolor='lightblue')
        ax.set_xlabel(r"$\log_2(\mathrm{%s})$" % target_gene)
        ax.set_ylabel(r"$\log_2(\mathrm{%s})$" % g)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "pdl1_vs_%s_scatter.png" % g), dpi=200)










