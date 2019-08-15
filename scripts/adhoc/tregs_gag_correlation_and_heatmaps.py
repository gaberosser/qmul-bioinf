import os
import pickle

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from hgic_consts import NH_ID_TO_PATIENT_ID_MAP
from plotting import common
from rnaseq import loader
from scripts.hgic_final import consts, two_strategies_grouped_dispersion as tsgd
from settings import HGIC_LOCAL_DIR, INTERMEDIATE_DIR
from stats import nht
from utils import output, setops, reference_genomes

if __name__ == "__main__":
    outdir = output.unique_output_dir()
    pids = consts.PIDS
    DE_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'de')

    eps = .1  # offset for log transform

    target_gene = 'CD274'
    target_ens = reference_genomes.gene_symbol_to_ensembl(target_gene)

    # load Salmon data

    obj_cc = loader.load_by_patient(pids, source='salmon')
    ix = obj_cc.meta.index.isin(consts.S1_RNASEQ_SAMPLES)
    obj_cc.filter_samples(ix)

    # add PID to cell culture metadata
    obj_cc.meta.insert(
        0,
        'pid',
        obj_cc.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>')
    )

    obj_ff = loader.load_by_patient(pids, source='salmon', type='ffpe')
    obj_ff.filter_by_sample_name(consts.FFPE_RNASEQ_SAMPLES)

    # add PID to FFPE metadata
    nh_id = obj_ff.meta.index.str.replace(r'(_?)(DEF|SP).*', '')
    p_id = [NH_ID_TO_PATIENT_ID_MAP[t.replace('_', '-')] for t in nh_id]
    obj_ff.meta.insert(0, 'nh_id', nh_id)
    obj_ff.meta.insert(0, 'patient_id', p_id)

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
    lookup_ens = reference_genomes.gene_symbol_to_ensembl(lookup_gene)

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
        cmap='RdBu_r',
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
    y = y[x.index]

    corr_pdl1_pvals = {}
    corr_pdl1_r = {}

    for gene, vals in y.iterrows():
        corr_pdl1_r[gene], corr_pdl1_pvals[gene] = nht.spearman_exact(x, vals.loc[x.index], nperm=10000)
        # corr_pdl1_r[gene], corr_pdl1_pvals[gene] = stats.pearsonr(x, vals.loc[x.index])
        # corr_pdl1_r[gene], corr_pdl1_pvals[gene] = stats.spearmanr(x, vals.loc[x.index])

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

    # plot scatterplots for the more correlated ones
    cmap = common.get_best_cmap(len(consts.PIDS))
    to_plot = corr_pdl1_r.index[corr_pdl1_r.abs() > 0.25]

    # to_plot = ['CHSY3', 'HS6ST3']
    for g in to_plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i, pid in enumerate(pids):
            ix = obj_cc.meta.pid[x.index] == pid
            ax.scatter(
                x.loc[ix],
                y.loc[g, ix],
                facecolors=cmap[i],
                s=30,
                edgecolor='k',
                linewidths=1.,
                label=pid
            )
        ax.set_xlabel(r"$\log_2(\mathrm{%s})$" % target_gene)
        ax.set_ylabel(r"$\log_2(\mathrm{%s})$" % g)
        ax.legend(frameon=True, facecolor='w', framealpha=0.6)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "pdl1_vs_%s_scatter.png" % g), dpi=200)

    plt.close('all')

    # expression of Tregs signature genes across our FFPE samples

    # load the signatures
    xcell_sign_fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current/input_data/xcell',
        'esm3_cell_type_signatures.xlsx'
    )

    xcell_s = pd.read_excel(xcell_sign_fn, header=0, index_row=0)
    xcell_signatures = {}
    for i, row in xcell_s.iterrows():
        xcell_signatures[row.Celltype_Source_ID] = set(row.iloc[2:].dropna().values)

    tregs_signatures = dict([(k, v) for k, v in xcell_signatures.items() if 'tregs' in k.lower()])
    # since there's a lot of complementarity, reduce this down to a single list
    tregs_signature_combined = sorted(setops.reduce_union(*tregs_signatures.values()))
    tregs_signature_combined_ens = reference_genomes.gene_symbol_to_ensembl(tregs_signature_combined)

    # heatmap
    dat = log2_tpm_ff.loc[tregs_signature_combined_ens]
    dat.columns = obj_ff.meta.loc[dat.columns, 'patient_id']
    dat.index = tregs_signature_combined_ens.index
    fig = plt.figure(figsize=(5.5, 8))
    ax = fig.add_subplot(111)
    sns.heatmap(
        dat,
        cmap='Reds',
        vmin=0,
        vmax=8,
        cbar_kws={'label': r'$\log_2(\mathrm{TPM}+0.01)$'}
    )
    ax.xaxis.label.set_visible(False)
    ax.set_ylabel('Treg signature gene')
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "ffpe_xcell_tregs_signature_genes_logtpm.png"), dpi=200)



