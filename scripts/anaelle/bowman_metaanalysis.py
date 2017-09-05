import pandas as pd
from settings import DATA_DIR_NON_GIT, DATA_DIR
import os
import references
import datetime
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np


def get_de_tissue_tumour():

    # load tables S1A, S1B: lists of DE genes in healthy vs TAM
    fn_s1a = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S1A.csv')
    fn_s1b = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S1B.csv')

    s1a = pd.read_csv(fn_s1a, header=0, index_col=None)
    s1b = pd.read_csv(fn_s1b, header=0, index_col=None)

    # manual corrections
    s1a.replace('AI414108', 'Igsf9b', inplace=True)
    s1a.replace('Fam101a', 'Rflna', inplace=True)
    s1b.replace('Fam176b', 'Eva1b', inplace=True)
    s1b.replace('Gm14047', 'Il1bos', inplace=True)
    s1b.replace('Gpr114', 'Adgrg5', inplace=True)


    # get DE genes in (MG vs TAM-MG) and (monocytes vs TAM-BMDM)
    mg = s1a.MG.dropna()
    bmdm = s1b.BMDM.dropna()

    # convert gene symbols to ENS ID
    fn = os.path.join(DATA_DIR, 'ensembl', 'mouse', 'mart_export.txt.gz')
    ref = pd.read_csv(fn, sep='\t', index_col=None, header=0).set_index('Gene name')

    mg_ens = ref.loc[mg.values, 'Gene stable ID'].dropna().unique()
    bmdm_ens = ref.loc[bmdm.values, 'Gene stable ID'].dropna().unique()

    # now can paste these into DAVID / gProfile
    print '\n'
    print "DE in MG / TAM-MG"
    print '\n'.join(mg_ens)
    print '\n'
    print "DE in BMDM / TAM-BMDM"
    print '\n'.join(bmdm_ens)
    print '\n'


def ssgsea(sample_data, gene_set, alpha=0.25, norm_by_gene_count=True, return_ecdf=False):
    """
    Run single sample gene set enrichment analysis (ssGSEA) on the supplied data, following the details given in:

    Barbie, D.A., Tamayo, P., Boehm, J.S., Kim, S.Y., Moody, S.E., Dunn, I.F., Schinzel, A.C., Sandy, P., Meylan, E.,
    Scholl, C., et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1.
    Nature 462, 108-112.

    See R/stats/ssgsea.R for a further example (found online)

    :param sample_data: Pandas dataframe, rows are genes and columns are samples
    :param gene_set: Dictionary. Each entry has a key giving the name of the gene set and value giving a list of genes.
    :param alpha: The weighting used to compute the ECDF. Default is 0.25 following Barbie et al. (2009)
    :param return_ecdf: If True, also return the two ECDFs being considered. Useful for plotting?
    :return:
    """
    s = sample_data

    # rank the sample in descending order (1 corresponds to highest expression)
    rs = s.rank(method='average', ascending=False)  # ties resolved by averaging

    # sort in decreasing order
    rs = rs.sort_values(ascending=False)

    # boolean vector for inclusion in gene set
    in_gene_set = rs.index.isin(gene_set).astype(float)
    out_gene_set = (~rs.index.isin(gene_set)).astype(float)

    # ECDF
    x_in = (rs * in_gene_set) ** alpha
    ecdf_in = (x_in.cumsum() / x_in.sum()).values

    # the ECDF for samples out is NOT weighted, which is strange
    ecdf_out = out_gene_set.cumsum() / out_gene_set.sum()

    # if we were to weight it, it would look like:
    # x_out = (rs * out_gene_set) ** alpha
    # ecdf_out = x_out.cumsum() / x_out.sum()

    # enrichment score is the difference in the integrals
    es = (ecdf_in - ecdf_out).sum()

    if norm_by_gene_count:
        es /= float(rs.shape[0])

    if return_ecdf:
        return (es, ecdf_in, ecdf_out)
    else:
        return es



if __name__ == "__main__":

    rnaseq_type = 'counts'
    remove_idh1 = False

    # cutoff for discarding genes
    fpkm_cutoff = 1.
    fpkm_min_samples = 10

    ### load RNA-Seq data annotated by Brennan

    if rnaseq_type == 'counts':
        rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count')
        rnaseq_dat_fn = os.path.join(rnaseq_dir, 'counts.csv')
    elif rnaseq_type == 'fpkm':
        rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count_fpkm')
        rnaseq_dat_fn = os.path.join(rnaseq_dir, 'fpkm.csv')
    else:
        raise NotImplementedError("Unrecognised rnaseq data type")
    rnaseq_meta_fn = os.path.join(rnaseq_dir, 'sources.csv')

    rnaseq_dat_raw = pd.read_csv(rnaseq_dat_fn, header=0, index_col=0)
    rnaseq_meta = pd.read_csv(rnaseq_meta_fn, header=0, index_col=0)

    if remove_idh1:
        # filter IDH1 mutants
        idh1_wt = (~rnaseq_meta.idh1_status.isnull()) & (rnaseq_meta.idh1_status == 'WT')

        rnaseq_meta = rnaseq_meta.loc[idh1_wt]
        rnaseq_dat = rnaseq_dat_raw.loc[:, idh1_wt.values]
    else:
        rnaseq_dat = rnaseq_dat_raw.loc[:, rnaseq_dat_raw.columns.str.contains('TCGA')]

    # add gene symbols for gene signature scoring?
    gs = references.ensembl_to_gene_symbol(rnaseq_dat.index).dropna()
    rnaseq_dat = rnaseq_dat.loc[gs.index]
    rnaseq_dat.index = gs.values

    if rnaseq_type == 'counts':
        # convert to CPM
        rnaseq_dat = rnaseq_dat.divide(rnaseq_dat.sum(axis=0), axis=1) * 1e6

    # plot a histogram showing distribution of expression values
    xx = rnaseq_dat.values.flatten()
    xx = xx[xx != 0.]  # required to avoid MASSIVE spike at precisely 0
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.distplot(np.log10(xx + 1.), bins=200, kde=False, ax=ax)
    if rnaseq_type == 'fpkm':
        ax.set_xlabel('log10(FPKM + 1)')
    elif rnaseq_type == 'counts':
        ax.set_xlabel('log10(CPM + 1)')
    ax.set_ylabel('Density')
    ax.axvline(np.log10(fpkm_cutoff + 1.), ls='--', c='k')

    if rnaseq_dat.index.duplicated().any():
        print "Warning: some gene symbols are duplicated."
        print ', '.join(rnaseq_dat.index[rnaseq_dat.index.duplicated()].tolist())

    ### load microarray data annotated by Verhaak

    marr_dir = os.path.join(DATA_DIR_NON_GIT, 'microarray', 'tcga_gbm', 'verhaak')
    marr_dat_fn = os.path.join(marr_dir, 'scaled_expression.txt')
    marr_meta_fn = os.path.join(marr_dir, 'idh1_status.txt')

    marr_dat = pd.read_csv(marr_dat_fn, sep='\t', header=0, index_col=0)
    marr_meta = pd.read_csv(marr_meta_fn, sep='\t', header=0, index_col=0)

    # rename data columns (manual check for duplicates passed)
    cols = marr_dat.columns.str.replace(r'(TCGA-[0-9]{2}-[0-9]{4})-.*', r'\1')
    marr_dat.columns = cols

    # collapse to intersecting set (only keep samples for which we know the IDH1 status)
    # TODO: we can actually gain 70 more samples if we skip this step and use a different annotation file
    marr_dat = marr_dat.loc[:, marr_meta.index]

    if remove_idh1:
        # filter IDH1 mutants
        marr_dat = marr_dat.loc[:, (marr_meta.IDH1mut == 'wt').values]
        marr_meta = marr_meta.loc[marr_meta.IDH1mut == 'wt']

    # load signatures (Bowman et al.)
    fn_s1a = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S1A.csv')
    fn_s1b = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S1B.csv')
    fn_s4 = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'GSE86573', 'table_S4.csv')

    s1a = pd.read_csv(fn_s1a, header=0, index_col=None)
    s1b = pd.read_csv(fn_s1b, header=0, index_col=None)
    s4 = pd.read_csv(fn_s4, header=0, index_col=None)

    # generate list of orthologs of the relevant gene signatures
    from scripts.agdex_mouse_human_mb_microarray import generate_ortholog_table as got
    orth = got.homologs(got.mouse_tid, got.human_tid)
    orth.set_index('gene_symbol_%d' % got.mouse_tid, inplace=True)

    # use this to generate human gene lists
    all_genes_in_set = set()
    s4_hu = {}
    rna_s4_hu = {}
    marr_s4_hu = {}
    for c in s4.columns:
        s4_hu[c] = orth.loc[s4.loc[:, c].dropna()].dropna().iloc[:, 0].values

    for c in s4_hu:
        this_geneset = set(s4_hu[c].tolist()).intersection(rnaseq_dat.index)
        removed = set(s4_hu[c].tolist()).difference(rnaseq_dat.index)
        if len(removed):
            print "%d genes were removed from RNA-Seq geneset %s as they are not found in the expression data. " \
                  "%d remaining of original %d.\n" % (
                len(removed),
                c,
                len(this_geneset),
                s4[c].dropna().size
            )
            # print  ', '.join(list(removed))
        rna_s4_hu[c] = list(this_geneset)
        all_genes_in_set.update(this_geneset)

        this_geneset = set(s4_hu[c].tolist()).intersection(marr_dat.index)
        removed = set(s4_hu[c].tolist()).difference(marr_dat.index)
        if len(removed):
            print "%d genes were removed from microarray geneset %s as they are not found in the expression data. " \
                  "%d remaining of original %d.\n" % (
                len(removed),
                c,
                len(this_geneset),
                s4[c].dropna().size
            )
            # print ', '.join(list(removed))
        marr_s4_hu[c] = list(this_geneset)

    # remove genes that have no appreciable expression level
    # >=10 samples must have FPKM >= 1
    to_keep = ((rnaseq_dat > fpkm_cutoff).sum(axis=1) > fpkm_min_samples) | (rnaseq_dat.index.isin(all_genes_in_set))
    print "Keeping %d / %d genes that are sufficiently abundant" % (to_keep.sum(), to_keep.size)
    rnaseq_dat = rnaseq_dat.loc[to_keep]

    # run ssGSEA
    # we'll also store the ecdfs
    rna_es = pd.DataFrame(index=s4_hu.keys(), columns=rnaseq_dat.columns)
    rna_ecdf_in = dict()
    rna_ecdf_out = dict()
    for s_name in rnaseq_dat.columns:
        rna_ecdf_in[s_name] = dict()
        rna_ecdf_out[s_name] = dict()
        for g_name in rna_s4_hu:
            rna_es.loc[g_name, s_name], rna_ecdf_in[s_name][g_name], rna_ecdf_out[s_name][g_name] = ssgsea(
                rnaseq_dat.loc[:, s_name],
                rna_s4_hu[g_name],
                alpha=0.25,
                return_ecdf=True
            )

    # scale using the Z transform
    # z = es.subtract(es.mean(axis=1), axis=0).divide(es.std(axis=1), axis=0)
    rna_z = (rna_es - rna_es.values.flatten().mean()) / rna_es.values.flatten().std()

    fig = plt.figure(num="TCGA RNA-Seq")
    ax = fig.add_subplot(111)
    for g_name in s4_hu:
        sns.kdeplot(rna_z.loc[g_name], ax=ax)
    ax.set_xlabel("Normalised ssGSEA score")
    ax.set_ylabel("Density")

    # repeat with the microarray data
    marr_es = pd.DataFrame(index=s4_hu.keys(), columns=marr_dat.columns)
    marr_ecdf_in = dict()
    marr_ecdf_out = dict()
    for s_name in marr_dat.columns:
        marr_ecdf_in[s_name] = dict()
        marr_ecdf_out[s_name] = dict()
        for g_name in s4_hu:
            marr_es.loc[g_name, s_name], marr_ecdf_in[s_name][g_name], marr_ecdf_out[s_name][g_name] = ssgsea(
                marr_dat.loc[:, s_name],
                marr_s4_hu[g_name],
                alpha=0.25,
                return_ecdf=True
            )

    # scale using the Z transform
    # z = es.subtract(es.mean(axis=1), axis=0).divide(es.std(axis=1), axis=0)
    marr_z = (marr_es - marr_es.values.flatten().mean()) / marr_es.values.flatten().std()

    fig = plt.figure(num="TCGA microarray")
    ax = fig.add_subplot(111)
    for g_name in s4_hu:
        sns.kdeplot(marr_z.loc[g_name], ax=ax)
    ax.set_xlabel("Normalised ssGSEA score")
    ax.set_ylabel("Density")

    # now split by subgroup
    subgroups = rnaseq_meta.groupby('expression_subclass').groups
    if remove_idh1:
        try:
            subgroups.pop('G-CIMP')
        except Exception:
            pass

    # boxplot
    bplot = {}
    for g_name in rna_s4_hu:
        bplot[g_name] = {}
        for sg in subgroups:
            bplot[g_name][sg] = rna_z.loc[g_name, subgroups[sg]].values


    tmp = [list(t) for t in bplot['TAM BMDM'].values()]
    plt.figure(num='TAM BMDM')
    plt.boxplot(tmp)
    ax = plt.gca()
    ax.set_xticklabels(subgroups.keys(), rotation=45)

    tmp = [list(t) for t in bplot['TAM MG'].values()]
    plt.figure(num='TAM MG')
    plt.boxplot(tmp)
    ax = plt.gca()
    ax.set_xticklabels(subgroups.keys(), rotation=45)

    tmp = [list(t) for t in bplot['Core BMDM'].values()]
    plt.figure(num='Core BMDM')
    plt.boxplot(tmp)
    ax = plt.gca()
    ax.set_xticklabels(subgroups.keys(), rotation=45)

    tmp = [list(t) for t in bplot['Core MG'].values()]
    plt.figure(num='Core MG')
    plt.boxplot(tmp)
    ax = plt.gca()
    ax.set_xticklabels(subgroups.keys(), rotation=45)

    plt.show()