import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster import hierarchy
import os
from settings import DATA_DIR
from matplotlib import pyplot as plt, colors
import seaborn as sns
import multiprocessing as mp
from sklearn.decomposition import PCA
from plotting import heatmap, clustering

from utils import setops, output


def one_vs_many_correlation(one, many, method='pearson'):
    if method == 'pearson':
        the_func = stats.pearsonr
    elif method == 'spearman':
        the_func = stats.spearmanr
    else:
        raise NotImplementedError("Unrecognised method %s" % method)

    cor = pd.Series(index=many.index)
    pval = pd.Series(index=many.index)

    for idx, vals in many.iterrows():
        cor[idx], pval[idx] = the_func(one, vals)

    return cor, pval


def pairwise_correlation(many, method='pearson'):
    if method == 'pearson':
        the_func = stats.pearsonr
    elif method == 'spearman':
        the_func = stats.spearmanr
    else:
        raise NotImplementedError("Unrecognised method %s" % method)

    cor = pd.DataFrame(index=many.index, columns=many.index)
    pval = pd.DataFrame(index=many.index, columns=many.index)

    for i, row in enumerate(many.index):
        for j in range(i, len(many.index)):
            col = many.index[j]
            c, p = the_func(many.loc[row], many.loc[col])
            cor.loc[row, col] = c
            cor.loc[col, row] = c
            pval.loc[row, col] = p
            pval.loc[col, row] = p

    return cor, pval


def standardise(data, axis=1):
    """
    Standardise the RMA data in the exon array
    Defined as median centred and divided by the IQR.
    :param data:
    :param axis: The axis to run this - default is rows (probes)
    :return:
    """
    other_axis = 0 if axis == 1 else 1
    return data.subtract(data.median(axis=axis), axis=other_axis).divide(stats.iqr(data, axis=axis), axis=other_axis)


if __name__ == '__main__':
    myc_probe = 3115514
    myc_gene = 'MYC'
    corr_method = 'pearson'
    cross_corr_threshold = 0.5
    within_corr_threshold = 0.75
    alpha = 0.05
    outlier_max_dist = 100. # maximum distance from the centroid before a sample is declared an outlier
    remove_outliers = False

    # these genes were validated, so we'd like to retain them!
    validated_genes = [
        'CXCL1',
        'CXCL2',
        'IL6',
        'IL1RL2',
    ]

    outdir = output.unique_output_dir("cpt_myc")

    fn = os.path.join(DATA_DIR, 'exon_array', 'GSE60982', 'raw', 'GSE60892_HuEx-ALL40-Upload-Transposed.txt.gz')
    rma_data = pd.read_csv(fn, sep='\t', comment='#', header=0, index_col=0)

    ann_fn = os.path.join(DATA_DIR, 'exon_array', 'GSE60982', 'HuEx-1_0-st-v2.na36.hg19.probeset.csv.gz')
    ann = pd.read_csv(ann_fn, sep=',', comment='#', header=0, index_col=0, na_values='---')

    # annotation files
    # latest version from Affymetrix / Thermofisher
    # every probe in this annotation is defined as 'core'
    ann = ann.loc[rma_data.index]
    ann = ann.loc[~ann.gene_assignment.isnull()]
    ann.loc[:, 'gene_assignment'] = ann.loc[:, 'gene_assignment'].str.replace(r' /+ ', '|')
    the_symbols = ann.gene_assignment.str.split('|').apply(lambda x: x[1])

    # get a standardised dataset
    the_data = standardise(rma_data)

    # note to future self: the samples do *not* separate as expected when using the full array representation
    # (see below) - instead we need to run a decomposition method later, once we have selected probes / genes
    # pca = PCA(n_components=10)
    # y = pca.fit_transform(the_data.transpose())

    # get all MYC probes
    myc_probes = the_symbols.index[the_symbols == myc_gene]
    n_myc_probes_orig = len(myc_probes)

    # reduce down to those that correlate with one another
    # form the square pairwise corr matrix, remove rows that don't correlate
    myc_probe_corr, myc_probe_pval = pairwise_correlation(the_data.loc[myc_probes], method=corr_method)
    drop_idx = ((myc_probe_corr < within_corr_threshold) & (myc_probe_pval > alpha)).sum(axis=0) == (n_myc_probes_orig - 1)
    myc_probes = myc_probes[~drop_idx]

    print "Of the original %d MYC probes, we retained %d that correlate well: %s." % (
        n_myc_probes_orig,
        len(myc_probes),
        ', '.join(myc_probes.astype(str))
    )

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.heatmap(myc_probe_corr.astype(float), cmap='RdBu_r', ax=ax)
    ax.set_title('MYC probe %s correlation' % corr_method.capitalize())
    fig.savefig(os.path.join(outdir, 'myc_probe_cross_correlation.png'), dpi=200)

    # find all probes correlated with them (individually, taking the union)
    # this can be slow so we can run in parallel
    jobs = {}
    pool = mp.Pool()
    myc_corr_probes = []

    for p in myc_probes:
        print "Computing correlation with MYC probe %s" % str(p)
        base = the_data.loc[p]
        jobs[p] = pool.apply_async(
            one_vs_many_correlation,
            args=(the_data.loc[p], the_data),
            kwds=dict(method=corr_method)
        )
        # cor, pval = one_vs_many_correlation(the_data.loc[p], the_data, method=corr_method)
        # these_probes = cor.index[(cor.abs() > cross_corr_threshold) & (pval < alpha)]
        # myc_corr_probes.append(these_probes)

    pool.close()
    pool.join()
    for p in myc_probes:
        cor, pval = jobs[p].get(1e4)
        these_probes = cor.index[(cor.abs() > cross_corr_threshold) & (pval < alpha)]
        myc_corr_probes.append(these_probes)

    #  out of interest, what is the overlap between these? (presumably quite high?)
    vs, vc = setops.venn_from_arrays(*myc_corr_probes)

    # union of probes
    keep_probes = setops.reduce_union(*myc_corr_probes)

    print "After comparing all data against each MYC probe, we are left with %d correlated probes" % len(keep_probes)

    genes_corr_with_myc = the_symbols.loc[keep_probes].dropna()
    print "These correspond to %d unique genes." % len(genes_corr_with_myc.unique())

    # check the overlap with validated genes
    overlap = pd.Index(validated_genes).intersection(genes_corr_with_myc.unique())
    if len(overlap) == len(validated_genes):
        print "Good news: all %d validated genes are in the genes shortlist." % len(validated_genes)
    else:
        print "Oh dear, %d of the %d validated genes are NOT in the shortlist (%s)" % (
            len(validated_genes) - len(overlap),
            len(validated_genes),
            ', '.join(pd.Index(validated_genes).difference(genes_corr_with_myc.unique()))
        )

    # PCA decomposition and centroid distance calculcation to spot outliers
    pca = PCA(n_components=2)
    y = pca.fit_transform(the_data.loc[keep_probes].transpose())

    pca_dist_from_centroid = (y ** 2).sum(axis=1) ** .5
    x = the_data.loc[keep_probes]
    dist_from_centroid = (x.subtract(x.mean(axis=1), axis=0) ** 2).sum(axis=0) ** .5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(y[:, 0], y[:, 1], c=dist_from_centroid, cmap='RdYlGn_r')

    # outliers = np.where(dist_from_centroid > 100)[0]
    outliers = np.where(dist_from_centroid > outlier_max_dist)[0]

    for i in outliers:
        ax.text(y[i, 0], y[i, 1], the_data.columns[i], horizontalalignment='right')
    ax.set_xlabel("Principal component 1 (%.2f%%)" % (pca.explained_variance_ratio_[0] * 100.))
    ax.set_ylabel("Principal component 2 (%.2f%%)" % (pca.explained_variance_ratio_[1] * 100.))
    fig.tight_layout()

    fig.savefig(os.path.join(outdir, 'myc_genes_pca_for_outliers.png'), dpi=200)


    if remove_outliers and len(outliers):
        print "Removing %d sample(s) that ha(s/ve) been branded outlier(s): %s" % (
            len(outliers),
            ', '.join(the_data.columns[outliers])
        )
        # make a backup in case we want to go back
        the_data_before_filt = the_data.copy()
        the_data = standardise(rma_data.loc[:, dist_from_centroid <= outlier_max_dist])
    elif len(outliers):
        print "%d sample(s) is/are outlier(s), but we will retain them: %s" % (
            len(outliers),
            ', '.join(the_data.columns[outliers])
        )

    # aggregate by gene (only within the pre-selected probes)
    dat_corr_with_myc_aggr = the_data.loc[keep_probes].groupby(genes_corr_with_myc, axis=0).mean()

    # re-run correlation analysis at the gene level
    base = dat_corr_with_myc_aggr.loc[myc_gene]
    cor_gene, pval_gene = one_vs_many_correlation(base, dat_corr_with_myc_aggr, method=corr_method)

    # reduce to significant and relevant
    keep_genes = cor_gene.index[(cor_gene.abs() > cross_corr_threshold) & (pval_gene < alpha)]

    # remove MYC itself when reporting
    print "Having aggregated these by gene, %d are correlated with %s" % (len(keep_genes.drop(myc_gene)), myc_gene)

    # cluster using this representation of the data
    # force the order of the columns to match the correlation with MYC

    keep_genes_sorted = cor_gene.loc[keep_genes].sort_values(ascending=False).index

    # define a column colour band
    # base_n = 2 * (base - base.min()) / (base.max() - base.min()) - 1.
    # cmap = plt.cm.RdBu_r
    # norm = colors.Normalize(vmin=-1, vmax=0.)
    # sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    # vals = [colors.rgb2hex(sm.to_rgba(t)) for t in base_n]

    norm = colors.Normalize(vmin=base.min(), vmax=base.max())
    sm = plt.cm.ScalarMappable(cmap=plt.cm.gray_r, norm=norm)

    vals = [colors.rgb2hex(sm.to_rgba(t)) for t in base]

    col_colours = pd.DataFrame(vals, index=the_data.columns, columns=['MYC'])

    cg = clustering.plot_clustermap(
        dat_corr_with_myc_aggr.loc[keep_genes_sorted],
        cmap='RdBu_r',
        metric='euclidean',
        method='ward',
        row_cluster=False,
        col_colors=col_colours,
        vmin=-8, vmax=8
    )
    cg.fig.set_size_inches((7., 7.))
    cg.cax.set_ylabel("Gene expression")
    cg.cax.yaxis.set_label_coords(-.7, 0.5)
    cg.gs.update(bottom=0.12, left=0.04, top=0.97, right=0.93)
    cg.savefig(os.path.join(outdir, 'myc_genes_clustermap.png'), dpi=200)
    cg.savefig(os.path.join(outdir, 'myc_genes_clustermap.tiff'), dpi=200)

    # does the clustering partition by MYC expression level?
    dend = cg.dendrogram_col.calculate_dendrogram()
    lkg = cg.dendrogram_col.linkage

    n_clust = 4
    # get cluster ID
    cluster_id = hierarchy.fcluster(lkg, n_clust, criterion='maxclust')
    myc_group_mean = base.groupby(cluster_id).mean()

    grouped_means = pd.Series(index=base.index)
    for i, m in myc_group_mean.iteritems():
        grouped_means.loc[cluster_id == i] = m

    # add this as a second column colour bar
    vals_2 = [colors.rgb2hex(sm.to_rgba(t)) for t in grouped_means]
    col_colours.insert(1, 'Mean MYC', vals_2)

    cg = clustering.plot_clustermap(
        dat_corr_with_myc_aggr.loc[keep_genes_sorted],
        cmap='RdBu_r',
        metric='euclidean',
        method='ward',
        row_cluster=False,
        col_colors=col_colours,
        vmin=-8, vmax=8
    )
    cg.fig.set_size_inches((7., 7.))
    cg.cax.set_ylabel("Gene expression")
    cg.cax.yaxis.set_label_coords(-.7, 0.5)
    cg.gs.update(bottom=0.12, left=0.04, top=0.97, right=0.88)
    cg.savefig(os.path.join(outdir, 'myc_genes_clustermap_with_group_means.png'), dpi=200)
    cg.savefig(os.path.join(outdir, 'myc_genes_clustermap_with_group_means.tiff'), dpi=200)

    # export to excel for downstream analysis
    to_export = pd.DataFrame.from_dict(
        {'pearson_correlation': cor_gene.sort_values(ascending=False)}
    )
    # add probes
    tmp = pd.Series(genes_corr_with_myc.index, index=genes_corr_with_myc.values)
    probe_col = tmp.groupby(tmp.index).apply(lambda x: ','.join([str(t) for t in x])).loc[to_export.index]
    to_export.insert(1, 'probes', probe_col)

    to_export.to_excel(os.path.join(outdir, 'GSE60982_genes_correlated_with_myc.xlsx'))

    jennie_list = [
        'MYC',
        'CXCL2',
        'CXCL1',
        'TNFAIP3',
        'IL8',
        'C17orf47',
        'RAET1L',
        'TEX14',
        'SERPINE1',
        'CCL2',
        'LIF',
        'TMEM49',
        'IL1RL2',
        'ASTL',
        'MSC',
        'CXCL3',
        'NFKBIZ',
        'FOSB',
        'PLK3',
        'G0S2',
        'CHST7',
        'GADD45B',
        'MAP2K3',
        'IL28A',
        'CXCL6',
        'ZNF165',
        'C10orf55',
        'PCK1',
        'VCAM1',
        'PLAUR',
        'B3GNT5',
        'CCNL1',
        'ATP13A3',
        'TIFA',
        'ATF3',
        'SEC14L4',
        'EIF2AK3',
        'IL6',
        'PPP1R15A',
        'WTAP',
        'ADM',
        'ZC3H12A',
        'JUN',
        'ALOXE3',
        'BAG3', # 'NBAG3',
        'FOSL1',
        'BIRC3',
        'CCL20',
        'PLAU',
        'RRAD',
        'KLF4',
        'NFKBID',
        'RHOB',
        'SFRP4',
        'CTNNA3',
        'DUSP2',
        'MAFF',
        'HSPA5',
        'SLC4A7',
        'SPRY2',
        'STX5',
        'CREB5',
        'SLC27A6',
        'WDR74',
        'ANKRD1',
        'NAMPT',
        'IER5',
        'ICAM1',
        'CLCF1', # 'CCLF1',
        'TFPI2', # 'TFP12',
        'SOD2',
        'REL',
        'SLFN5',
        'C14orf43',
        'WEE1',
        'CH25H',
        'NFIL3',
        'NGEF',
        'EIF2A',
        'IRF1',
        'EDN1',
        'CD44',
        'CD69',
        'HSPA6',
        'IL13RA2',
        'EMP1',
        'TGIF1',
        'ELF3',
        'C15orf48',
        'CSRNP1',
        'NFKBIA',
        'DMAP1',
        'NCOA7',
        'IL28A',
        'IL1B',
        'CCL4',
        'SLC2A3',
        'TACSTD2',
        'MT1M',
        'DUSP4',
        'DNAJB1',
        'C1orf94',
        'BACH2',
        'PTMA',
        'NLRP3',
        'CCL4',
        'KLF10',
        'HSPA7',
        'OFD1',
        'SLC2A14',
        'PRMT3',
        'ILDR1',
        'LRRC2',
        'ME1',
        'PTGS2',
        'DDX26B',
        'MT2A',
        'KBTBD11',
        'LDHC',
        'UBC',
        'BHLHE41',
        'ULBP2',
        'CCL3',
        'PIGA', # 'P1GA',
        'FHL1',
        'TSPYL2',
        'NR4A1',
        'JUNB',
        'FAM46B',
        'HES6',
        'DMGDH',
        'ZFP36',
        'EIF1AX',
        'PMP22',
        'MT1X',
        'DUSP1',
        'MT1B',
        'GPR183',
        'HBEGF',
        'FLJ43826',
        'KCNJ5',
        'ENDOU',
        'STAT5A',
        'IER2',
        'VPS54',
        'MT1L',
        'ARRDC3',
        'TMEM130', # 'TMEM',
        'KDM6A',
        'THOC6',
        'NCEH1',
        'TRIB1',
        'SYCP2',
        'KLF6',
        'AOX1',
        'PAMR1',
        'FERMT3',
        'DGKZ',
        'ARL14',
        'C6orf141',
        'RND1',
        'TXK',
        'EHD1',
        'CAST',
        'LAMC2',
        'GK',
        'TAGAP',
        'NAV2',
        'SKIL',
        'MAGEA3',
        'CCDC132',
        'CCDC106',
        'LRP5L',
        'NFIC',
        'GOLT1A',
        'RCVRN',
        'PEX19',
        'LRRC20',
        'ABHD1',
        'RWDD1',
        'LRRN2',
        'PPFIA3',
        'CAPRIN2',
        'HPSE',
        'HDGFRP2',
        'BTN2A1',
        'SFTPB',
        'ADD1',
        'HIC1',
        'NRG2',
        'SDCCAG8',
        'NOVA2',
        'SCAI',
        'CLPP',
        'C9orf86',
        'PPP6R1',
        'ICK',
        'ZNF442',
        'HCG8',
        'EPB41L4B',
        'SIGLEC1',
        'CNNM2',
        'APLN',
        'SERF2',
        'RAB6B',
        'FLJ42875',
        'MLL2',
        'DSC3',
        'SMUG1',
        'ZBTB40',
        'SCAP',
        'DNAJC16',
        'RAB33A',
    ]

    # bizarrely, not all of these genes have a match - in fact, some don't even exist online?!
    # anyway, do the best I can
    j_not_ours = pd.Index(jennie_list).difference(keep_genes)
    ours_not_j = keep_genes.difference(jennie_list)