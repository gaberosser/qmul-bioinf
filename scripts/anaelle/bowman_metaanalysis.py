import pandas as pd
from settings import DATA_DIR_NON_GIT, DATA_DIR
import os
import references
import datetime


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


if __name__ == "__main__":

    ### load RNA-Seq data annotated by Brennan

    # rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count')
    rnaseq_dir = os.path.join(DATA_DIR_NON_GIT, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'htseq-count_fpkm')
    # rnaseq_dat_fn = os.path.join(rnaseq_dir, 'counts.csv')
    rnaseq_dat_fn = os.path.join(rnaseq_dir, 'fpkm.csv')
    rnaseq_meta_fn = os.path.join(rnaseq_dir, 'sources.csv')

    rnaseq_dat = pd.read_csv(rnaseq_dat_fn, header=0, index_col=0)
    rnaseq_meta = pd.read_csv(rnaseq_meta_fn, header=0, index_col=0)

    # filter IDH1 mutants
    idh1_wt = (~rnaseq_meta.idh1_status.isnull()) & (rnaseq_meta.idh1_status == 'WT')

    rnaseq_meta = rnaseq_meta.loc[idh1_wt]
    rnaseq_dat = rnaseq_dat.loc[:, idh1_wt.values]

    # add gene symbols for gene signature scoring?
    gs = references.ensembl_to_gene_symbol(rnaseq_dat.index).dropna()
    rnaseq_dat = rnaseq_dat.loc[gs.index]
    rnaseq_dat.index = gs.values

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
    s4_hu = {}
    for c in s4.columns:
        s4_hu[c] = orth.loc[s4.loc[:, c].dropna()].dropna().iloc[:, 0].values

    # mg_mo = s1a.MG.dropna()
    # bmdm_mo = s1b.BMDM.dropna()

    # mg_hu = orth.loc[mg_mo].dropna().iloc[:, 0].values
    # bmdm_hu = orth.loc[bmdm_mo].dropna().iloc[:, 0].values

    # create a temporary GMT file to hold the gene lists
    # only keep genes that are also in the data themselves
    timestamp = datetime.datetime.now().isoformat()
    tmp_fn = 'gene_list.%s.gmt' % timestamp
    with open(tmp_fn, 'wb') as f:
        for c in s4_hu:
            this_geneset = set(s4_hu[c].tolist()).intersection(rnaseq_dat.index)
            removed = set(s4_hu[c].tolist()).difference(rnaseq_dat.index)
            if len(removed):
                print "%d genes were removed from geneset %s as they are not found in the expression data: %s" % (
                    len(removed), c, ', '.join(list(removed))
                )

            row = [c, c] + list(this_geneset)
            f.write('\t'.join(row) + '\n')

    # run ssGSEA
    # weirdly, this requires gene symbols to be in the first column
    rr1 = pd.DataFrame(rnaseq_dat.index)
    rr2 = rnaseq_dat.copy()
    rr2.index = rr1.index
    rr = pd.concat((rr1, rr2), axis=1)

    # let's test z scores, too
    zz = rr2.subtract(rr2.mean(axis=0), axis=1)
    zz = zz.divide(rr2.std(axis=0), axis=1)
    zz = pd.concat((rr1, zz), axis=1)

    import gseapy as gp
    # ss = gp.ssgsea(rr, tmp_fn, outdir='ssGSEA_out.%s' % timestamp, permutation_num=10)  # deliberately low perm number to speed things up

    # manual ssGSEA, described by Barbie et al., Nature (2009)

    # gene set and name
    g_name = s4_hu.keys()[0]
    g = s4_hu[g_name]

    # sample and name
    s_name = rnaseq_dat.columns[0]
    s = rnaseq_dat.iloc[:, 0]

    # rank the sample
    rs = s.rank(ascending=False)

    # compute 2 x ecdf at every position (for each gene)

