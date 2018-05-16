import pandas as pd
from methylation import loader
import numpy as np
import os
import collections
from utils import output, setops
from methylation import dmr, process
from settings import DATA_DIR_NON_GIT
import multiprocessing as mp


def compute_median_betas_one_sample(the_dat, probes_by_gene):
    missing_probes = set()
    missing_genes = []
    res = {}
    for g, probes in probes_by_gene.items():
        try:
            res[g] = the_dat.loc[probes].median(axis=0)
        except KeyError:
            missing_probes.update(probes)
            missing_genes.append(g)
    return res, missing_genes, missing_probes


if __name__ == "__main__":
    outdir = output.unique_output_dir("report_beta_values")

    ######################
    # 1: PRIMARY TUMOUR  #
    ######################

    # in this case, we want the median beta value over all probes that are associated with a given gene
    # we'll exclude those associated with gene body only
    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'tcga_gbm', 'primary_tumour')
    meta_fn = os.path.join(indir, 'methylation.450k.meta.csv')
    dat_fn = os.path.join(indir, 'methylation.450k.csv.gz')

    meta = pd.read_csv(meta_fn, header=0, index_col=0)
    dat = pd.read_csv(dat_fn, header=0, index_col=0, skiprows=[1])

    print "Primary tumour (%d samples)" % meta.shape[0]

    anno = loader.load_illumina_methylation450_annotation()

    # reduce anno to (probe ID, gene, relation)
    probe_tups = set()
    for i, row in anno.iterrows():
        if pd.isnull(row.UCSC_RefGene_Name):
            continue
        genes = row.UCSC_RefGene_Name.split(';')
        rels = row.UCSC_RefGene_Group.split(';')
        for g, r in zip(genes, rels):
            probe_tups.add(
                (i, g, r)
            )

    probe_tups = list(probe_tups)

    # only keep those that are TSS-related

    tss_probe_tups = [t for t in probe_tups if 'TSS' in t[2]]
    probes_by_gene = {}
    for i, g, _ in tss_probe_tups:
        probes_by_gene.setdefault(g, []).append(i)


    # extract beta values

    # this is VERY slow because of the number of samples:

    median_beta = pd.DataFrame(index=probes_by_gene.keys(), columns=meta.index)
    missing_probes = set()
    missing_genes = []


    for g, probes in probes_by_gene.items():
        try:
            median_beta.loc[g] = dat.loc[probes].median(axis=0)
        except KeyError:
            # print "Gene %s: none of the %d probes are found in the data" % (g, len(probes))
            missing_probes.update(probes)
            missing_genes.append(g)
    median_beta = median_beta.dropna()

    # solution (?): split into samples and run multiprocessing
    # pool = mp.Pool()
    # jobs = {}
    #
    # for sn in dat.columns:
    #     the_dat = dat[sn]
    #     jobs[sn] = pool.apply_async(compute_median_betas_one_sample, args=(the_dat, probes_by_gene))
    #
    # pool.close()
    # pool.join()
    #
    # median_beta = pd.DataFrame(index=probes_by_gene.keys(), columns=meta.index)
    # missing_genes = {}
    # missings_probes = {}
    #
    # for sn in jobs:
    #     this_res, missing_genes[sn], missings_probes[sn] = jobs[sn].get(10)
    #     median_beta.loc[this_res.keys(), sn] = this_res.values()

    # median_beta = median_beta.dropna()

    print "Skipped %d genes due to %d missing probes." % (len(missing_genes), len(missing_probes))
    median_beta.to_excel(os.path.join(outdir, 'tcga450K_primary_tumour_median_beta_by_gene_tss_only.xlsx'))

    ######################
    # 2: SOLID NORMAL    #
    ######################

    indir = os.path.join(DATA_DIR_NON_GIT, 'methylation', 'tcga_gbm', 'solid_tissue_normal')
    meta_fn = os.path.join(indir, 'methylation_normal.450k.meta.csv')
    dat_fn = os.path.join(indir, 'methylation_normal.450k.csv.gz')

    meta = pd.read_csv(meta_fn, header=0, index_col=0)
    dat = pd.read_csv(dat_fn, header=0, index_col=0, skiprows=[1])

    print "Solid normal tissue (%d samples)" % meta.shape[0]

    median_beta_normal = pd.DataFrame(index=probes_by_gene.keys(), columns=meta.index)
    missing_probes = set()
    missing_genes = []


    for g, probes in probes_by_gene.items():
        try:
            median_beta_normal.loc[g] = dat.loc[probes].median(axis=0)
        except KeyError:
            # print "Gene %s: none of the %d probes are found in the data" % (g, len(probes))
            missing_probes.update(probes)
            missing_genes.append(g)
    median_beta_normal = median_beta_normal.dropna()

    print "Skipped %d genes due to %d missing probes." % (len(missing_genes), len(missing_probes))
    median_beta_normal.to_excel(os.path.join(outdir, 'tcga450K_solid_tissue_normal_median_beta_by_gene_tss_only.xlsx'))