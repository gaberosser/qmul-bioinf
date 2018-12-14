import pandas as pd
import os
import re
from settings import HGIC_LOCAL_DIR, GIT_LFS_DATA_DIR
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import clustering, common
from utils import output
from scipy import stats
from scipy.cluster import hierarchy as hc
import numpy as np
from utils import setops, excel
from scripts.hgic_final import consts
from stats import nht
import multiprocessing as mp


XCELL_SIGNATURE_FN = os.path.join(GIT_LFS_DATA_DIR, 'xcell', 'ESM3_signatures.xlsx')
IPA_PATHWAY_DIR = os.path.join(
    HGIC_LOCAL_DIR,
    'current/core_pipeline/rnaseq/merged_s1_s2/ipa/pathways'
)
IPA_PATHWAY_FN = os.path.join(
    IPA_PATHWAY_DIR,
    'ipa_results_significant.xlsx'
)


if __name__ == '__main__':

    outdir = output.unique_output_dir()
    pids = consts.PIDS

    # load syngeneic deregulated genes from raw reports
    file_patt = 'de_s2_{pid}_{cmp}.txt'
    genes_dereg = {}
    genes_dereg_intersect = {}
    genes_dereg_union = {}
    z_scores = {}
    logp = {}

    for pid in pids:
        genes_dereg[pid] = {}
        fn = os.path.join(IPA_PATHWAY_DIR, file_patt.format(pid=pid, cmp='syngeneic'))
        this = pd.read_csv(fn, sep='\t', skiprows=2, header=0, index_col=0)
        this.columns = ['-logp', 'ratio', 'z', 'genes']
        this.index = [x.decode('utf-8') for x in this.index]
        the_genes = this.genes.apply(lambda x: x.split(','))
        z_scores[pid] = this.z
        logp[pid] = this['-logp']
        for pw, arr in the_genes.items():
            genes_dereg[pid][pw] = ','.join(arr)
            if pw in genes_dereg_intersect:
                genes_dereg_intersect[pw] = set(genes_dereg_intersect[pw]).intersection(arr)
            else:
                genes_dereg_intersect[pw] = set(arr)

            if pw in genes_dereg_union:
                genes_dereg_union[pw] = set(genes_dereg_union[pw]).union(arr)
            else:
                genes_dereg_union[pw] = set(arr)


    for k, v in genes_dereg_intersect.items():
        genes_dereg_intersect[k] = ','.join(v)

    for k, v in genes_dereg_union.items():
        genes_dereg_union[k] = ','.join(v)

    # intersection
    gdi = pd.Series(genes_dereg_intersect)
    gdi.to_excel(os.path.join(outdir, 'ipa_deregulated_genes_intersection_across_patients.xlsx'))

    # union
    gdu = pd.Series(genes_dereg_union)
    gdu.to_excel(os.path.join(outdir, 'ipa_deregulated_genes_union_across_patients.xlsx'))

    # also export full lists
    gd_all = pd.DataFrame(genes_dereg)
    gd_all.to_excel(os.path.join(outdir, 'ipa_deregulated_genes_patient_specific.xlsx'))

    # z scores
    z_scores = pd.DataFrame(z_scores).dropna(axis=0, how='all')
    z_scores.to_excel(os.path.join(outdir, 'ipa_pathways_z_scores.xlsx'))

    # p values
    p_values = pd.DataFrame(logp).dropna(axis=0, how='all')
    p_values.to_excel(os.path.join(outdir, 'ipa_pathways_-log10_p.xlsx'))
