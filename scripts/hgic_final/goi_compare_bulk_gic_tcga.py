"""
For specified gene(s) of interest, generate a plot comparing the expression levels observed in GIC, FFPE and
across the GBM TCGA cohort
"""

import os
from decimal import Decimal
from settings import DATA_DIR
from utils import output
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from rnaseq import loader
from scripts.hgic_final import consts
from hgic_consts import NH_ID_TO_PATIENT_ID_MAP
from utils import reference_genomes


if __name__ == '__main__':
    outdir = output.unique_output_dir()
    eps = 0.1
    seps = Decimal('0.1')
    gois = ['PTGER4']
    ens_gois = reference_genomes.gene_symbol_to_ensembl(gois)

    # pids = consts.PIDS
    pids = ['018']

    meta_fn = os.path.join(DATA_DIR, 'rnaseq', 'tcga_gbm', 'primary_tumour/htseq-count_fpkm/sources.csv')
    dat_fn = os.path.join(DATA_DIR, 'rnaseq', 'tcga_gbm', 'primary_tumour/htseq-count_fpkm/fpkm.csv')
    tcga_meta = pd.read_csv(meta_fn, header=0, index_col=0)
    tcga_dat = pd.read_csv(dat_fn, header=0, index_col=0)

    # filter: primary GBM only
    ix = (tcga_meta.idh1_status == 'WT')
    tcga_meta = tcga_meta[ix]
    tcga_dat = tcga_dat[tcga_meta.index]
    tcga_dat = tcga_dat.divide(tcga_dat.sum(axis=0), axis=1) * 1e6

    # load our data
    gic_obj = loader.load_by_patient(consts.PIDS, source='salmon', include_control=False, type='cell_culture')
    ffpe_obj = loader.load_by_patient(consts.PIDS, source='salmon', include_control=False, type='ffpe')

    # add NH ID and patient ID to FFPE
    nh_id = ffpe_obj.meta.index.str.replace(r'(_?)(DEF|SP).*', '')
    p_id = [NH_ID_TO_PATIENT_ID_MAP[t.replace('_', '-')] for t in nh_id]
    ffpe_obj.meta.insert(0, 'nh_id', nh_id)
    ffpe_obj.meta.insert(0, 'patient_id', p_id)

    # ditto GIC
    gic_obj.meta.insert(
        0,
        'patient_id',
        gic_obj.meta.index.str.replace(r'GBM(?P<pid>[0-9]*)_.*', r'\g<pid>')
    )

    # filter: GIC only
    ix = gic_obj.meta['type'] == 'GBM'
    gic_obj.filter_samples(ix)

    # filter: only best FFPE samples
    ffpe_obj.filter_by_sample_name(consts.FFPE_RNASEQ_SAMPLES, exact=True)

    ## PLOT
    for g in gois:
        e = ens_gois[g]
        fig, ax = plt.subplots(1)
        the_tcga_dat = np.log10(tcga_dat.loc[e].squeeze() + eps)
        the_gic_dat = gic_obj.data.loc[e, gic_obj.meta.patient_id.isin(pids)]
        the_gic_dat = np.log10(the_gic_dat + eps)
        the_gic_dat.index = the_gic_dat.index.str.replace('GBM', 'GIC').str.replace(r'_(?P<a>P[0-9]*)', r' (\g<a>)')
        the_ffpe_dat = ffpe_obj.data.loc[e, ffpe_obj.meta.patient_id.isin(pids)]
        the_ffpe_dat = np.log10(the_ffpe_dat + eps)
        ax.hist(the_tcga_dat, 20, label='TCGA')
        for col, t in the_gic_dat.items():
            ax.axvline(t, c='k')
            ax.text(t, 20, col, rotation=90, ha='right', va='top', color='k')
            # ax.scatter([t], [1], label=col)
        for col, t in the_ffpe_dat.items():
            # ax.scatter([t], [1], label='FFPE')
            ax.axvline(t, c='0.3')
            ax.text(t, 20, "GBM%s (bulk)" % ffpe_obj.meta.patient_id[col], rotation=90, ha='right', va='top', color='0.3')
        ax.set_title(g)
        ax.set_xlabel(r'$\log_{10}(\mathrm{TPM} + %s)$' % seps)
        ax.set_ylabel('Frequency')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, '%s_tcga_gic_ffpe.png' % g), dpi=200)
        fig.savefig(os.path.join(outdir, '%s_tcga_gic_ffpe.tiff' % g), dpi=200)
        fig.savefig(os.path.join(outdir, '%s_tcga_gic_ffpe.pdf' % g))



