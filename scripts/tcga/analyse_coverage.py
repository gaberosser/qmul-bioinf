import pandas as pd
from utils import output
from settings import OUTPUT_DIR, DATA_DIR
import os
from plotting import venn, clustering
from matplotlib import pyplot as plt


if __name__ == "__main__":
    outdir = output.unique_output_dir("tcga_gbm_analysis", reuse_empty=True)
    # load meta files
    meta_fn = {
        'rnaseq': os.path.join(DATA_DIR, 'rnaseq', 'tcga_gbm', 'primary_tumour', 'rnaseq.meta.csv'),
        'marr_u133': os.path.join(DATA_DIR, 'microarray', 'tcga_gbm', 'primary_tumour', 'microarray.meta.ht_hg_u133a.csv'),
        'marr_agilent1': os.path.join(DATA_DIR, 'microarray', 'tcga_gbm', 'primary_tumour', 'microarray.meta.agilentg4502a_07_1.csv'),
        'marr_agilent2': os.path.join(DATA_DIR, 'microarray', 'tcga_gbm', 'primary_tumour', 'microarray.meta.agilentg4502a_07_2.csv'),
        'meth_450k': os.path.join(DATA_DIR, 'methylation', 'tcga_gbm', 'primary_tumour', 'methylation.450k.meta.csv'),
        'meth_27k': os.path.join(DATA_DIR, 'methylation', 'tcga_gbm', 'primary_tumour', 'methylation.27k.meta.csv'),
    }

    meta = {}
    for k, fn in meta_fn.items():
        meta[k] = pd.read_csv(fn, header=0, index_col=0)

    # microarray only
    fig = plt.figure()
    ax = fig.add_subplot(111)
    venn.venn_diagram(
        meta['marr_u133'].case_id,
        meta['marr_agilent1'].case_id,
        meta['marr_agilent2'].case_id,
        set_labels=('U133', 'Agilent 1', 'Agilent 2'),
        ax=ax
    )
    plt.tight_layout(rect=(0, 0, 1., 1.05))
    plt.savefig(os.path.join(outdir, "venn_by_caseid_microarray.png"), dpi=200)

    # rnaseq / methylation
    fig = plt.figure()
    ax = fig.add_subplot(111)
    venn.venn_diagram(
        meta['rnaseq'].index,
        meta['meth_450k'].index,
        meta['meth_27k'].index,
        set_labels=('RNA-Seq', 'Methylation 450K', 'Methylation 27K'),
        ax=ax
    )
    plt.tight_layout(rect=(0, 0, 1.05, 1.05))
    plt.savefig(os.path.join(outdir, "venn_rnaseq_methylation.png"), dpi=200)

    # check that matching by case ID makes no difference (same case, different portion)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    venn.venn_diagram(
        meta['rnaseq'].case_id,
        meta['meth_450k'].case_id,
        meta['meth_27k'].case_id,
        set_labels=('RNA-Seq', 'Methylation 450K', 'Methylation 27K'),
        ax=ax
    )
    plt.tight_layout(rect=(0, 0, 1., 1.05))
    plt.savefig(os.path.join(outdir, "venn_by_caseid_rnaseq_methylation.png"), dpi=200)

    # microarray / methylation
    fig = plt.figure()
    ax = fig.add_subplot(111)
    venn.venn_diagram(
        meta['marr_u133'].case_id,
        meta['meth_450k'].case_id,
        meta['meth_27k'].case_id,
        set_labels=('U133', 'Methylation 450K', 'Methylation 27K'),
        ax=ax
    )
    ax.set_facecolor('w')
    ax.figure.tight_layout(rect=(0, 0, 1.05, 1.05))
    ax.figure.savefig(os.path.join(outdir, "venn_by_caseid_u133_methylation.png"), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    venn.venn_diagram(
        meta['meth_450k'].index,
        meta['marr_u133'].index,
        meta['marr_agilent1'].index,
        meta['marr_agilent2'].index,
        set_labels=('Methylation 450K', 'U133', 'Agilent 1', 'Agilent 2'),
        ax=ax
    )
    ax.set_facecolor('w')
    ax.figure.tight_layout(rect=(0, 0, 1.05, 1.05))
    ax.figure.savefig(os.path.join(outdir, "venn_microarray_methylation450.png"), dpi=200)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    venn.venn_diagram(
        meta['meth_27k'].index,
        meta['marr_u133'].index,
        meta['marr_agilent1'].index,
        meta['marr_agilent2'].index,
        set_labels=('Methylation 27K', 'U133', 'Agilent 1', 'Agilent 2'),
        ax=ax
    )
    ax.set_facecolor('w')
    ax.figure.tight_layout(rect=(0, 0, 1.05, 1.05))
    ax.figure.savefig(os.path.join(outdir, "venn_microarray_methylation27.png"), dpi=200)

    # just the U133, RNA-Seq and methylation
    fig = plt.figure()
    ax = fig.add_subplot(111)
    venn.venn_diagram(
        meta['rnaseq'].index,
        meta['meth_450k'].index,
        meta['meth_27k'].index,
        meta['marr_u133'].index,
        set_labels=('RNA-Seq', 'Methylation 450K', 'Methylation 27K', 'U133'),
        ax=ax
    )
    ax.set_facecolor('w')
    ax.figure.tight_layout(rect=(0, 0, 1.05, 1.05))
    ax.figure.savefig(os.path.join(outdir, "venn_rna_microarray_methylation.png"), dpi=200)

    # by Case ID
    fig = plt.figure()
    ax = fig.add_subplot(111)
    venn.venn_diagram(
        meta['rnaseq'].case_id,
        meta['meth_450k'].case_id,
        meta['meth_27k'].case_id,
        meta['marr_u133'].case_id,
        set_labels=('RNA-Seq', 'Methylation 450K', 'Methylation 27K', 'U133'),
        ax=ax
    )
    ax.set_facecolor('w')
    ax.figure.tight_layout(rect=(0, 0, 1.05, 1.05))
    ax.figure.savefig(os.path.join(outdir, "venn_by_caseid_rna_microarray_methylation.png"), dpi=200)

    # methylation classification where available
    m = meta['meth_450k'].loc[~meta['meth_450k'].match.isnull()]
    counts = m.groupby('match').size()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(counts.size), counts.values)
    ax.xaxis.set_ticks(range(counts.size))
    ax.xaxis.set_ticklabels(counts.index, rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "sturm_classification_frequency.png"), dpi=200)

