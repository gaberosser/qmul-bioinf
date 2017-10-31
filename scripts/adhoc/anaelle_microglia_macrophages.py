import os
import pandas as pd
from matplotlib import pyplot as plt
from plotting import venn
from settings import GIT_LFS_DATA_DIR
from utils.output import unique_output_dir


if __name__ == '__main__':
    outdir = unique_output_dir("mg_bmdm_venn")
    indir = os.path.join(GIT_LFS_DATA_DIR, 'GSE86573_bowman_de')

    gl261_mg = pd.read_csv(os.path.join(indir, 'gl261_mg_vs_healthy_mg.csv'), header=0, index_col=0)
    gl261_bmdm = pd.read_csv(os.path.join(indir, 'gl261_bmdm_vs_healthy_monocyte.csv'), header=0, index_col=0)
    gemm_mg = pd.read_csv(os.path.join(indir, 'gemm_mg_vs_healthy_mg.csv'), header=0, index_col=0)
    gemm_bmdm = pd.read_csv(os.path.join(indir, 'gemm_bmdm_vs_healthy_monocyte.csv'), header=0, index_col=0)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    v, sets, counts = venn.venn_diagram(gl261_mg.index, gemm_mg.index, set_labels=("GL261 MG", "GEMM MG"), ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "gl261_mg-gemm_mg.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "gl261_mg-gemm_mg.pdf"))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    v, sets, counts = venn.venn_diagram(gl261_bmdm.index, gemm_bmdm.index, set_labels=("GL261 BMDM", "GEMM BMDM"), ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "gl261_bmdm-gemm_bmdm.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "gl261_bmdm-gemm_bmdm.pdf"))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    v, sets, counts = venn.venn_diagram(gl261_mg.index, gl261_bmdm.index, set_labels=("GL261 MG", "GL261 BMDM"), ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "gl261_mg-gl261_bmdm.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "gl261_mg-gl261_bmdm.pdf"))