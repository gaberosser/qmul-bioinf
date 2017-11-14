from rnaseq import differential_expression
from load_data import rnaseq_data
from scripts.integrate_rnaseq_methylation import rtkii_integrated_analysis
from plotting import venn
import numpy as np
import os
from utils import output, setops
from matplotlib import pyplot as plt
import seaborn as sns


if __name__ == "__main__":
    # all n=2 samples
    pids = ['018', '044', '049', '050', '052', '054', '061']
    outdir = output.unique_output_dir("compare_edger_models")

    obj = rnaseq_data.load_by_patient(pids, annotate_by='Ensembl Gene ID')
    obj.data = obj.data.loc[obj.data.index.str.contains('ENSG')]
    de_qlglm = rtkii_integrated_analysis.compute_de(obj, pids, method='QLGLM')
    de_glm = rtkii_integrated_analysis.compute_de(obj, pids, method='GLM')

    figwidth = 3 * len(pids) ** 0.9
    set_labels = ('Isogenic', 'Reference')
    for res, lbl in [(de_glm, 'GLM'), (de_qlglm, 'QL GLM')]:
        fig, axs = plt.subplots(ncols=len(pids), nrows=3, figsize=[figwidth, 8])
        for i, pid in enumerate(pids):

            _, venn_sets, venn_counts = venn.venn_diagram(
                res['de_matched'][pid].index,
                res['de_gibco'][pid].index,
                set_labels=set_labels,
                ax=axs[0, i]
            )
            sign_m = np.sign(res['de_matched'][pid].loc[venn_sets['11'], 'logFC'])
            sign_r = np.sign(res['de_gibco'][pid].loc[venn_sets['11'], 'logFC'])
            print "GBM%s, model %s: %d genes have discordant FC direction." % (pid, lbl, (sign_m != sign_r).sum())

            up_idx_m = res['de_matched'][pid].logFC > 0
            up_idx_r = res['de_gibco'][pid].logFC > 0
            venn.venn_diagram(
                res['de_matched'][pid].index[up_idx_m],
                res['de_gibco'][pid].index[up_idx_r],
                set_labels=set_labels,
                ax=axs[1, i]
            )
            down_idx_m = res['de_matched'][pid].logFC < 0
            down_idx_r = res['de_gibco'][pid].logFC < 0
            venn.venn_diagram(
                res['de_matched'][pid].index[down_idx_m],
                res['de_gibco'][pid].index[down_idx_r],
                set_labels=set_labels,
                ax=axs[2, i]
            )
            axs[0, i].set_title("GBM%s all" % pid)
            axs[1, i].set_title("Up in GBM%s vs NSC" % pid)
            axs[2, i].set_title("Down in GBM%s vs NSC" % pid)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "de_gene_counts_%s.png" % lbl), dpi=200)
        fig.savefig(os.path.join(outdir, "de_gene_counts_%s.pdf" % lbl))

    # look at the overlap between methods (choose matched comparison for this)
    set_labels = ['LR test', 'QL F test']
    fig, axs = plt.subplots(ncols=len(pids), nrows=1, figsize=[figwidth, 3.5])
    for i, pid in enumerate(pids):
        _, venn_sets, venn_counts = venn.venn_diagram(
            de_glm['de_matched'][pid].index,
            de_qlglm['de_matched'][pid].index,
            set_labels=set_labels,
            ax=axs[i]
        )
        axs[i].set_title("GBM%s" % pid)
        sign_lr = np.sign(de_glm['de_matched'][pid].loc[venn_sets['11'], 'logFC'])
        sign_ql = np.sign(de_qlglm['de_matched'][pid].loc[venn_sets['11'], 'logFC'])
        print "GBM%s: %d genes have discordant FC direction." % (pid, (sign_lr != sign_ql).sum())
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "de_gene_counts_different_tests.png"), dpi=200)
    fig.savefig(os.path.join(outdir, "de_gene_counts_different_tests.pdf"))
