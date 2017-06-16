import os
from load_data import rnaseq_data
import pandas as pd
from rnaseq import filter
from matplotlib import pyplot as plt
from plotting import clustering
from stats.transformations import median_absolute_deviation, variance_stabilizing_transform
import numpy as np
from utils.output import unique_output_dir


if __name__ == "__main__":
    outdir = unique_output_dir("compare_gibco_h9")
    loader_hgic = rnaseq_data.all_hgic_loader(annotate_by='Ensembl Gene ID')
    loader_h9 = rnaseq_data.gse61794(annotate_by='Ensembl Gene ID')

    genes = loader_hgic.data.index.intersection(loader_h9.data.index)
    genes = genes[genes.str.contains('ENSG')]

    # collapse H9 replicates
    h9_data = loader_h9.data.sum(axis=1).loc[genes]
    h9_data.name = 'H9_NSC'
    h9_meta = pd.Series({
        'type': 'NSC',
        'read_count': sum(loader_h9.meta.read_count),
        'sample': 'H9_NSC',
        'disease_subgroup': 'control',
    }, name='H9_NSC')

    data = pd.concat((loader_hgic.data.loc[genes], h9_data), axis=1)
    meta = pd.concat((loader_hgic.meta, pd.DataFrame(h9_meta).transpose()), axis=0)

    # filter
    data = filter.filter_by_cpm(data, min_cpm=1, min_n_samples=3, unless_cpm_gt=10)

    # normalise by read count
    cpm = (data + 1).divide((data + 1).sum(axis=0), axis=1) *  1e6

    # transform
    log_data = np.log2(cpm)
    vst_data = variance_stabilizing_transform(cpm)
    mad_log_srt = median_absolute_deviation(log_data).sort_values(ascending=False)
    mad_vst_srt = median_absolute_deviation(vst_data).sort_values(ascending=False)

    for NGENE in [500, 1000, 1500, 2000, 2500]:

        cg = clustering.plot_correlation_clustermap(log_data.loc[mad_log_srt.index[:NGENE]], vmin=0., vmax=1., metric='correlation')
        cg.gs.update(bottom=0.3, right=0.7)
        cg.fig.savefig(os.path.join(outdir, "gbm_nsc_correlation_clustermap_logtransform_top%d.png" % NGENE), dpi=200)
        cg.fig.savefig(os.path.join(outdir, "gbm_nsc_correlation_clustermap_logtransform_top%d.pdf" % NGENE))

        cg = clustering.plot_correlation_clustermap(vst_data.loc[mad_vst_srt.index[:NGENE]], vmin=0., vmax=1., metric='correlation')
        cg.gs.update(bottom=0.3, right=0.7)
        cg.fig.savefig(os.path.join(outdir, "gbm_nsc_correlation_clustermap_vsttransform_top%d.png" % NGENE), dpi=200)
        cg.fig.savefig(os.path.join(outdir, "gbm_nsc_correlation_clustermap_vsttransform_top%d.pdf" % NGENE))

    log_nsc_data = log_data.loc[:, log_data.columns.str.contains('NSC')]
    vst_nsc_data = vst_data.loc[:, log_data.columns.str.contains('NSC')]
    mad_log_nsc_srt = median_absolute_deviation(log_nsc_data).sort_values(ascending=False)
    mad_vst_nsc_srt = median_absolute_deviation(vst_nsc_data).sort_values(ascending=False)

    for NGENE in [500, 1000, 1500, 2000, 2500]:
        cg = clustering.plot_correlation_clustermap(log_nsc_data.loc[mad_log_nsc_srt.index[:NGENE]], vmin=0., vmax=1., metric='correlation')
        cg.gs.update(bottom=0.3, right=0.7)
        cg.fig.savefig(os.path.join(outdir, "nsc_correlation_clustermap_logtransform_top%d.png" % NGENE), dpi=200)
        cg.fig.savefig(os.path.join(outdir, "nsc_correlation_clustermap_logtransform_top%d.pdf" % NGENE))
        cg = clustering.plot_correlation_clustermap(vst_nsc_data.loc[mad_vst_nsc_srt.index[:NGENE]], vmin=0., vmax=1., metric='correlation')
        cg.gs.update(bottom=0.3, right=0.7)
        cg.fig.savefig(os.path.join(outdir, "nsc_correlation_clustermap_vsttransform_top%d.png" % NGENE), dpi=200)
        cg.fig.savefig(os.path.join(outdir, "nsc_correlation_clustermap_vsttransform_top%d.pdf" % NGENE))

        cg = clustering.plot_clustermap(
            log_nsc_data.loc[mad_log_nsc_srt.index[:NGENE]],
            show_gene_labels=False,
            rotate_xticklabels=True,
            cmap='RdBu_r',
            metric='correlation',
        )
        cg.gs.update(bottom=0.2)
        cg.fig.savefig(os.path.join(outdir, "nsc_expression_clustermap_logtransform_top%d.png" % NGENE), dpi=200)
        cg.fig.savefig(os.path.join(outdir, "nsc_expression_clustermap_logtransform_top%d.pdf" % NGENE))
        cg = clustering.plot_clustermap(
            vst_nsc_data.loc[mad_vst_nsc_srt.index[:NGENE]],
            show_gene_labels=False,
            rotate_xticklabels=True,
            cmap='RdBu_r',
            metric='correlation',
        )
        cg.gs.update(bottom=0.2)
        cg.fig.savefig(os.path.join(outdir, "nsc_expression_clustermap_vsttransform_top%d.png" % NGENE), dpi=200)
        cg.fig.savefig(os.path.join(outdir, "nsc_expression_clustermap_vsttransform_top%d.pdf" % NGENE))