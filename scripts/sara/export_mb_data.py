import pandas as pd
from rnaseq import loader, filter, general
import os
from plotting import clustering, common, pca
from sklearn.decomposition.pca import PCA
from stats import transformations
import numpy as np
from utils import output
import os
import pandas as pd
from matplotlib import pyplot as plt


if __name__ == "__main__":
    eps = 1e-2

    outdir = output.unique_output_dir("export_sb_data")
    obj_star = loader.load_by_patient(['ICb1299', '3021'], source='star', type='cell_culture', include_control=False)
    obj_salmon = loader.load_by_patient(['ICb1299', '3021'], source='salmon', type='cell_culture', include_control=False)

    # cluster plot
    tpm = filter.filter_by_cpm(obj_salmon.data, min_cpm=1, min_n_samples=4)

    batch_colours = common.COLOUR_BREWERS[len(obj_salmon.meta.batch.unique())]
    line_colours = common.COLOUR_BREWERS[2]
    cc = pd.DataFrame(line_colours[0], index=tpm.columns, columns=['Batch', 'Cell line'])

    aa, bb = obj_salmon.meta.batch.factorize()
    for i in range(aa.max()):
        cc.loc[aa == i, 'Batch'] = batch_colours[i]
    cc.loc[cc.index.str.contains('3021'), 'Cell line'] = line_colours[1]

    cg = clustering.dendrogram_with_colours(
        np.log2(tpm + eps),
        cc,
    )
    cg['fig'].savefig(os.path.join(outdir, "dendrogram_pearson_log_tpm_all_genes.png"), dpi=200)

    # pca plot
    p = PCA()
    y = p.fit_transform(np.log2(tpm + eps).transpose())
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for a, b in enumerate(bb):
        ax.scatter(
            y[aa == a, 0],
            y[aa == a, 1],
            facecolor=batch_colours[a],
            edgecolor='k',
            s=30,
            label=b
        )
    ax.set_xlabel("PC1 (%.1f %%)" % (p.explained_variance_ratio_[0] * 100))
    ax.set_ylabel("PC2 (%.1f %%)" % (p.explained_variance_ratio_[1] * 100))
    ax.legend(loc='lower right', frameon=True, facecolor='w', framealpha=0.7)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "pca_log_tpm.png"), dpi=200)

    # this looks almost the same!
    # log_tpm_qn = transformations.quantile_normalisation(np.log2(tpm + eps))
    # cg = clustering.dendrogram_with_colours(
    #     log_tpm_qn,
    #     cc,
    # )

    # export data
    cpm = obj_star.data.divide(obj_star.data.sum(), axis=1) * 1e6
    counts = obj_star.data.copy()
    tpm = obj_salmon.data.copy()
    general.add_gene_symbols_to_ensembl_data(counts, tax_id=9606)
    general.add_gene_symbols_to_ensembl_data(cpm, tax_id=9606)
    general.add_gene_symbols_to_ensembl_data(tpm, tax_id=9606)

    counts.to_excel(os.path.join(outdir, "ICb1299_CRL3021.counts.xlsx"))
    cpm.to_excel(os.path.join(outdir, "ICb1299_CRL3021.cpm.xlsx"))
    tpm.to_excel(os.path.join(outdir, "ICb1299_CRL3021.tpm.xlsx"))
    obj_star.meta.to_excel(os.path.join(outdir, "ICb1299_CRL3021.meta.xlsx"))