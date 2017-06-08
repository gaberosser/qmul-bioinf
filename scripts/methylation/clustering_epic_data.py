from load_data import methylation_array
from plotting import clustering, pca
from sklearn.decomposition import PCA


if __name__ == "__main__":
    data, meta = methylation_array.hgic_methylationepic(norm_method='swan')
    # Don't know why, but some probes (~2000) are only present in one OR the other sample
    # Therefore, remove those
    data = data.dropna()

    # add some extra meta information
    meta.loc[:, 'cell_type'] = 'NSC'
    meta.loc[meta.index.str.contains('GBM'), 'cell_type'] = 'GBM'
    meta.loc[:, 'subgroup'] = 'RTK I'
    meta.loc[meta.index.str.contains('024'), 'subgroup'] = 'Unknown'
    meta.loc[meta.index.str.contains('026'), 'subgroup'] = 'Unknown'
    meta.loc[meta.index.str.contains('044'), 'subgroup'] = 'Mesenchymal'
    meta.loc[meta.index.str.contains('GIBCO'), 'subgroup'] = 'NSC'

    p = PCA(n_components=3)
    p.fit(data.transpose())
    y = p.transform(data.transpose())

    pca.pca_plot_by_group_3d(y, meta.cell_type, plot_ellipsoids=False)

    s = data.std(axis=1).sort_values(ascending=False)

    # plot a correlation clustermap of GBM without 024 (which is a definite outlier)
    clustering.plot_correlation_clustermap(
        data.loc[s.index[:8000], (meta.cell_type == 'GBM') & (~meta.index.str.contains('024'))]
    )
