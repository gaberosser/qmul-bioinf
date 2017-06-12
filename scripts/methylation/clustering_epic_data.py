from load_data import methylation_array
from plotting import clustering, pca
from clustering import lda
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np


if __name__ == "__main__":
    REF_META_SUBGRP_LABEL = 'dna methylation subgroup'

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

    ref_data, ref_meta = methylation_array.gse36278()

    # only keep intersecting probes
    probe_idx = ref_data.index.intersection(data.index)

    sample_names = np.random.permutation(ref_data.columns)
    n_train = 120
    n_test = ref_data.shape[1] - n_train

    ref_train = ref_data.loc[probe_idx, sample_names[:n_train]]
    ref_test = ref_data.loc[probe_idx, sample_names[n_train:]]
    ref_meta_train = ref_meta.loc[sample_names[:n_train]]
    ref_meta_test = ref_meta.loc[sample_names[n_train:]]

    grp_ids, grp_names = ref_meta.loc[sample_names, REF_META_SUBGRP_LABEL].factorize()
    grp_ids_train = grp_ids[:n_train]
    grp_ids_test = grp_ids[n_train:]
    n_grp = len(np.unique(grp_ids))

    # construct a LDA classifier
    # my simple implementation
    deltas = np.concatenate((
        np.linspace(0., 3., 30),
        np.arange(3.5, 10.5, 0.5)
    ))

    ld = lda.NearestCentroidClassifier(ref_train.loc[probe_idx], ref_meta.loc[sample_names, REF_META_SUBGRP_LABEL])

    nerr_train, nerr_test, clas = lda.run_validation(
        deltas,
        ref_train.loc[probe_idx],
        ref_meta_train.loc[:, REF_META_SUBGRP_LABEL],
        test_data=ref_test.loc[probe_idx],
        test_labels=ref_meta_test.loc[:, REF_META_SUBGRP_LABEL],
        flat_prior=True
    )

    # shrinkage consumes too much memory
    l = LinearDiscriminantAnalysis()
    l = l.fit(ref_train.loc[probe_idx].transpose(), grp_ids_train)
    score_train = l.score(ref_train.loc[probe_idx].transpose(), grp_ids_train)
    score_test = l.score(ref_test.loc[probe_idx].transpose(), grp_ids_test)

    # can reduce dimensionality like this:
    red_ref = l.transform(ref_data.loc[probe_idx].transpose())  # 142 x 6
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cmap = plt.cm.get_cmap('Set1')

    for i, the_grp in enumerate(np.unique(grp_ids[grp_ids != -1])):
        this_idx = np.where(grp_ids == the_grp)[0]
        this_c = cmap.colors[i]
        this_lbl = grp_names[i]
        ax.scatter(red_ref[this_idx, 0], red_ref[this_idx, 1], c=this_c, label=this_lbl)
    this_idx = np.where(grp_ids == -1)[0]
    this_c = 'gray'
    this_lbl = 'Unknown'
    ax.scatter(red_ref[this_idx, 0], red_ref[this_idx, 1], c=this_c, label=this_lbl)
    ax.legend()

    rf = RandomForestClassifier(n_estimators=1000)
    rf = rf.fit(ref_data.loc[probe_idx].transpose(), grp_ids)
