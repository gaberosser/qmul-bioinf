import pandas as pd
import numpy as np
from methylation import dmr, annotation_gene_to_ensembl


## TODO: test this

def compute_joint_de_dmr(
        dmr_results,
        de_results,
        de_gene_column='Gene Symbol',
        dmr_include='significant',
):
    """

    :param dmr_results: Dictionary containing DmrResults objects. Keyed by pid.
    :param de_results: Dictionary keyed by PID giving pd.DataFrame instances containing results.
    :param de_gene_column:
    :param dmr_include: Which results to include from the DMR results
    :return:
    """
    # quick check that the two results dictionaries have matching keys
    if set(dmr_results.keys()) != set(de_results.keys()):
        raise ValueError("Keys in the dmr_results and de_results dictionaries do not match.")

    # DMR and DE column exclusion list
    dmr_col_exclusion = [
        'cluster_id', 'gene', 'median_1', 'median_2'
    ]
    de_col_exclusion = [
        'Gene Symbol', 'unshrunk.logFC', 'logCPM', 'PValue', 'Direction'
    ]

    res = {}

    for sid in dmr_results:
        res[sid] = {}

        this_dmr = dmr_results[sid].to_table(include=dmr_include, expand=True)
        this_de = de_results[sid]

        # use apply to get the matching DE lines for each DMR entry to filter
        ens_ix = annotation_gene_to_ensembl.gene_to_ens(this_dmr.gene)
        ens_ix.index = this_dmr.index

        de_match = this_de.reindex(ens_ix)
        de_match.index = ens_ix.index

        dmr_match = this_dmr.loc[~de_match.isnull().all(axis=1)].drop(dmr_col_exclusion, axis=1)
        de_match = de_match.dropna(axis=0, how='all').drop(de_col_exclusion, axis=1)
        ens_ix = ens_ix[dmr_match.index]

        dmr_match.columns = ["dmr_%s" % t for t in dmr_match.columns]
        de_match.columns = ["de_%s" % t for t in de_match.columns]

        dmr_match.insert(
            dmr_match.shape[1],
            'dmr_direction',
            ['Hyper' if t > 0 else 'Hypo' for t in dmr_match['dmr_median_delta']]
        )
        de_match.insert(
            de_match.shape[1],
            'de_direction',
            ['U' if t > 0 else 'D' for t in de_match['de_logFC']]
        )

        matched = pd.concat((de_match, dmr_match), axis=1)
        matched.insert(0, 'gene', [t[1] for t in matched.index])
        matched.insert(0, 'cluster_id', [t[0] for t in matched.index])

        res[sid] = matched

    return res