import pandas as pd
import numpy as np
from methylation import dmr


def compute_joint_de_dmr(
        dmr_results,
        de_results,
        de_gene_column='Gene Symbol',
):
    """

    :param dmr_results: Dictionary containing DmrResults objects. Keyed by pid.
    :param de_results: Dictionary keyed by PID giving pd.DataFrame instances containing results.
    :param de_gene_column:
    :return:
    """

    # DE map: (name in combined_list, name in original DMR list, datatype)
    de_map = [
        ('gene', 'Gene Symbol', 'object'),
        ('de_logfc', 'logFC', 'float'),
        ('de_logcpm', 'logCPM', 'float'),
        ('de_padj', 'FDR', 'float'),
    ]
    de_cols = [t[0] for t in de_map]
    de_col_dtypes = pd.Series([t[2] for t in de_map], index=de_cols)

    res = {}

    for sid in dmr_results:
        res[sid] = {}

        this_classes = sorted(dmr_results[sid].classes)
        this_dmr = dmr_results[sid].to_table(include='significant')
        this_de = de_results[sid]

        # use apply to get the matching DE lines for each DMR entry to filter
        m = this_dmr.genes.apply(lambda x: this_de.loc[this_de.loc[:, de_gene_column].isin(x)])
        n_match = m.apply(lambda x: x.shape[0])
        this_dmr = this_dmr.loc[n_match > 0]
        m = m.loc[n_match > 0]

        # DMR map: (name in combined_list, name in original DMR list, datatype)
        dmr_map = [
            ('dmr_chr', 'chr', 'object'),
            ('dmr_median1', 'median_1', 'float'),
            ('dmr_median2', 'median_2', 'float'),
            ('dmr_median_delta', 'median_delta', 'float'),
            ('dmr_padj', 'padj', 'float'),
        ] + [('dmr_class_%s' % t, 'class_%s' % t, 'bool') for t in this_classes]
        dmr_cols = ['dmr_cid'] + [t[0] for t in dmr_map]
        dmr_col_dtypes = pd.Series(
            ['int'] + [t[2] for t in dmr_map],
            index=dmr_cols
        )

        this_df = pd.DataFrame(columns=de_cols + dmr_cols)
        this_df = this_df.astype(pd.concat((de_col_dtypes, dmr_col_dtypes)).to_dict())

        for cluster_id, row in this_dmr.iterrows():
            # matching entry in DE (by gene name)
            # select and change the DE column names
            de_match = m.loc[cluster_id].loc[:, [t[1] for t in de_map]]
            de_match.columns = de_cols
            n = de_match.shape[0]

            # form the DMR data block by repeating the same row
            row_for_rpt = [cluster_id] + row.loc[[t[1] for t in dmr_map]].tolist()
            dmr_match = pd.DataFrame([row_for_rpt] * n, columns=dmr_cols, index=de_match.index).astype(dmr_col_dtypes.to_dict())

            this_match = pd.concat((de_match, dmr_match), axis=1)
            this_df = this_df.append(this_match, ignore_index=True)

        # add direction columns
        this_df.loc[:, 'de_direction'] = ['U' if t > 0 else 'D' for t in this_df.loc[:, 'de_logfc']]
        this_df.loc[:, 'dmr_direction'] = ['U' if t > 0 else 'D' for t in this_df.loc[:, 'dmr_median_delta']]
        res[sid] = this_df

    return res