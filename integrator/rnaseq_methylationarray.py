import pandas as pd
import numpy as np
from methylation import dmr


def compute_joint_de_dmr(
        dmr_results,
        de_results,
        de_gene_column='Gene Symbol',
):
    res = {}

    for sid in dmr_results:
        res[sid] = {}

        de_cols = de_results.values()[0].columns.tolist()
        de_col_dtypes = de_results.values()[0].dtypes

        meth_cols = ['chr', 'me_cid', 'me_mediandelta', 'me_median1', 'me_median2', 'me_fdr']
        meth_col_dtypes = pd.Series(['object', 'int', 'float', 'float', 'float', 'float'], index=meth_cols)

        master_df = pd.DataFrame(columns=de_cols + meth_cols)
        master_df = master_df.astype(pd.concat((de_col_dtypes, meth_col_dtypes)).to_dict())

        meth_attrs = ['median_change', 'median1', 'median2']

        for (chr, cls, cid), attrs in dmr.dict_iterator(dmr_results[sid], n_level=3):
            res[sid].setdefault(cls, master_df.copy())

            if len(attrs['genes']) == 0:
                continue

            try:
                # matching entry in DE (by gene name)
                de_match = de_results[sid].loc[de_results[sid].loc[:, de_gene_column].isin(attrs['genes'])]

                if de_match.shape[0] > 0:
                    # form the DMR data block by repeating the same row
                    me_data = np.tile(
                        [chr, cid] + [attrs[k] for k in meth_attrs] + [attrs['padj']],
                        (de_match.shape[0], 1)
                    )
                    me_match = pd.DataFrame(data=me_data, columns=meth_cols, index=de_match.index)\
                        .astype(meth_col_dtypes.to_dict())

                    this_match = pd.concat((de_match, me_match), axis=1)
                    res[sid][cls] = pd.concat(
                        (res[sid][cls], this_match), axis=0, ignore_index=True
                    )
            except Exception as exc:
                # cast to int reqd here in case we loaded from JSON (int keys become str)
                print "Failed to add data: (%s, %s, %d)" % (chr, cls, int(cid))
                print repr(exc)
                continue

        # combine all methylation cluster classes
        res[sid]['all'] = pd.concat(res[sid].values(), axis=0, ignore_index=True)

    return res