import pandas as pd
import numpy as np
from methylation import dmr


def compute_joint_de_dmr(dmr_results, de_results):
    res = {}

    for sid in dmr_results:
        print sid
        res[sid] = {}

        de_cols = ['genes', 'logFC', 'ensembl', 'direction', 'FDR', 'logCPM']
        meth_cols = ['me_genes', 'chr', 'me_cid', 'me_mediandelta', 'me_median1', 'me_median2', 'me_fdr']
        meth_attrs = ['median_change', 'median1', 'median2']

        for (chr, cls, cid), attrs in dmr.dict_iterator(dmr_results[sid], n_level=3):
            res[sid].setdefault(cls, pd.DataFrame(columns=de_cols + meth_cols))

            if len(attrs['genes']) == 0:
                continue

            try:
                # matching entry in DE (by gene name)
                de_match = de_results[sid].loc[de_results[sid].loc[:, 'genes'].isin(attrs['genes'])]

                if de_match.shape[0] > 0:
                    # form the DMR data block by repeating the same row
                    me_data = np.tile(
                        [chr, cid] + [attrs[k] for k in meth_attrs] + [attrs['padj']],
                        (de_match.shape[0], 1)
                    )
                    me_data = np.concatenate(
                        (
                            np.reshape(de_match.genes.values, (de_match.shape[0], 1)),
                            me_data
                        ),
                        axis=1
                    )
                    me_match = pd.DataFrame(data=me_data, columns=meth_cols, index=de_match.index)

                    this_match = pd.concat((de_match, me_match), axis=1)
                    res[sid][cls] = pd.concat(
                        (res[sid][cls], this_match), axis=0, ignore_index=True
                    )
            except Exception as exc:
                print "Failed to add data: (%s, %s, %d)" % (chr, cls, cid)
                print repr(exc)
                continue

        # combine all methylation cluster classes
        res[sid]['all'] = pd.concat(res[sid].values(), axis=0, ignore_index=True)

    return res