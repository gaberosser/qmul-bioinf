import os
import pandas as pd
from settings import GIT_LFS_DATA_DIR, HGIC_LOCAL_DIR


if __name__ == '__main__':
    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']
    comparisons = ['syngeneic', 'gibco', 'h9']
    ipa_indir = os.path.join(
        HGIC_LOCAL_DIR,
        'current/core_pipeline/rnaseq/merged_s1_s2/ipa/pathways/'
    )
    xcell_sign_fn = os.path.join(GIT_LFS_DATA_DIR, 'xcell', 'ESM3_signatures.xlsx')

    xcell_s = pd.read_excel(xcell_sign_fn, header=0, index_row=0)
    xcell_signatures = {}
    for i, row in xcell_s.iterrows():
        xcell_signatures[row.Celltype_Source_ID] = row.iloc[2:].dropna().values

    