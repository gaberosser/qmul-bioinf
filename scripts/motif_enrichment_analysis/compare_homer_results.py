import pandas as pd
from matplotlib import pyplot as plt
from plotting import venn

import seaborn as sns
import os
from settings import OUTPUT_DIR
import numpy as np
from utils import excel, powerpoint


def process_results(df, fdr=0.01):
    """
    Modify raw Homer results dataframe in-place.
    :param df:
    :return:
    """
    tmp = df.loc[:, 'Motif Name'].str.split('/')
    nm = tmp.apply(lambda x: x[0])
    src = tmp.apply(lambda x: x[1])
    df.insert(0, 'source', src)
    df.insert(0, 'name', nm)

    k_bg = df.columns[df.columns.str.contains("% of Background Sequences")][0]
    df.insert(2, 'target_pct', df["% of Target Sequences with Motif"].str.replace('%', '').astype(float))
    df.insert(3, 'bg_pct', df[k_bg].str.replace('%', '').astype(float) + 0.01)
    df.insert(2, 'log2_enrichment_factor', np.log2(df.target_pct / df.bg_pct))

    df = df.loc[df['q-value (Benjamini)'] < fdr]
    return df


if __name__ == '__main__':
    fdr = 0.01
    base_dir = os.path.join(OUTPUT_DIR, 'dmr_without_classes.3')
    outdir = base_dir
    pids = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']
    required_comparisons = {
        '017': 'hypo',
        '018': 'hyper',
        '019': 'hypo',
        '030': 'hypo',
        '031': 'hypo',
        '050': 'hyper',
        '054': 'hyper',
        '061': 'hyper',
        '026': 'hyper',
        '052': 'hyper',
    }

    res = {}
    idx_for_upset = {}

    for p in pids:
        res[p] = {}
        # for typ in required_comparisons[p]:
        for typ in ['hypo', 'hyper', 'dmr']:
            # load patient-specific result; use this for results tbl
            # load full result for comparison; add this in showing enrichment
            # create plot (?) with motif (?)
            subdir = os.path.join(base_dir, "%s_%s_oligo_mappings" % (p, typ))
            subdir_full = os.path.join(base_dir, "%s_%s_oligo_mappings_full" % (p, typ))
            if not os.path.isdir(subdir):
                raise Exception("No directory %s" % subdir)
            if not os.path.isdir(subdir_full):
                raise Exception("No directory %s" % subdir_full)

            fn = os.path.join(subdir, "knownResults.txt")
            df = process_results(pd.read_csv(fn, sep='\t', header=0, index_col=None), fdr=fdr)

            fn_full = os.path.join(subdir_full, "knownResults.txt")
            df_full = process_results(pd.read_csv(fn_full, sep='\t', header=0, index_col=None), fdr=fdr)

            this_dat = []
            for _, row in df.iterrows():
                this_dat.append({
                    'name': row['name'],
                    'source': row['source'],
                    'log2(enrichment)': row.log2_enrichment_factor,
                    'fdr': row['q-value (Benjamini)'],
                    'in_full_list': row['Motif Name'] in df_full['Motif Name'].values
                })

            if len(this_dat) > 0:
                this_res = pd.DataFrame.from_dict(this_dat)[
                    ['name', 'source', 'log2(enrichment)', 'fdr', 'in_full_list']
                ]
            else:
                this_res = pd.DataFrame(columns=['name', 'source', 'log2(enrichment)', 'fdr', 'in_full_list'])

            res[p][typ] = this_res

            if typ == required_comparisons[p]:
                idx_for_upset[p] = df.loc[:, 'Motif Name']

    to_xls = {}
    for p, typ in required_comparisons.items():
        to_xls['%s_%s' % (p, typ)] = res[p][typ]
        powerpoint.df_to_powerpoint(os.path.join(outdir, '%s_%s_table.pptx' % (p, typ)), res[p][typ])

    excel.pandas_to_excel(to_xls, os.path.join(outdir, 'patient_specific_direction_specific_dmr.xlsx'))

    # upset
    idx_for_upset = [idx_for_upset[p] for p in pids]
    venn.upset_set_size_plot(idx_for_upset, pids, n_plot=15)
    set_colours = []

