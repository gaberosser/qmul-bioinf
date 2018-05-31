import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import colors
from plotting import venn, common
from methylation import loader

import seaborn as sns
import os
import re
from settings import OUTPUT_DIR
import numpy as np
from utils import excel, powerpoint


def process_results(df, fdr=None):
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

    if fdr is not None:
        df = df.loc[df['q-value (Benjamini)'] < fdr]

    return df


def load_one(indir, pid, typ, cpg_status=None, fdr=None):
    """
    :param indir:
    :param pid:
    :param typ:
    :param cpg_status:
    :return:
    """
    if cpg_status is None:
        subdir = os.path.join(indir, "%s_%s_oligo_mappings" % (pid, typ))
    else:
        subdir = os.path.join(indir, "%s_%s_%s_oligo_mappings" % (pid, typ, cpg_status))
    subdir_full = "%s_full" % subdir

    if not os.path.isdir(subdir):
        raise Exception("No directory %s" % subdir)
    if not os.path.isdir(subdir_full):
        raise Exception("No directory %s" % subdir_full)

    fn = os.path.join(subdir, "knownResults.txt")
    df = process_results(pd.read_csv(fn, sep='\t', header=0, index_col=None), fdr=fdr)

    fn_full = os.path.join(subdir_full, "knownResults.txt")
    df_full = process_results(pd.read_csv(fn_full, sep='\t', header=0, index_col=None), fdr=fdr)

    return df, df_full


def get_n_target_sequences(df):
    """
    Extract the number of target sequences from the supplied results DataFrame.
    :param df:
    :return:
    """
    the_col = df.columns[df.columns.str.contains("# of Target Sequences with Motif")][0]
    r = re.search(
        r'\(of (?P<n>[0-9]*)\)',
        the_col
    )
    return int(r.group('n'))


def motif_results_to_summary(df, df_full):
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
    return this_res


def interval_names_from_bin_edges(edges, add_infty=True, lbl_format='n'):
    """
    Given the iterable of bin edges, generate interval names
    :param edges:
    :param add_infty:
    :param lbl_format:
    :return:
    """
    hi_sub1 = False
    if lbl_format == 'n':
        lbl_str = "{lo} <= n < {hi}"
        lbl_str_top = "n >= {lo}"
    elif lbl_format == 'interval':
        lbl_str = "[{lo}, Infty)"
        lbl_str_top = lbl_str
    elif lbl_format == 'hyphen':
        lbl_str = "{lo} - {hi}"
        lbl_str_top = "{lo} -"
        hi_sub1 = True

    lbl = []
    for i in range(len(edges) - 1):
        lo = edges[i]
        hi = edges[i + 1]
        if hi_sub1:
            hi -= 1
        lbl.append(lbl_str.format(lo=lo, hi=hi))

    if add_infty:
        lbl.append(lbl_str_top.format(lo=edges[-1]))

    return lbl


if __name__ == '__main__':
    fdr = 0.01
    base_dir = os.path.join(OUTPUT_DIR, 'dmr_without_classes')
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
    cpg_statuses = ['island', 'shore', 'shelf', 'open_sea']

    # get background probe distribution across CpG island statuses
    k_open_sea = 'open_sea'
    cats = {
        'N_Shore': 'shore',
        'S_Shore': 'shore',
        'Island': 'island',
        'N_Shelf': 'shelf',
        'S_Shelf': 'shelf',
        k_open_sea: 'open_sea',
    }

    anno = loader.load_illumina_methylationepic_annotation()
    this_counts = anno.loc[:, 'Relation_to_UCSC_CpG_Island'].fillna(k_open_sea).value_counts().to_dict()
    island_counts_all = dict([
        (v, this_counts.get(k, 0)) for k, v in cats.items()
    ])
    bg_dist = pd.Series(island_counts_all)
    bg_dist = bg_dist.divide(bg_dist.sum())

    res = {}
    idx_for_upset = {}

    for p in pids:
        res[p] = {}
        # for typ in required_comparisons[p]:
        for typ in ['hypo', 'hyper', 'dmr']:
            # load patient-specific result; use this for results tbl
            # load full result for comparison; add this in showing enrichment
            # create plot (?) with motif (?)

            df, df_full = load_one(base_dir, p, typ, fdr=fdr)
            res[p][typ] = motif_results_to_summary(df, df_full)

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

    # repeat but now also separate by CpG island status
    fdr = 0.01
    res2 = {}
    motif_counts = {}
    motif_counts_nif = {}

    probe_counts = {}
    probe_counts_full = {}
    probe_dist = {}
    probe_dist_full = {}
    probe_dist_rel_bg = {}
    probe_dist_full_rel_bg = {}

    for p in pids:
        res2[p] = {}
        motif_counts[p] = {}
        motif_counts_nif[p] = {}

        probe_counts[p] = {}
        probe_counts_full[p] = {}

        probe_dist[p] = {}
        probe_dist_full[p] = {}

        probe_dist_rel_bg[p] = {}
        probe_dist_full_rel_bg[p] = {}

        for typ in ['hypo', 'hyper', 'dmr']:
            res2[p][typ] = {}
            the_counts = {}
            the_counts_full = {}

            motif_counts[p][typ] = pd.Series(0, index=cpg_statuses)
            motif_counts_nif[p][typ] = pd.Series(0, index=cpg_statuses)

            for st in cpg_statuses:
                try:
                    df, df_full = load_one(base_dir, p, typ, cpg_status=st, fdr=fdr)
                    res2[p][typ][st] = motif_results_to_summary(df, df_full)
                    the_counts[st] = get_n_target_sequences(df)
                    the_counts_full[st] = get_n_target_sequences(df_full)
                except IOError:
                    # no result - just insert an empty result
                    res2[p][typ][st] = pd.DataFrame(columns=['name', 'source', 'log2(enrichment)', 'fdr', 'in_full_list'])
                    the_counts[st] = 0.
                    the_counts_full[st] = 0.
                motif_counts[p][typ][st] = res2[p][typ][st].shape[0]
                motif_counts_nif[p][typ][st] = (~res2[p][typ][st]['in_full_list']).sum()

            the_counts = pd.Series(the_counts)
            the_counts_full = pd.Series(the_counts_full)

            probe_counts[p][typ] = the_counts
            probe_counts_full[p][typ] = the_counts_full

            probe_dist[p][typ] = the_counts.divide(the_counts.sum())
            probe_dist_full[p][typ] = the_counts_full.divide(the_counts_full.sum())

            probe_dist_rel_bg[p][typ] = probe_dist[p][typ] / bg_dist
            probe_dist_full_rel_bg[p][typ] = probe_dist_full[p][typ] / bg_dist

    motif_count_breaks = [0, 2, 5, 10, 20, 30, 100]
    motif_count_breaks_colnames = interval_names_from_bin_edges(motif_count_breaks, add_infty=True)
    motif_count_breaks_nif = [0, 2, 4, 6, 8, 10]
    motif_count_breaks_nif_colnames = interval_names_from_bin_edges(motif_count_breaks_nif, add_infty=True)

    motif_counts_binned = {}
    motif_counts_nif_binned = {}
    for p in pids:
        motif_counts_binned[p] = {}
        motif_counts_nif_binned[p] = {}
        for typ in ['hypo', 'hyper']:
            motif_counts_binned[p][typ] = pd.Series(index=cpg_statuses)
            motif_counts_nif_binned[p][typ] = pd.Series(index=cpg_statuses)
            tmp = np.digitize(motif_counts[p][typ].values, motif_count_breaks)
            tmp_nif = np.digitize(motif_counts_nif[p][typ].values, motif_count_breaks_nif)
            for i, st in enumerate(cpg_statuses):
                motif_counts_binned[p][typ][st] = motif_count_breaks_colnames[tmp[i] - 1]
                motif_counts_nif_binned[p][typ][st] = motif_count_breaks_nif_colnames[tmp_nif[i] - 1]

    cset = common.COLOUR_BREWERS[len(motif_count_breaks_colnames)]
    colours_motif = pd.Series(dict([
        (motif_count_breaks_colnames[i], cset[i]) for i in range(len(cset))
    ]))

    ## TODO: may wish to make the colour map continuous to make it easier to discern different levels?

    fig, axs = plt.subplots(len(pids), 2, figsize=(4, 8), sharex=True, sharey=False)
    ymax = 0
    x = range(len(cpg_statuses))
    for i, p in enumerate(pids):
        for j, typ in enumerate(['hypo', 'hyper']):
            ax = axs[i][j]
            y = probe_dist_rel_bg[p][typ].loc[cpg_statuses]
            c = colours_motif[motif_counts_binned[p][typ].values].values
            # y = probe_counts[p][typ].loc[cpg_statuses]
            ymax = max(ymax, y.max())
            ax.bar(x, y, color=c)
            ax.set_xticks(x)
            ax.set_xticklabels(cpg_statuses, rotation=90)
            if j == 0:
                ax.set_ylabel(p)
            if i == 0:
                ax.set_title(typ.capitalize())

    ymax = int(np.ceil(1.05 * ymax))

    for i, p in enumerate(pids):
        for j, typ in enumerate(['hypo', 'hyper']):
            # text: number of probes
            ax = axs[i][j]
            ax.text(
                0,
                ymax * 0.95,
                "%d probes" % probe_counts[p][typ].sum(),
                horizontalalignment='left', verticalalignment='top'
            )

    for ax in axs.flat:
        ax.set_ylim([0, ymax])
        ax.set_yticks([0, ymax])
        ax.axhline(1., c='k', ls='--', alpha=0.4)

    fig.subplots_adjust(left=0.1, right=0.99, top=0.95, bottom=0.1, wspace=0.2, hspace=0.4)

