import os
from settings import HGIC_LOCAL_DIR
import pandas as pd
from utils import output, setops, excel
from plotting import common
import consts
from matplotlib import pyplot as plt, patches, collections, gridspec
import seaborn as sns
import numpy as np


if __name__ == '__main__':
    top_n = 30
    indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq/s2_syngeneic_vs_reference/ipa/pathways')
    outdir = output.unique_output_dir()

    # keys are the term used in the filename, values are those used in the columns
    comps = {
        'syngeneic': 'syngeneic',
        'h9': 'H9',
        'gibco': 'GIBCO'
    }
    pids = consts.PIDS
    ipa_res = {}
    for c, kc in comps.items():
        p = pd.read_csv(
            os.path.join(indir, 'full_de_%s_logp.txt' % c),
            sep='\t',
            skiprows=2,
            index_col=0,
            header=0,
            usecols=range(11)
        )
        z = pd.read_csv(
            os.path.join(indir, 'full_de_%s_z.txt' % c),
            sep='\t',
            skiprows=2,
            index_col=0,
            header=0,
            usecols=range(11)
        ).loc[p.index]
        p.columns = p.columns.str.replace('_logFC', '')
        z.columns = z.columns.str.replace('_logFC', '')
        for pid in pids:
            t = pd.DataFrame(index=p.index, columns=['-logp', 'z'])
            t.loc[:, '-logp'] = p["%s_%s" % (pid, kc)]
            t.loc[:, 'z'] = z["%s_%s" % (pid, kc)]
            # remove results with no significance at all
            t = t.loc[t['-logp'] > 1e-3]
            # unicode index
            t.index = [x.decode('utf-8') for x in t.index]
            ipa_res['%s_%s' % (pid, c)] = t

    excel.pandas_to_excel(ipa_res, os.path.join(outdir, "ipa_results_s2_de.xlsx"))

    # generate a 'top 30' heatmap for each of the PIDs

    alpha = 0.01

    # first pass: to obtain ordering of the pathways:
    # 1) 'Most shared' between reference and syngeneic
    # 2) z, where present

    all_p = {}
    all_z = {}
    all_in = {}
    tot_in = {}

    for i, pid in enumerate(pids):
        this_data = []
        this_p = []
        this_z = []
        tot_in[pid] = None
        for c in sorted(comps.keys()):
            # the_dat = ipa_res['%s_%s' % (pid, c)].sort_values(by='-logp', ascending=False)
            # # just keep z values
            # the_dat = the_dat.loc[the_dat.index[:top_n]]

            the_dat = ipa_res['%s_%s' % (pid, c)]
            the_dat = the_dat.loc[the_dat['-logp'] > -np.log10(alpha)]

            all_in["%s_%s" % (pid, c)] = (~the_dat['-logp'].isnull())
            all_p["%s_%s" % (pid, c)] = the_dat['-logp']
            all_z["%s_%s" % (pid, c)] = the_dat['z']

            if tot_in[pid] is None:
                tot_in[pid] = all_in["%s_%s" % (pid, c)].astype(int).fillna(0)
            else:
                tot_in[pid] = tot_in[pid].add(all_in["%s_%s" % (pid, c)].astype(int).fillna(0), fill_value=0.)

    tot_in = pd.DataFrame(tot_in).fillna(0)
    all_z = pd.DataFrame(all_z)
    all_p = pd.DataFrame(all_p)

    for_sorting = pd.concat((tot_in.sum(axis=1), all_z.sum(axis=1)), axis=1)
    for_sorting.columns = ['n_members', 'sum_z']
    p_order = for_sorting.sort_values(by=['n_members', 'sum_z'], ascending=False).index

    z_min = all_z.min().min()
    z_max = all_z.max().max()

    fig = plt.figure(figsize=(15, 5))
    gs_kw = dict(
        left=0.06,
        right=0.96,
        top=0.99,
        bottom=0.1,
        wspace=0.03,
        hspace=0.1,
    )

    # set up axis grid
    gs = gridspec.GridSpec(
        nrows=len(pids),
        ncols=2,
        width_ratios=[20, 1],
        **gs_kw
    )
    cax = fig.add_subplot(gs[1:-1, 1])
    # fig, axs = plt.subplots(nrows=len(pids), ncols=1, sharex=True, sharey=True, figsize=(15, 5))

    for i, pid in enumerate(pids):
        ax = fig.add_subplot(gs[i, 0])
        # ax = axs[i]

        this_p = []
        this_z = []
        for j, c in enumerate(sorted(comps.keys())):
            # embed this set of data in full sorted list
            the_z = all_z.loc[:, '%s_%s' % (pid, c)].reindex(p_order)
            the_p = all_p.loc[:, '%s_%s' % (pid, c)].reindex(p_order)
            this_z.append(the_z)
            this_p.append(the_p)
        this_z = pd.concat(this_z, axis=1).transpose().iloc[::-1]
        this_p = pd.concat(this_p, axis=1).transpose().iloc[::-1]
        # Z not null, P not null: pathway is enriched and has direction
        sns.heatmap(
            this_z,
            mask=this_p.isnull() | this_z.isnull(),
            cmap='RdBu_r',
            # square=True,
            linewidths=.2,
            linecolor='w',
            vmin=z_min,
            vmax=z_max,
            xticklabels=False,
            cbar=i == 0,
            cbar_ax=cax,
            ax=ax
        )
        # Z null, P not null: no direction information, but pathway is enriched
        tmp = this_z.fillna(1.)
        tmp[tmp != 1] = 0.
        sns.heatmap(
            tmp,
            mask=(~this_z.isnull()) | this_p.isnull(),
            cmap='binary',
            # square=True,
            linewidths=.2,
            linecolor='w',
            xticklabels=False,
            cbar=False,
            ax=ax
        )
        plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
        ax.yaxis.set_ticklabels(['Gibco', 'H9', 'Syngen.'], rotation=0.)
        ax.set_ylabel(pid)
        ax.set_facecolor('0.6')
        cax.set_ylabel('Z score')

        fig.savefig(os.path.join(outdir, "all_sign_pathways_by_zvalue.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "all_sign_pathways_by_zvalue.tiff"), dpi=200)

