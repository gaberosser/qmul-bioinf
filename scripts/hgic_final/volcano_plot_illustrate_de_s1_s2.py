"""
Created on 2018-09-24
Aim: recreate one of LG's figures: a correlation plot showing all DE genes for one patient, with patient-specific
and syngeneic-only genes highlighted, plus highlighting genes involved in a specific pathway.

I'm going to play around with this and will probably show a volcano plot instead of a correlation plot.

Rather than reload all samples and running lumped DE all over again, I'll load the pre-processed results.
"""
import os
import pickle

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, markers

from scripts.hgic_final import consts
from settings import HGIC_LOCAL_DIR, INTERMEDIATE_DIR
from utils import output, setops, reference_genomes


def de_results_hash(pids, de_params):
    hash_elements = tuple(sorted(pids)) + tuple([de_params[k] for k in sorted(de_params.keys())])
    return hash(hash_elements)


def scatter_plot(
    pid,
    de_res_s1,
    de_res_s1_full,
    de_res_s2,
    ipa_df,
    pw_colours,
    ms=30,
    s1_colour='darkblue',
    s2_colour='lightskyblue',
    external_ref_labels=('GIBCO', 'H9')
):
    external_ref_labels = list(external_ref_labels)

    dat_all = pd.DataFrame.from_dict({
        'logFC': de_res_s1_full[pid].logFC,
        'pFDR': -np.log10(de_res_s1_full[pid].FDR)
    })

    # patient specific

    vs, vc = setops.venn_from_arrays(*[de_res_s1[p].index for p in pids])
    set_ps = list(setops.binary_combinations_sum_eq(len(pids), 1))

    # remember that sets are defined in reverse!
    set_ps = dict([(pids[i], vs[set_ps[(-1 - i)]]) for i in range(len(pids))])

    dat_patient_specific = dat_all.loc[set_ps[pid]]

    # syngeneic only
    the_idx = [de_res_s2[(pid, k)].index for k in [pid] + external_ref_labels]
    vs, vc = setops.venn_from_arrays(*the_idx)

    dat_syngen_only = dat_all.loc[vs['100']]

    # pathway(s)

    in_pathways = {}
    pw_dat = {}

    for pw in pw_colours:
        the_genes = ipa_df.loc[pw].genes.split(',')
        the_ens_id = reference_genomes.gene_symbol_to_ensembl(the_genes).dropna().tolist()
        # resolve '/' weirdness
        for g in the_genes:
            if '/' in g:
                g_arr = g.split('/')
                e_arr = reference_genomes.gene_symbol_to_ensembl(g_arr).dropna().tolist()
                the_ens_id += e_arr
        # drop any genes not found in the data
        in_pathways[pw] = dat_all.index.intersection(the_ens_id)
        pw_dat[pw] = dat_all.loc[in_pathways[pw]]

    marker_full = markers.MarkerStyle(marker='o', fillstyle='full')
    marker_s1 = markers.MarkerStyle(marker='o', fillstyle='left')
    marker_s2 = markers.MarkerStyle(marker='o', fillstyle='right')

    fig = plt.figure(figsize=(7.1, 6.3))
    ax = fig.add_subplot(111)

    ax.scatter(
        dat_all['logFC'],
        dat_all['pFDR'],
        c='gray',
        s=ms,
        zorder=1,
        alpha=0.3,
        label='All DE'
    )
    # S1
    ax.scatter(
        dat_patient_specific.loc[:, 'logFC'],
        dat_patient_specific.loc[:, 'pFDR'],
        marker=marker_s1,
        c=s1_colour,
        s=ms,
        zorder=3,  # in front as there are fewer
        label='Patient specific'
    )
    ax.scatter(
        dat_patient_specific.loc[:, 'logFC'],
        dat_patient_specific.loc[:, 'pFDR'],
        marker=marker_full,
        edgecolor='k',
        linewidths=.5,
        c='none',
        s=ms,
        zorder=3,  # in front as there are fewer
        label=None,
    )
    # S2
    ax.scatter(
        dat_syngen_only['logFC'],
        dat_syngen_only['pFDR'],
        marker=marker_s2,
        c=s2_colour,
        s=ms,
        zorder=2,  # behind S1 as there are more
        label='Syngeneic only'
    )
    ax.scatter(
        dat_syngen_only['logFC'],
        dat_syngen_only['pFDR'],
        marker=marker_full,
        edgecolor='k',
        linewidths=.5,
        c='none',
        s=ms,
        zorder=2,  # behind S1 as there are more
        label=None,
    )

    for pw in pw_colours:
        ax.scatter(
            pw_dat[pw]['logFC'],
            pw_dat[pw]['pFDR'],
            c='none',
            edgecolor=pw_colours[pw],
            linewidths=1.5,
            s=ms + 2,
            zorder=4,
            label=pw
        )

    ax.legend(loc='upper right', frameon=True, facecolor='w', framealpha=0.7)
    ax.set_xlabel('logFC')
    ax.set_ylabel('-log10(FDR)')
    fig.tight_layout()

    return fig, ax


pids = consts.PIDS


DE_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'de')
outdir = output.unique_output_dir()

external_ref_labels = ['GIBCO', 'H9']

# params and sample names needed to ensure we load the right DE results
de_params = {
    'lfc': 1,
    'fdr': 0.01,
    'method': 'QLGLM'
}

sample_names_s1 = [
    'GBM018_P12',
    'GBM018_P10',
    'DURA018_NSC_N4_P4',
    'DURA018_NSC_N2_P6',
    'GBM019_P4',
    'GBM019_P3n6',
    'DURA019_NSC_N8C_P2',
    'DURA019_NSC_N5C1_P2',
    'GBM030_P9n10',
    'GBM030_P5',
    'DURA030_NSC_N16B6_P1',
    'DURA030_NSC_N9_P2',
    'GBM031_P7',
    'GBM031_P4',
    'DURA031_NSC_N44B_P2',
    'DURA031_NSC_N44_P3',
    'GBM017_P3',
    'GBM017_P4',
    'DURA017_NSC_N3C5_P4',
    'GBM050_P7n8',
    'GBM050_P9',
    'DURA050_NSC_N12_P3',
    'DURA050_NSC_N16_P4',
    'GBM054_P4',
    'GBM054_P6',
    'DURA054_NSC_N3C_P2',
    'DURA054_NSC_N2E_P1',
    'GBM061_P3',
    'GBM061_P5',
    'DURA061_NSC_N4_P2',
    'DURA061_NSC_N1_P3',
    'GBM026_P8',
    'GBM026_P3n4',
    'DURA026_NSC_N31D_P5',
    'GBM052_P6n7',
    'GBM052_P4n5',
    'DURA052_NSC_N4_P3',
    'DURA052_NSC_N5_P2'
]

sample_names_s2 =[
    'GBM018_P12',
    'GBM018_P10',
    'DURA018_NSC_N4_P4',
    'DURA018_NSC_N2_P6',
    'GBM019_P4',
    'GBM019_P3n6',
    'DURA019_NSC_N8C_P2',
    'DURA019_NSC_N5C1_P2',
    'GBM030_P9n10',
    'GBM030_P5',
    'DURA030_NSC_N16B6_P1',
    'DURA030_NSC_N9_P2',
    'GBM031_P7',
    'GBM031_P4',
    'DURA031_NSC_N44B_P2',
    'DURA031_NSC_N44_P3',
    'GBM017_P3',
    'GBM017_P4',
    'DURA017_NSC_N3C5_P4',
    'GBM050_P7n8',
    'GBM050_P9',
    'DURA050_NSC_N12_P3',
    'DURA050_NSC_N16_P4',
    'GBM054_P4',
    'GBM054_P6',
    'DURA054_NSC_N3C_P2',
    'DURA054_NSC_N2E_P1',
    'GBM061_P3',
    'GBM061_P5',
    'DURA061_NSC_N4_P2',
    'DURA061_NSC_N1_P3',
    'GBM026_P8',
    'GBM026_P3n4',
    'DURA026_NSC_N31D_P5',
    'GBM052_P6n7',
    'GBM052_P4n5',
    'DURA052_NSC_N4_P3',
    'DURA052_NSC_N5_P2',
    'GIBCO_NSC_P4',
    'H9_NSC_1',
    'H9_NSC_2'
]


# load DE results (or error if not found)

the_hash = de_results_hash(sample_names_s1, de_params)
filename = 'de_results_paired_comparison.%d.pkl' % the_hash
fn_s1 = os.path.join(DE_LOAD_DIR, filename)

if not os.path.isfile(fn_s1):
    raise IOError("No DE results file found: %s" % fn_s1)

with open(fn_s1, 'rb') as f:
    de_res_full_s1 = pickle.load(f)
    de_res_s1 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s1.items()])

the_hash = de_results_hash(sample_names_s2, de_params)
filename = 'de_results_cross_comparison.%d.pkl' % the_hash
fn_s2 = os.path.join(DE_LOAD_DIR, filename)

if not os.path.isfile(fn_s2):
    raise IOError("No DE results file found: %s" % fn_s2)

with open(fn_s2, 'rb') as f:
    de_res_full_s2 = pickle.load(f)
    de_res_s2 = dict([(k, v.loc[v.FDR < de_params['fdr']]) for k, v in de_res_full_s2.items()])

# load relevant pathways
ipa_fn = os.path.join(
    HGIC_LOCAL_DIR,
    'current',
    'core_pipeline',
    'rnaseq',
    's0_individual_patients_direct_comparison',
    'ipa',
    'pathways',
    'ipa_de_top_30_pathways.xlsx'
)

ipa_df = pd.read_excel(ipa_fn, header=0, index_col=0)

ms = 30
s1_colour= 'darkblue'
s2_colour = 'lightskyblue'

# Select one patient and one pathway for display

pid = '050'
# match pathway colours from previous plot
pw_colours = {
    'Aryl Hydrocarbon Receptor Signaling': '#E9605A',
    # 'cAMP-mediated signaling': '#25D86D'
}

fig, ax = scatter_plot(
    pid,
    de_res_s1,
    de_res_full_s1,
    de_res_s2,
    ipa_df,
    pw_colours,
    s1_colour=s1_colour,
    s2_colour=s2_colour,
    ms=ms,
    external_ref_labels=external_ref_labels
)

fig.savefig(os.path.join(outdir, "volcano_plot_patient_%s_S1_S2_ahr_pathway.png" % pid), dpi=200)
fig.savefig(os.path.join(outdir, "volcano_plot_patient_%s_S1_S2_ahr_pathway.tiff" % pid), dpi=200)

pw_colours['cAMP-mediated signaling'] = '#25D86D'

fig, ax = scatter_plot(
    pid,
    de_res_s1,
    de_res_full_s1,
    de_res_s2,
    ipa_df,
    pw_colours,
    s1_colour=s1_colour,
    s2_colour=s2_colour,
    ms=ms,
    external_ref_labels=external_ref_labels
)

fig.savefig(os.path.join(outdir, "volcano_plot_patient_%s_S1_S2_ahr_camp_pathway.png" % pid), dpi=200)
fig.savefig(os.path.join(outdir, "volcano_plot_patient_%s_S1_S2_ahr_camp_pathway.tiff" % pid), dpi=200)


pid = '052'
pw_colours = {
    'Tec Kinase Signaling': '#545ED4'
}

fig, ax = scatter_plot(
    pid,
    de_res_s1,
    de_res_full_s1,
    de_res_s2,
    ipa_df,
    pw_colours,
    s1_colour=s1_colour,
    s2_colour=s2_colour,
    ms=ms,
    external_ref_labels=external_ref_labels
)

fig.savefig(os.path.join(outdir, "volcano_plot_patient_%s_S1_S2_tec_pathway.png" % pid), dpi=200)
fig.savefig(os.path.join(outdir, "volcano_plot_patient_%s_S1_S2_tec_pathway.tiff" % pid), dpi=200)