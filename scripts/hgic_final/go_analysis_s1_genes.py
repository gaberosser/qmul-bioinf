import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from cytoscape import cyto
from scripts.hgic_final import consts
from settings import HGIC_LOCAL_DIR
from utils import log, output, excel, go_analysis, network, reference_genomes

logger = log.get_console_logger()


def reannotate(this_res):
    """
    Convert Ensembl IDs back to gene symbols for easier intepretability
    :param this_res:
    :return:
    """
    # get all Entrez gene IDs and convert in one go
    all_genes = set()
    for k, df in this_res.items():
        for t in df.study_items.str.split(', ').dropna():
            all_genes.update(t)
    gene_conv = reference_genomes.ensembl_to_gene_symbol(sorted(all_genes))

    new_res = {}

    for k in this_res.keys():
        df = this_res[k].copy()
        this_gene_symb = []
        for t in df.study_items:
            if pd.isnull(t):
                this_gene_symb.append('')
            else:
                this_gene_symb.append(','.join(gene_conv.loc[t.split(', ')].dropna().values))
        df.drop('study_items', axis=1, inplace=True)
        df.insert(df.shape[1], 'genes_in_term', this_gene_symb)
        new_res[k] = df
    return new_res


if __name__ == "__main__":
    """
    Aim:
    Load DE / DM genes previously generated in the context of the hGIC project and run GO analysis on them.
    """
    alpha = 0.01
    min_n_genes = 2
    namespaces = {'BP', 'MF'}

    patient_colours = {
        '018': '#ccffcc',
        '019': '#4dff4d',
        '030': '#00cc00',
        '031': '#004d00',
        '017': '#ffcccc',
        '050': '#ff4d4d',
        '054': '#cc0000',
        '061': '#660000',
        '026': '#ff80ff',
        '052': '#800080'
    }

    pids = consts.PIDS
    outdir = output.unique_output_dir()
    # load previously generated DE results
    fn = os.path.join(
        HGIC_LOCAL_DIR,
        'current',
        'core_pipeline',
        'rnaseq',
        'full_de_syngeneic_only.xlsx'
    )
    de_res = pd.read_excel(fn, header=0, index_col=0)

    # run GO
    bg_genes = de_res.index
    go_obj = go_analysis.GOAnalysis()
    gc = go_analysis.ens_to_entrez(bg_genes)
    go_obj.set_gene_conversion(gc)
    go_obj.set_bg_genes(bg_genes)

    # for each GIC / comparator, run GO on the 'syngeneic only' genes
    # this results in VERY few terms
    goea_res = {}
    for pid in pids:
        this_ids = de_res.index[de_res[pid] == 'Y']
        goea_res[pid] = go_obj.run_one(this_ids)

    # excel.pandas_to_excel(goea_res, os.path.join(outdir, "goea_res_all.xlsx"))

    # filter
    goea_res_filt = {}
    for k, df in goea_res.items():
        this_ix = (df.p_bonferroni <= alpha) & (df.NS.isin(namespaces)) & (df.enrichment == 'e')
        df_filt = df.loc[this_ix]
        # include bottom-most nodes only
        ix = []
        for go_id in df_filt.index:
            ix.append(len(go_obj.obo[go_id].get_all_children()) == 0)
        goea_res_filt[k] = df_filt.loc[ix]

    # minor manipulation of results, then save to a single Excel file
    # do this for full results and filtered

    all_res = reannotate(goea_res)
    all_res_filt = reannotate(goea_res_filt)

    excel.pandas_to_excel(all_res, os.path.join(outdir, "goea_de_all_results.xlsx"))
    excel.pandas_to_excel(all_res_filt, os.path.join(outdir, "goea_de_all_results_filtered.xlsx"))

    # create (mega)heatmap of all results
    tmp = pd.concat([v.name for v in all_res_filt.values()])
    tmp = tmp.loc[~tmp.duplicated()]

    for_plot = pd.DataFrame(index=tmp.values)

    for pid in pids[::-1]:
        this = all_res[pid].reindex(tmp.index)
        this.index = tmp.values
        for_plot.insert(0, pid, -np.log10(this['p_bonferroni']))

    # reorder based on mean logp
    for_plot = for_plot.loc[for_plot.mean(axis=1).sort_values(ascending=False).index]

    # abbreviate names
    abbr_rep = [
        ('development', 'dev.'),
        ('genesis', 'gen.'),
        ('response', 'resp.'),
        ('differentiation', 'diff.'),
        ('extracellular matrix', 'ECM'),
        ('regulation', 'reg.'),
        ('signaling pathway', 'signal. path.'),
        ('membrane', 'membr.'),
        ('RNA polymerase II', 'RNA Pol II'),
        ('cell adhesive protein', 'CAP'),
        ('calcium', 'Ca'),
        ('communication', 'comm.'),
    ]
    ix = for_plot.index
    for u, v in abbr_rep:
        ix = ix.str.replace(u, v)

    for_plot.index = ix

    fig = plt.figure(figsize=(6.5, 8.))
    ax = fig.add_subplot(111)
    h = sns.heatmap(
        for_plot,
        mask=for_plot.abs() < 1e-3,
        cmap='YlOrRd',
        linewidths=.2,
        linecolor='w',
        vmin=1.,
        vmax=8,
        yticklabels=True,
        cbar=True,
        cbar_kws={"shrink": 0.7},
        ax=ax,
    )
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0, fontsize=9)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90, fontsize=9)
    fig.subplots_adjust(left=0.7, bottom=0.1, top=0.98)

    # the cbar axis is the one that isn't the main axis!
    cax = [a for a in fig.get_axes() if a is not ax][0]
    cax.set_title(r'$-\log_{10}(p)$', fontsize=9, horizontalalignment='left')

    fig.savefig(os.path.join(outdir, "s1_de_go_terms_filtered_heatmap.png"), dpi=200)

    # Cytoscape
    min_edge_count = 10.
    attr_cols = ['p_bonferroni']
    p_to_g = {}
    node_attrs = {}

    for k, df in all_res_filt.items():
        p_to_g[k] = {}
        node_attrs[k] = {}
        for p, row in df.iterrows():
            nm = row['name']
            this_genes = row.genes_in_term.split(',')
            p_to_g[k][nm] = this_genes
            node_attrs[k][nm] = dict(df.loc[p][attr_cols])

    gg = network.nx_graph_from_multiple_member_dicts(
        p_to_g,
        member_key='genes',
        participants_key='patients',
        name='DE S1 GO',
        min_edge_count=min_edge_count,
        dict_of_node_attrs=node_attrs,
        use_namespace_for_individuals=False
    )

    # Cytoscape session
    cy_session = cyto.CytoscapeSession()
    cyto_nets = {}

    this_net = cy_session.add_networkx_graph(gg, name=gg.name)

    cyto_nets[gg.name] = this_net

    this_net.passthrough_node_label('name')
    this_net.passthrough_node_size_linear('n_genes')
    this_net.passthrough_edge_width_linear('n_genes', xmin=min_edge_count, ymin=0.4, ymax=5)
    this_net.set_node_border_width(0.)
    this_net.set_edge_colour('#b7b7b7')
    this_net.set_node_fill_colour('#ffffff')
    this_net.set_node_transparency(255)

    this_net.node_pie_charts(pids, colours=[patient_colours[p] for p in pids])

    layout_kwargs = dict(
        EdgeAttribute='n_genes',
        defaultSpringLength=20,
        defaultSpringCoefficient=20,
        maxWeightCutoff=max([v['n_genes'] for v in gg.nodes.values()]),
    )

    cy_session.apply_layout(**layout_kwargs)

    cy_session.cy_cmd.session.save(os.path.join(outdir, "go_cytoscape.cys"))
