from goatools import base, obo_parser, associations, go_enrichment
import wget
import os
import collections
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from settings import LOCAL_DATA_DIR, HGIC_LOCAL_DIR
from utils import log, output, excel, setops
from scripts.hgic_final import consts
import references
from cytoscape import cyto


logger = log.get_console_logger()


def ens_to_entrez(ens, genetoens_fn):
    entrezid_to_ens = pd.read_csv(genetoens_fn, sep='\t', index_col=None, header=0)
    entrezid_to_ens = entrezid_to_ens.loc[entrezid_to_ens['#tax_id'] == 9606, ['GeneID', 'Ensembl_gene_identifier']]
    entrezid_to_ens.columns = ['entrez_id', 'ensembl_id']
    entrezid_to_ens = entrezid_to_ens.loc[~entrezid_to_ens.entrez_id.duplicated()].set_index('entrez_id')

    ens_to_entrezid = entrezid_to_ens.copy()
    ens_to_entrezid.insert(0, 'entrez_id', ens_to_entrezid.index)
    ens_to_entrezid = ens_to_entrezid.loc[~ens_to_entrezid.ensembl_id.duplicated()]
    ens_to_entrezid.set_index('ensembl_id', inplace=True)

    return ens_to_entrezid.reindex(ens).dropna().astype(int)


if __name__ == "__main__":
    """
    Aim:
    Load DE / DM genes previously generated in the context of the hGIC project and run GO analysis on them.
    """
    alpha = 0.01
    min_n_genes = 2

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

    obo_fn = os.path.join(LOCAL_DATA_DIR, 'gene_ontology', 'current', 'go-basic.obo')
    genetogo_fn = os.path.join(LOCAL_DATA_DIR, 'gene_ontology', 'current', 'gene2go')
    genetoens_fn = os.path.join(LOCAL_DATA_DIR, 'gene_ontology', 'current', 'gene2ensembl.gz')
    genetoens_url = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz"

    # load conversion from Entrez ID to Ensembl
    if not os.path.isfile(genetoens_fn):
        logger.info("Downloading RefGene-Ensembl converter from %s, saving to %s.", genetoens_url, genetoens_fn)
        wget.download(genetoens_url, out=genetoens_fn)

    # entrezid_to_ens = pd.read_csv(genetoens_fn, sep='\t', index_col=None, header=0)
    # entrezid_to_ens = entrezid_to_ens.loc[entrezid_to_ens['#tax_id'] == 9606, ['GeneID', 'Ensembl_gene_identifier']]
    # entrezid_to_ens.columns = ['entrez_id', 'ensembl_id']
    # entrezid_to_ens = entrezid_to_ens.loc[~entrezid_to_ens.entrez_id.duplicated()].set_index('entrez_id')
    #
    # ens_to_entrezid = entrezid_to_ens.copy()
    # ens_to_entrezid.insert(0, 'entrez_id', ens_to_entrezid.index)
    # ens_to_entrezid = ens_to_entrezid.loc[~ens_to_entrezid.ensembl_id.duplicated()]
    # ens_to_entrezid.set_index('ensembl_id', inplace=True)

    # load ontologies
    obo_fn = base.download_go_basic_obo(obo_fn)
    obodag = obo_parser.GODag(obo_fn)

    # load associations
    genetogo_fn = base.download_ncbi_associations(genetogo_fn)
    geneid2gos_hu = associations.read_ncbi_gene2go(genetogo_fn, taxids=[9606])

    print("{N:,} annotated human genes".format(N=len(geneid2gos_hu)))

    # convert DE Ensembl IDs to Entrez
    ens_converted = ens_to_entrez(de_res.index, genetoens_fn)
    de_res_entrez = de_res.copy().loc[ens_converted.index]
    de_res_entrez.index = ens_converted.entrez_id

    # background gene set
    bg_genes = de_res_entrez.index

    # init
    goe_obj = go_enrichment.GOEnrichmentStudy(
        bg_genes,
        geneid2gos_hu,
        obodag,
    )

    # run GOEA for each patient
    # save full output for future use
    goea_res = {}
    for pid in pids:
        this_ids = de_res_entrez.index[de_res_entrez[pid] == 'Y']
        goea_res[pid] = goe_obj.run_study(this_ids)
        goe_obj.wr_xlsx(os.path.join(outdir, "goea_de_%s.xlsx" % pid), goea_res[pid])

    # filter to significant results and keep bottom 'leaves' of the GO hierarchy
    # only include the following namespaces
    namespaces = {'BP', 'MF'}

    goea_res_filt = {}
    for pid in pids:
        # also need to filter on number of genes involved - some are zero but have sign p value?!
        this_res = [
            t for t in goea_res[pid]
            if (t.p_bonferroni <= alpha)
            and (t.ratio_in_study[0] >= min_n_genes)
            and (t.NS in namespaces)
            and (t.enrichment == 'e')
        ]
        # include bottom-most nodes only
        goea_res_filt[pid] = [t for t in this_res if (len(t.goterm.get_all_children()) == 0)]
        goe_obj.wr_xlsx(os.path.join(outdir, "goea_de_%s_filtered.xlsx" % pid), goea_res_filt[pid])

    # reload results, minor manipulation, then save to a single Excel file
    # do this for full results and filtered

    def reload_and_annotate(fn):
        this_res = pd.read_excel(fn, header=0, index_col=0)
        # get all Entrez gene IDs and convert in one go
        all_genes = set()
        for t in this_res.study_items.str.split(', ').dropna():
            all_genes.update([int(x) for x in t])
        gene_conv = references.entrez_to_gene_symbol(sorted(all_genes))
        this_gene_symb = []
        for t in this_res.study_items:
            if pd.isnull(t):
                this_gene_symb.append('')
            else:
                this_gene_symb.append(','.join(gene_conv.loc[[int(x) for x in t.split(', ')]].dropna().values))
        this_res.drop('study_items', axis=1, inplace=True)
        this_res.insert(this_res.shape[1], 'genes_in_term', this_gene_symb)
        return this_res

    all_res = {}
    all_res_filt = {}
    for pid in pids:
        all_res[pid] = reload_and_annotate(os.path.join(outdir, "goea_de_%s.xlsx" % pid))
        all_res_filt[pid] = reload_and_annotate(os.path.join(outdir, "goea_de_%s_filtered.xlsx" % pid))

    excel.pandas_to_excel(all_res_filt, os.path.join(outdir, "goea_de_all_results.xlsx"))
    excel.pandas_to_excel(all_res_filt, os.path.join(outdir, "goea_de_all_results_filtered.xlsx"))

    # create (mega)heatmap of all results
    tmp = pd.concat([v.name for v in all_res.values()])
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
        ('calcium', 'Ca')
    ]
    ix = for_plot.index
    for u, v in abbr_rep:
        ix = ix.str.replace(u, v)

    for_plot.index = ix

    fig = plt.figure(figsize=(6.5, 8.))
    ax = fig.add_subplot(111)
    h = sns.heatmap(
        for_plot,
        mask=for_plot.isnull(),
        cmap='YlOrRd',
        linewidths=.2,
        linecolor='w',
        vmin=1.,
        vmax=8,
        yticklabels=True,
        cbar=True,
        # cbar_ax=cax,
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

