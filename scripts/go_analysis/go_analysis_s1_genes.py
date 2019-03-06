from goatools import base, obo_parser, associations, go_enrichment
import wget
import os
import collections
import pandas as pd

from settings import LOCAL_DATA_DIR, HGIC_LOCAL_DIR
from utils import log, output
from scripts.hgic_final import consts


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
    goea_res = {}
    for pid in pids:
        this_ids = de_res_entrez.index[de_res_entrez[pid] == 'Y']
        goea_res[pid] = goe_obj.run_study(this_ids)
        goe_obj.wr_xlsx(os.path.join(outdir, "goea_de_%s.xlsx" % pid), goea_res[pid])