import os
import gffutils
import pandas as pd
from settings import LOCAL_DATA_DIR
from utils.log import get_console_logger
logger = get_console_logger(__name__)


gtf_dir_human = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.release87', 'gtf')
gtf_fn_human = os.path.join(gtf_dir_human, 'Homo_sapiens.GRCh38.87.gtf')

gtf_dir_mouse = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'mouse', 'ensembl', 'GRCm38.p5', 'gtf')

# GTF input by taxonomy ID
GTF_INPUTS_BY_TAX_ID = {
    9606: os.path.join(gtf_dir_human, 'Homo_sapiens.GRCh38.87.gtf'),
    10090: os.path.join(gtf_dir_mouse, 'Mus_musculus.GRCm38.88.gtf'),
}


def gene_lengths(gtf_fn=gtf_fn_human, overwrite=False):
    gtf_dir, filename = os.path.split(gtf_fn)
    filestem, _ = os.path.splitext(filename)
    out_fn = os.path.join(gtf_dir, "%s.gene_lengths.csv" % filestem)

    if os.path.isfile(out_fn) and not overwrite:
        logger.info("Loading gene lengths from existing file %s.", out_fn)
        return pd.read_csv(out_fn, index_col=0, header=None)

    db_fn = os.path.join(gtf_dir, 'gtf.db')
    if os.path.exists(db_fn):
        logger.info("Using existing DB at %s.", db_fn)
        db = gffutils.FeatureDB(db_fn, keep_order=True)
    else:
        logger.info("Creating DB at %s using input %s. This could take several minutes.", db_fn, gtf_fn)
        db = gffutils.create_db(gtf_fn, db_fn, disable_infer_genes=True, disable_infer_transcripts=True)
        logger.info("Successfully created DB.")
    genes = list(db.features_of_type('gene'))
    gene_lengths = dict([(g.id, db.children_bp(g, child_featuretype='exon', merge=True)) for g in genes])
    gene_lengths = pd.Series(gene_lengths)

    gene_lengths.to_csv(out_fn)
    logger.info("Wrote gene lengths to %s.", out_fn)
    ## FIXME: the first time we run this, "some bug" results in an Exception being raised
    ## the "some bug" goes away on a re-run.
    return gene_lengths


def gene_length_by_tax_id(tax_id, **kwargs):
    gtf_fn = GTF_INPUTS_BY_TAX_ID[tax_id]
    return gene_lengths(gtf_fn=gtf_fn, **kwargs)