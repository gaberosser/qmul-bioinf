import os
import gffutils
from settings import DATA_DIR, DATA_DIR_NON_GIT


if __name__ == '__main__':
    fn = '/home/gabriel/local_data/reference_genomes/ensembl/GRCh38/gtf/Homo_sapiens.GRCh38.87.gtf'
    db = gffutils.create_db(fn, ':memory:', disable_infer_genes=True, disable_infer_transcripts=True)
    exons_by_gene = {}
    # iterate over exons
    for ex in db.features_of_type('exon'):
        t = {
            'chrom': ex.chrom,
            'start': ex.start,
            'stop': ex.stop,
            'strand': ex.strand,
        }
        for gid in ex.attributes['gene_id']:
            exons_by_gene.setdefault(gid, []).append(t)