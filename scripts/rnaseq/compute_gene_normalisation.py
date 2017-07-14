import os
import gffutils
from settings import DATA_DIR_NON_GIT


if __name__ == '__main__':
    outdir = os.path.join(DATA_DIR_NON_GIT, 'reference_genomes', 'ensembl', 'GRCh38', 'gtf')
    fn = os.path.join(outdir, 'Homo_sapiens.GRCh38.87.gtf')
    outfn = os.path.join(outdir, 'gtf.db')
    if os.path.exists(outfn):
        db = gffutils.FeatureDB(outfn, keep_order=True)
    else:
        db = gffutils.create_db(fn, outfn, disable_infer_genes=True, disable_infer_transcripts=True)
    genes = list(db.features_of_type('gene'))
    gene_lengths = dict([(g, db.children_bp(g, child_featuretype='exon')) for g in genes])
#    exons_by_gene = {}
    # iterate over exons
#    for ex in db.features_of_type('exon'):
#        t = {
#            'chrom': ex.chrom,
#            'start': ex.start,
#            'stop': ex.stop,
#            'strand': ex.strand,
#        }
#        for gid in ex.attributes['gene_id']:
#            exons_by_gene.setdefault(gid, []).append()
