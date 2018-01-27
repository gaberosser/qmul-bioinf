import pandas as pd
import os
import warnings
import glob
import random
from settings import LOCAL_DATA_DIR


def reference_genome_chrom_lengths(tax_id=9606, version=None, fn=None):
    """
    Get the lengths of each chromosome in the reference genome.
    :param tax_id: Taxonomy ID
    :param version: Not currently supported, but will allow the specification of more than one version.
    :param fn: If supplied, this directly overrides any other inputs and is used.
    """
    if fn is None:
        if tax_id == 9606:
            if version is None:
                indir = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.p10.release90', 'fa')
            else:
                raise NotImplementedError("Not yet supporting multiple versions.")
        elif tax_id == 10090:
            if version is None:
                indir = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'mouse', 'ensembl', 'GRCm38.p5.r88', 'fa')
            else:
                raise NotImplementedError("Not yet supporting multiple versions.")
        else:
            raise ValueError("Unsupported tax_id: %d" % tax_id)

        flist = glob.glob(os.path.join(indir, '*.fai'))
        if len(flist) == 0:
            raise AttributeError("No .fai files found in directory %s." % indir)
        elif len(flist) > 1:
            warnings.warn("Found >1 .fai files. Using an arbitrary choice.")
        fn = flist[0]

    fai = pd.read_csv(fn, sep='\t', header=None, index_col=0)
    return fai.iloc[:, 0]


def random_genomic_interval(chrom_lengths, n_bp):
    """
    Generate a random genomic interval based on the input chromosome lengths
    """
    the_chr = chrom_lengths.index[random.randint(0, chrom_lengths.size - 1)]
    l = chrom_lengths.loc[the_chr]
    the_start = random.randint(0, l - n_bp)
    the_end = the_start + n_bp - 1
    return (the_chr, the_start, the_end)