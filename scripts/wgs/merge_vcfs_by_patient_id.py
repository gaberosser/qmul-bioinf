import os
import gzip
from glob import glob
import pickle
import re
import itertools
import vcf.utils
import collections
import pandas as pd
from utils import setops, output, log, dictionary
from settings import DATA_DIR
from scripts.hgic_final import consts
import subprocess


logger = log.get_console_logger()


if __name__ == '__main__':
    """
    Here we are going to iterate over a VCF file containing all variants called in our WGS dataset.
    We are looking specifically for any alterations in or in the vicinity of a pre-defined list of candidate genes
    known to be involved in DNA methylation.

    While we are at it, we may generate some overview plots, too.
    """
    pids = consts.PIDS
    contigs = set(['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])
    contigs_str = ','.join(contigs)

    ####### V2: iterate over pre-made short files and store data in memory
    base_indir = os.path.join(DATA_DIR, 'wgs', 'x17067/2017-12-12/')
    meta_fn = os.path.join(DATA_DIR, 'wgs', 'x17067/2017-12-12/', 'sources.csv')

    meta = pd.read_csv(meta_fn, header=0, index_col=0)
    meta.loc[:, 'patient_id'] = ["%03d" % t for t in meta.patient_id]

    outdir = output.unique_output_dir()

    for pid in pids:
        logger.info("Patient %s", pid)
        this_meta = meta.loc[meta.patient_id == pid]
        in_fns = []
        types = []
        for t in this_meta.index:
            the_fn = os.path.join(base_indir, t, "%s.vcf.gz" % meta.loc[t, 'sample'])
            in_fns.append(the_fn)
            types.append(this_meta.loc[t, 'type'])
        logger.info("Found %d conditions to combine: %s", len(types), ', '.join(types))

        out_fn = os.path.join(outdir, "%s.vcf.gz" % pid)

        cmd = "bcftools merge -r {regions} -O z -o {out_fn} {files}".format(
            regions=contigs_str,
            out_fn=out_fn,
            files=' '.join(in_fns)
        )

        subprocess.call(cmd, shell=True)
