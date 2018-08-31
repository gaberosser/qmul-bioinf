import requests
import pandas as pd
import os
from utils import output, log
from load_data import geo_repo, sra
import ftplib

from settings import RNASEQ_DIR

"""
AIM: retrieve metadata (combined SRA and GEO) for a number of GEO datasets in one go
"""

GEO_FTP = "ftp.ncbi.nlm.nih.gov"
GEO_BASE = "/geo/series/{stripped_id}/{full_id}/matrix/{full_id}_series_matrix.txt.gz"

logger = log.get_console_logger("sra_batch_get_metadata")

inputs = [
    ('GSE116124', 'SRP151040'),
    ('GSE97265', 'SRP102810'),
    ('GSE89056', 'SRP091957'),
    ('GSE107654', 'SRP126289'),
    ('GSE97904', 'SRP104149'),
    ('GSE97619', 'SRP103788'),
    ('GSE85839', 'SRP082406'),
    ('GSE53094', 'SRP033569'),
    ('GSE67915', 'SRP057205'),
    # ('GSE62772', 'SRP049340'),  # need to dl two geo files
    ('GSE73211', 'SRP063867'),
]

ftp = ftplib.FTP(host=GEO_FTP, user='', passwd='')
ftp.login()

for gs_id, sr_id in inputs:
    outdir = os.path.join(RNASEQ_DIR, gs_id)
    if not os.path.exists(outdir):
        logger.info("Creating output directory %s" % outdir)
        os.makedirs(outdir)
    elif not os.path.isdir(outdir):
        raise Exception("%s exists but is not a directory." % outdir)

    out_fn = os.path.join(outdir, 'sources.csv')
    geo_out_fn = os.path.join(outdir, '%s_series_matrix.txt.gz' % sr_id)

    srr_tbl = sra.get_runinfo_table(sr_id)

    # generate GEO series file link and retrieve
    stripped_id = gs_id[:-3] + 'nnn'
    geo_dir = GEO_BASE.format(stripped_id=stripped_id, full_id=gs_id)

    with open(geo_out_fn, 'wb') as f:
        ftp.retrbinary('RETR %s' % geo_dir, f.write)

    # combine with SRA data
    df = geo_repo.combine_geo_sra_metadata([geo_out_fn], sr_id)
    df.to_csv(out_fn)