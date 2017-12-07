import ftplib
import urllib2
from utils.output import unique_output_dir
from utils.log import get_console_logger
import pandas as pd
import os
import re
from gzip import GzipFile
logger = get_console_logger(__name__)


"""
For a simple method to get all the data associated with a GEO project, use scripts/bash/download_sra_project_fastq.sh
Get associated SRA attributes with sra.get_runinfo_table
"""


def download_from_ftp(dl_paths, outdir=None, host='ftp-trace.ncbi.nlm.nih.gov', user='', passwd=''):
    """
    :param dl_paths: List of paths to download
    :param host:
    :param user:
    :param passwd:
    :return:
    """
    if outdir is None:
        outdir = __name__
    outdir = unique_output_dir(outdir)

    ftp = ftplib.FTP(host=host, user=user, passwd=passwd)
    ftp.login()

    for ff in dl_paths:

        # get download filename
        sp = [t for t in ff.split('/') if len(t)]
        sp = sp[-1]
        outfile = os.path.join(outdir, sp)
        logger.info("Attempting to download %s to %s", ff, outfile)
        with open(outfile, 'wb') as f:
            try:
                ftp.retrbinary('RETR %s' % ff, f.write)
            except Exception:
                logger.exception("Download failed")


def load_metadata_from_series_matrix(infile):
    """
    Load metadata from the specified text matrix (available for many GEO datasets)
    :param infile: .txt or .txt.gz series matrix
    :return:
    """
    meta_map = {
        'Sample_title': 'title',
        'Sample_geo_accession': 'accession',
        'Sample_source_name_ch1': 'description',
    }
    nested_headers = (
        'Sample_characteristics_ch1',
    )

    if re.search(re.compile(r'\.gz$', re.IGNORECASE), infile):
        opener = lambda x: GzipFile(filename=x, mode='rb')
    else:
        opener = lambda x: open(x, 'rb')

    meta = pd.DataFrame()

    with opener(infile) as f:
        while True:
            line = f.readline().strip('\n')
            if len(line) == 0:
                continue
            header = re.match(r'^!(?P<hd>[^\t]*)', line).group('hd')
            if header in nested_headers:
                the_data = [t.strip('"') for t in line.split('\t')[1:]]
                hdr = re.match(r'^(?P<hd>[^:]*)', the_data[0]).group('hd')
                the_data = [re.sub(r'%s: ' % hdr, '', t) for t in the_data]
                meta[hdr] = the_data
            elif header in meta_map:
                the_data = [t.strip('"') for t in line.split('\t')[1:]]
                meta[meta_map[header]] = the_data
            if line == '!series_matrix_table_begin':
                break

    return meta


def combine_metadata_from_series_matrices(*infiles):
    """
    Combine GEO metadata from multiple series matrices
    :param infiles:
    :return:
    """
    mg = [load_metadata_from_series_matrix(gsf) for gsf in infiles]
    return pd.concat(mg, axis=0, ignore_index=True)


def combine_geo_sra_metadata(geo_series_files, srp_id, add_all_sra=False):
    """
    Create a combined metadata table
    :param geo_series_files: Iterable of .txt or .txt.gz  filenames for series matrices, or a single string.
    will be passed to load_metadata_from_series_matrix
    :param srp_id: SRPxxxxx project ID (or PRJNAxxxx)
    :param add_all_sra: If True, add ALL columns from the SRA runinfo table. Otherwise, just add the SRA
    accessions (default)
    :return:
    """
    from sra import get_runinfo_table
    if not hasattr(geo_series_files, '__iter__'):
        geo_series_files = [geo_series_files]
    meta_geo = combine_metadata_from_series_matrices(*geo_series_files)
    meta_geo.set_index('accession', inplace=True)
    meta_sra = get_runinfo_table(srp_id).dropna(axis=1, how='all')

    # there may be multiple runs for a single sample (e.g. separate lanes)
    # in these cases, duplicate the values in the GEO repo

    g = meta_geo.loc[meta_sra.SampleName]
    g.insert(0, 'accession', g.index)
    g.index = meta_sra.index

    if add_all_sra:
        combined = pd.concat((g, meta_sra), axis=1)
        # rename columns for compatibility
        # cols = combined.columns.tolist()
        # combined.columns = cols
    else:
        combined = g

    combined.index.name = 'filename'
    if 'title' in combined.columns:
        cols = combined.columns.tolist()
        cols[cols.index('title')] = 'sample'
        combined.columns = cols
    return combined
