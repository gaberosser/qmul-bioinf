import ftplib
import urllib2
from utils.output import unique_output_dir
from utils.log import get_console_logger
import os
logger = get_console_logger(__name__)



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
