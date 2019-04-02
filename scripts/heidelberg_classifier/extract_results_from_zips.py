import pandas as pd
import zipfile
import shutil
import os
import numpy as np
import re
from PyPDF2 import PdfFileReader, PdfFileWriter


from utils import log
logger = log.get_console_logger()


def pdf_merger(input_pdf_fns, out_fn):
    """
    Source: https://www.blog.pythonlibrary.org/2018/04/11/splitting-and-merging-pdfs-with-python/
    :param input_pdf_fns: Iterable of paths to pdf files
    :param out_fn: Output filename
    :return:
    """
    pdf_writer = PdfFileWriter()

    for fn in input_pdf_fns:
        pdf_reader = PdfFileReader(fn)
        for page in range(pdf_reader.getNumPages()):
            pdf_writer.addPage(pdf_reader.getPage(page))

    with open(out_fn, 'wb') as fh:
        pdf_writer.write(fh)


def extract_calibrated_scores_from_directory(the_dir):
    cal_scores = {}
    for root_dir, _, filenames in os.walk(the_dir):
        for fn in filenames:
            the_name, the_ext = os.path.splitext(fn)
            if the_ext.lower() == '.zip':
                ff = os.path.join(root_dir, fn)
                with zipfile.ZipFile(ff, mode='r') as z:
                    flist = [x for x in z.namelist() if re.search(r'scores_cal.csv$', x)]
                    if len(flist) != 1:
                        logger.error("Could not find a unique .scores_cal.csv member in archive %s", ff)
                    else:
                        with z.open(flist[0]) as f:
                            dat = pd.read_csv(f, header=0, index_col=0).squeeze().sort_values(ascending=False)
                        cal_scores[the_name] = dat
    return pd.DataFrame(cal_scores)


def extract_report_pdfs_from_directory(the_dir, merged_fn='merged_reports.pdf'):
    """

    :param the_dir:
    :param merged_fn: If supplied, merge all reports in each subdirectory into a single pdf with this filename.
    :return:
    """
    for root_dir, _, filenames in os.walk(the_dir):
        saved = []
        for fn in filenames:
            the_name, the_ext = os.path.splitext(fn)
            if the_ext.lower() == '.zip':
                ff = os.path.join(root_dir, fn)
                with zipfile.ZipFile(ff, mode='r') as z:
                    flist = [
                        x for x in z.namelist() if re.search(r'idat_reportBrain', x) and not re.search(r'_sample_', x)
                    ]
                    if len(flist) != 1:
                        logger.error("Could not find a unique PDF report member in archive %s", ff)
                    else:
                        src = z.extract(flist[0], root_dir)
                        # rename
                        dest = os.path.join(root_dir, "%s_report.pdf" % the_name)
                        shutil.move(src, dest)
                        saved.append(dest)
        if len(saved) > 1 and merged_fn is not None:
            pdf_merger(saved, os.path.join(root_dir, merged_fn))




if __name__ == "__main__":
    """
    Tools for extracting results from the 'full download' zip files served by the Heidelberg classifier.
    Usage example
    dir = '/path/to/heidelberg_classifier/results/'
    cal_scores = extract_calibrated_scores_from_directory(dir)
    cal_scores.to_csv(os.path.join(dir, 'calibrated_scores.csv'))
    extract_report_pdfs_from_directory(dir)
    """
