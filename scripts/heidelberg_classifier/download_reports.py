import os

import requests
from bs4 import BeautifulSoup
from scripts.heidelberg_classifier import api
from settings import DATA_DIR_NON_GIT, HEIDELBERG_CLASSIFIER_CONFIG
from utils.log import get_console_logger

logger = get_console_logger(__name__)


# the_login_url = 'https://www.molecularneuropathology.org/mnp/authenticate'
# the_sample_url = 'https://www.molecularneuropathology.org/mnp/sample/{sid}'
# the_report_url = 'https://www.molecularneuropathology.org/mnp/sample/{sid}/run/{rid}/report'
# the_analysis_results_url = 'https://www.molecularneuropathology.org/mnp/sample/{sid}/run/{rid}/analysisData'


# def read_table(tbl):
#     res = []
#     for row in tbl.findAll('tr'):
#         this_row = []
#         for col in row.findAll(['td', 'th']):
#             this_row.append(col.text)
#         res.append(this_row)
#     return res

# sample_ids = range(2248, 2298)
sample_ids = range(2285, 2298)
outdir = os.path.join(
    DATA_DIR_NON_GIT,
    'methylation',
    'tcga_gbm',
    'heidelberg_classifier',
    '2017-03'
)

if __name__ == "__main__":
    sample_ids = range(6857, 6977)
    # sample_ids = range(6858, 6977)
    obj = api.Heidelberg()
    res = []
    success = []
    error = []
    reset = []
    for sid in sample_ids:
        try:
            this_res = obj.get_result(sid)
            res.append(this_res)
            if this_res is None:
                reset.append(sid)
            success.append(sid)
        except Exception as exc:
            print "Failed to retrieve sample %d: %s" % (sid, repr(exc))
            error.append(sid)




# if not os.path.exists(outdir):
#     logger.info("Creating output directory %s", outdir)
#     os.makedirs(outdir)
# else:
#     logger.info("Using existing output directory %s", outdir)


# s = requests.Session()
# s.post(the_login_url, data=HEIDELBERG_CLASSIFIER_CONFIG)
#
# for sid in sample_ids:
#     try:
#         the_url = the_sample_url.format(sid=sid)
#         resp = s.get(the_url)
#         soup = BeautifulSoup(resp.content, "html.parser")
#
#         # get the report table (it's the 1st table)
#         tbl = soup.find('table')
#         tbl_cont = dict(read_table(tbl))
#         sname = tbl_cont['Sample identifier:']
#
#         # get the most recent run number (it's in the 2nd table)
#         # FIXME: we're assuming that this table has only one row, will break otherwise
#         tbl = soup.findAll('table')[1]
#         tbl_cont = dict(zip(*read_table(tbl)))
#         rid = int(tbl_cont['Run'])
#
#         # download the pdf report
#         the_url = the_report_url.format(sid=sid, rid=rid)
#         resp = s.get(the_url)
#
#         outfile = os.path.join(outdir, "%s.pdf" % sname)
#         logger.info("Saving PDF file to %s", outfile)
#         with open(outfile, 'wb') as f:
#             f.write(resp.content)
#
#         # download the full analysis results
#         the_url = the_analysis_results_url.format(sid=sid, rid=rid)
#         resp = s.get(the_url)
#
#         outfile = os.path.join(outdir, "%s.zip" % sname)
#         logger.info("Saving zip file to %s", outfile)
#         with open(outfile, 'wb') as f:
#             f.write(resp.content)
#
#     except Exception as exc:
#         logger.exception("Sample ID %d failed", sid)
#         print repr(exc)
