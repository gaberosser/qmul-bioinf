import requests
import os
from bs4 import BeautifulSoup

from settings import DATA_DIR_NON_GIT
from log import get_console_logger
logger = get_console_logger(__name__)


pw = 'Nucleus1'
email = 'g.rosser@qmul.ac.uk'
the_login_url = 'https://www.molecularneuropathology.org/mnp/authenticate'
the_sample_url = 'https://www.molecularneuropathology.org/mnp/sample/{sid}'
the_report_url = 'https://www.molecularneuropathology.org/mnp/sample/{sid}/run/{rid}/report'
the_analysis_results_url = 'https://www.molecularneuropathology.org/mnp/sample/{sid}/run/{rid}/analysisData'


def read_table(tbl):
    res = []
    for row in tbl.findAll('tr'):
        this_row = []
        for col in row.findAll(['td', 'th']):
            this_row.append(col.text)
        res.append(this_row)
    return res

# sample_ids = [784, 785] + range(787, 800) + [801]
# outdir = os.path.join(
#     DATA_DIR_NON_GIT,
#     'methylation',
#     '2016-06-10_brandner',
#     'heidelberg_classifier',
#     '2017-01'
# )

# sample_ids = range(768, 784) + [701, 702]
# outdir = os.path.join(
#     DATA_DIR_NON_GIT,
#     'methylation',
#     '2016-12-19_ucl_genomics',
#     'heidelberg_classifier',
#     '2017-01'
# )

sample_ids = range(802, 810)
outdir = os.path.join(
    DATA_DIR_NON_GIT,
    'methylation',
    '2016-09-21_dutt',
    'heidelberg_classifier',
    '2017-01'
)

if not os.path.exists(outdir):
    logger.info("Creating output directory %s", outdir)
    os.makedirs(outdir)
else:
    logger.info("Using existing output directory %s", outdir)


s = requests.Session()
s.post(the_login_url, data={"email": "g.rosser@qmul.ac.uk", "password": pw})

for sid in sample_ids:
    try:
        the_url = the_sample_url.format(sid=sid)
        resp = s.get(the_url)
        soup = BeautifulSoup(resp.content, "html.parser")

        # get the report table (it's the 1st table)
        tbl = soup.find('table')
        tbl_cont = dict(read_table(tbl))
        sname = tbl_cont['Sample identifier:']

        # get the most recent run number (it's in the 2nd table)
        # FIXME: we're assuming that this table has only one row, will break otherwise
        tbl = soup.findAll('table')[1]
        tbl_cont = dict(zip(*read_table(tbl)))
        rid = int(tbl_cont['Run'])

        # download the pdf report
        the_url = the_report_url.format(sid=sid, rid=rid)
        resp = s.get(the_url)

        outfile = os.path.join(outdir, "%s.pdf" % sname)
        logger.info("Saving PDF file to %s", outfile)
        with open(outfile, 'wb') as f:
            f.write(resp.content)

        # download the full analysis results
        the_url = the_analysis_results_url.format(sid=sid, rid=rid)
        resp = s.get(the_url)

        outfile = os.path.join(outdir, "%s.zip" % sname)
        logger.info("Saving zip file to %s", outfile)
        with open(outfile, 'wb') as f:
            f.write(resp.content)

    except Exception as exc:
        logger.exception("Sample ID %d failed", sid)
        print repr(exc)