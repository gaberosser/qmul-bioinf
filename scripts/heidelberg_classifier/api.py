import requests
from utils.log import get_console_logger
import os
import re
from settings import HEIDELBERG_CLASSIFIER_CONFIG
from bs4 import BeautifulSoup
import pandas as pd

logger = get_console_logger(__name__)


DEFAULT_PAYLOAD = {
    'age': '',
    'location': '',
    'diagnosis': '',
    'notes': '',
    'gender': 'NA',
    'chipType': 'NA',
    'sampleType': 'KRYO DNA',
    'confirmation': 'true',
    'keepFiles': 'false',
}


def read_samplesheet(fn):
    with open(fn, 'rb') as f:
        n_skip = None
        for i, line in enumerate(f):
            if re.match(r'\[Data\]', line):
                n_skip = i + 1
                break
    if n_skip is None:
        raise AttributeError("Unable to find the start of the data block in sample sheet %s." % fn)
    return pd.read_csv(fn, header=0, index_col=None, skiprows=n_skip)


def paired_idat_file(one_file):
    """
    Given the supplied path, find the paired idat file.
    :param one_file:
    :return:
    """
    if not os.path.isfile(one_file):
        raise IOError("Unable to find the input file %s" % one_file)
    patterns = (
        (r'_Grn\.idat', '_Red.idat'),
        (r'_Red\.idat', '_Grn.idat'),
    )
    for patt, repl in patterns:
        other = re.sub(patt, repl, one_file)
        if other != one_file:
            if os.path.isfile(other):
                return other
    raise AttributeError("Unable to find paired file for input file %s" % one_file)


def read_table(tbl):
    res = []
    for row in tbl.findAll('tr'):
        this_row = []
        for col in row.findAll(['td', 'th']):
            this_row.append(col.text)
        res.append(this_row)
    return res


class HeidelbergClassifier(object):
    LOGIN_URL = 'https://www.molecularneuropathology.org/mnp/authenticate'
    SUBMISSION_URL = 'https://www.molecularneuropathology.org/mnp/sample/add'
    SAMPLE_URL = 'https://www.molecularneuropathology.org/mnp/sample/{sid}'
    REPORT_URL = 'https://www.molecularneuropathology.org/mnp/sample/{sid}/run/{rid}/report'
    ANALYSIS_RESULTS_URL = 'https://www.molecularneuropathology.org/mnp/sample/{sid}/run/{rid}/analysisData'

    def __init__(self, chiptype=None, sampletype=None):
        self.default_payload = dict(DEFAULT_PAYLOAD)
        if chiptype is not None:
            self.default_payload['chipType'] = chiptype
        if sampletype is not None:
            self.default_payload['sampleType'] = sampletype
        self.session = None
        self.establish_session()
        self.submitted = []

    def establish_session(self):
        self.session = requests.Session()
        resp = self.session.post(self.LOGIN_URL, data=HEIDELBERG_CLASSIFIER_CONFIG)
        if resp.status_code != 200:
            raise AttributeError("Login failed")

    def submit_sample(self, name, file, **payload_kwargs):
        """
        Submit a sample for classification
        :param name: Sample name
        :param file: Path to either the red or the green file. The other file is inferred from this name.
        :param payload_kwargs: Other payload kwargs. Defaults are used when these are missing - see DEFAULT_PAYLOAD.
        :return: Response object
        """
        file2 = paired_idat_file(file)

        with open(file, 'rb') as f1, open(file2, 'rb') as f2:
            files = {
                'file1': f1,
                'file2': f2
            }
            payload = dict(payload_kwargs)
            payload['name'] = name
            for f, v in self.default_payload.items():
                payload.setdefault(f, v)
            resp = self.session.post(self.SUBMISSION_URL, data=payload, files=files)
        self.submitted.append({'payload': payload, 'files': files, 'response': resp, 'response_code': resp.status_code})
        if resp.status_code != 200:
            raise AttributeError("Upload failed for sample %s: %s" % (name, resp.content))
        return resp

    def get_result(self, sample_id, outdir, filestem=None):
        # TODO: test this.
        # TODO: Can we get the sample_id automatically from self.submitted?
        # TODO: Support optional filestem argument, required when names themselves are too long / contain forbidden
        # characters
        the_url = self.SAMPLE_URL.format(sid=sample_id)
        resp = self.session.get(the_url)
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
        the_url = self.REPORT_URL.format(sid=sample_id, rid=rid)
        resp = self.session.get(the_url)

        outfile = os.path.join(outdir, "%s.pdf" % sname)
        logger.info("Saving PDF file to %s", outfile)
        with open(outfile, 'wb') as f:
            f.write(resp.content)

        # download the full analysis results
        the_url = self.ANALYSIS_RESULTS_URL.format(sid=sample_id, rid=rid)
        resp = self.session.get(the_url)

        outfile = os.path.join(outdir, "%s.zip" % sname)
        logger.info("Saving zip file to %s", outfile)
        with open(outfile, 'wb') as f:
            f.write(resp.content)