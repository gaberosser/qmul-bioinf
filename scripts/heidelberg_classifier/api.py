import requests
from utils.log import get_console_logger
import os
import re
import datetime
from utils.output import unique_output_dir
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


class Heidelberg(object):
    ROOT_URL = 'https://www.molecularneuropathology.org'
    LOGIN_URL = 'https://www.molecularneuropathology.org/mnp/authenticate'
    SUBMISSION_URL = 'https://www.molecularneuropathology.org/mnp/sample/add'
    SAMPLE_URL = 'https://www.molecularneuropathology.org/mnp/sample/{sid}'
    REPORT_URL = 'https://www.molecularneuropathology.org/mnp/sample/{sid}/run/{rid}/report'
    ANALYSIS_RESULTS_URL = 'https://www.molecularneuropathology.org/mnp/sample/{sid}/run/{rid}/analysisData'
    RESTART_ANALYSIS_URL = 'https://www.molecularneuropathology.org/mnp/sample/{sid}/{rid}/'

    def __init__(self):
        self.session = None
        self.establish_session()
        self.submitted = []
        self.outdir = None

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
        payload = dict(DEFAULT_PAYLOAD)
        payload.update(payload_kwargs)
        payload['name'] = name

        with open(file, 'rb') as f1, open(file2, 'rb') as f2:
            files = {
                'file1': f1,
                'file2': f2
            }
            resp = self.session.post(self.SUBMISSION_URL, data=payload, files=files)
        self.submitted.append(
            {'payload': payload, 'files': files, 'response': resp, 'response_code': resp.status_code})
        if resp.status_code != 200:
            raise AttributeError("Upload failed for sample %s: %s" % (name, resp.content))
        return resp

    def get_result_soup(self, sample_id):
        the_url = self.SAMPLE_URL.format(sid=sample_id)
        resp = self.session.get(the_url)
        if resp.status_code != 200:
            raise AttributeError("Error %d accessing %s: %s" % (resp.status_code, the_url, resp.content))
        return BeautifulSoup(resp.content, "html.parser")

    def get_summary_data(self, sample_id=None, soup=None):
        """
        :param sample_id: If supplied, we use this to retrieve the relevant results page
        :param soup: If supplied, sample_id is ignored and we use the BeautifulSoup data here to get our info.
        """
        if soup is None:
            if sample_id is None:
                raise AttributeError("Must supply either soup or sample_id")
            soup = self.get_result_soup(sample_id)

        # get the summary table (it's the 1st table)
        tbl = soup.find('table')
        if tbl is None:
            raise ValueError("Unable to find any table elements")
        tbl_cont = dict(read_table(tbl))
        batch, sentrix_id, sname = tbl_cont['Sample identifier:'].split(';')

        # get the most recent run number and creation time (it's in the 2nd table)
        # FIXME: we're assuming that this table has only one row, will break otherwise
        tbl = soup.findAll('table')[1]
        tbl_cont = dict(zip(*read_table(tbl)))

        run_id = int(tbl_cont['Run'])
        created_at = datetime.datetime.strptime(tbl_cont['Created at '], "%Y-%m-%d %H:%M")

        res = {
            'batch': batch,
            'run_id': run_id,
            'sentrix_id': sentrix_id,
            'sample_name': sname,
            'created_at': created_at,
        }
        if sample_id is not None:
            res['soup'] = soup

        return res

    def restart_analysis(self, sample_id, run_id=None, soup=None):
        """
        Restart a previously submitted analysis run.
        This is required if it has stalled
        """
        # FIXME: POST request returns a 404 error - why?
        if run_id is None:
            if soup is None:
                soup = self.get_result_soup(sample_id)
            summary = self.get_summary_data(soup=soup)
            run_id = summary['run_id']

        the_url = self.RESTART_ANALYSIS_URL.format(sid=sample_id, rid=run_id)
        return self.session.post(the_url)

    def get_result(self, sample_id, outdir=None, sample_name=None, run_id=None):
        """
        Retrieve results relating to a sample and save to disk.
        :param sample_name: If supplied, this overrides the submitted sample name
        :param run_id: If supplied, this is used, otherwise the latest run is automatically determined
        """

        if self.outdir is None:
            if outdir is None:
                self.outdir = unique_output_dir('heidelberg_classifier', reuse_empty=True)
            else:
                self.outdir = outdir
            print "Data will be downloaded to %s" % self.outdir

        the_url = self.SAMPLE_URL.format(sid=sample_id)
        resp = self.session.get(the_url)
        soup = BeautifulSoup(resp.content, "html.parser")

        summary = self.get_summary_data(soup=soup)
        batch = summary['batch']
        if sample_name is None:
            sample_name = summary['sample_name'].strip()
        # ensure sample name is a valid identifier
        sample_name = re.sub(r' +', '_', sample_name)
        sample_name = sample_name.replace('/', '-')
        sample_name = sample_name.replace('\\','-')

        if run_id is None:
            run_id = summary['run_id']
        created_at = summary['created_at']

        logger.info("Sample %s, run ID %d, batch %s", sample_name, run_id, batch)

        # one of three situations:
        # 1) Classifier has not finished any modules. Probably needs restarting.
        # 2) Classification has completed but full report not available. Retrieve classification scores.
        # 3) Full report available. Download all data.

        t = soup.findAll(text=re.compile(r'.*Classifier script not finished.*'))
        if len(t) > 0:
            # situation (1)
            # get creation time
            # this is fragile, but easier than trawlind through tables!
            logger.info("Sample ID %d (%s). Classification script is not finished. Nothing to do.", sample_id, sample_name)
            dt = (datetime.datetime.utcnow() - created_at).total_seconds()
            if dt > 18000:
                logger.warn("Submitted more than 5 hours ago. Consider restarting.")
            return

        # create the output subdir if necessary
        out_subdir = os.path.join(self.outdir, batch)
        if not os.path.isdir(out_subdir):
            try:
                os.makedirs(out_subdir)
            except OSError as exc:
                logger.error("Failed to create output directory %s", out_subdir)


        # try getting the pdf report
        aa = soup.findAll('a')
        aa = [t for t in aa if re.search(r'Download *idat_', t.get_text())]
        if len(aa) == 0:
            logger.error("No download link found")
        else:
            the_url = '/'.join([self.ROOT_URL, aa[0]['href']])

            resp = self.session.get(the_url)
            if resp.status_code == 200:
                # if this works, we know we're in situation (3)
                outfile = os.path.join(self.outdir, batch, "%s.pdf" % sample_name)
                if os.path.isfile(outfile):
                    logger.error("File already exists: %s", outfile)
                logger.info("Saving PDF file to %s", outfile)
                with open(outfile, 'wb') as f:
                    f.write(resp.content)

            # download the full analysis results
            the_url = self.ANALYSIS_RESULTS_URL.format(sid=sample_id, rid=run_id)
            logger.info("Downloading zipped results file for sample %s", sample_id)
            resp = self.session.get(the_url)

            outfile = os.path.join(self.outdir, batch, "%s.zip" % sample_name)
            if os.path.isfile(outfile):
                logger.error("File already exists: %s", outfile)
            logger.info("Saving zip file to %s", outfile)
            with open(outfile, 'wb') as f:
                f.write(resp.content)

        # situation (2) OR (3)
        # Either way, get the classifier results
        try:
            raw_scores = read_table(soup.find(attrs={'id': 'rawScores'}))
            raw_scores = self.save_scores(raw_scores, sample_name, batch, 'raw_scores')
        except Exception:
            logger.exception("Failed to retrieve raw scores.")

        try:
            cal_scores = read_table(soup.find(attrs={'id': 'calibratedScores'}))
            cal_scores = self.save_scores(cal_scores, sample_name, batch, 'calibrated_scores')
        except Exception:
            logger.exception("Failed to retrieve calibrated scores.")

        ## TODO: return filename(s) where result(s) are stored

    def save_scores(self, table, sample_name, batch, typ):
        scores = pd.DataFrame(table)
        scores.columns = scores.iloc[0].str.strip()
        scores = scores.iloc[1:]
        scores.sort_values('Score', ascending=False, inplace=True)
        outfile = os.path.join(self.outdir, batch, "%s.%s.csv" % (sample_name, typ))

        if os.path.isfile(outfile):
            logger.error("File already exists: %s", outfile)

        scores.to_csv(outfile, index=False)
        print "Saved scores to %s" % outfile
        return scores



class HeidelbergUploader(object):
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
