import os
from settings import OUTPUT_DIR
import requests
import gzip
import glob
import datetime
from utils import log


go_dir = os.path.join(OUTPUT_DIR, 'go_analysis')
datestr_fmt = '%Y-%m-%d'


def get_gaf_file(tax_id=9606):
    #TODO
    pass


def get_obo_file():
    url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    resp = requests.get(url)
    resp.raise_for_status()
    return resp.content


class GOAnalysis(object):

    def __init__(self, tax_id=9606, logger=None, force_obo_update=False):
        self.logger = logger or log.get_console_logger(self.__name__)
        self.tax_id = tax_id
        if not os.path.isdir(go_dir):
            self.logger.warn("Creating master GO directory at %s.", go_dir)
            os.makedirs(go_dir)
        else:
            self.logger.info("Using existing GO directory at %s.", go_dir)
        self.check_and_get_obo(force_update=force_obo_update)
        self.check_and_get_gaf()

    def check_and_get_obo(self, force_update=False):
        root_dir = os.path.join(go_dir, "obo")
        if not os.path.isdir(root_dir):
            self.logger.warn("Creating OBO directory at %s.", root_dir)
            os.makedirs(root_dir)
        flist = glob.glob(os.path.join(root_dir, "*.obo"))
        files_seen = {}
        for full_fn in flist:
            fn = os.path.split(full_fn[1])
            try:
                d = datetime.datetime.strptime(fn, "%s.obo" % datestr_fmt)
            except ValueError:
                self.logger.warn("Failed to parse version of OBO file %s. Skipping.", full_fn)
            else:
                files_seen[d] = full_fn
        if force_update or len(files_seen) == 0:
            fn_out = os.path.join(root_dir, "%s.obo" % datetime.date.today().strftime(datestr_fmt))
            obo_dat = get_obo_file()
            with open(fn_out, 'wb') as fout:
                fout.write(obo_dat)
            self.logger.info("Downloaded new OBO file and saved it at %s", fn_out)
            self.obo_fn = fn_out
        else:
            latest_date = max(files_seen.keys())
            self.logger.info("Using existing OBO file %s.", files_seen[latest_date])
            self.obo_fn = files_seen[latest_date]

    def check_and_get_gaf(self):
        pass




