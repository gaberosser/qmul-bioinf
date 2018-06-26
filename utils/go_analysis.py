import os
from settings import OUTPUT_DIR
import requests
import urlparse
import gzip
import glob
import datetime
from utils import log
from Bio.UniProt.GOA import gafiterator
from goatools import obo_parser, go_enrichment


go_dir = os.path.join(OUTPUT_DIR, 'go_analysis')
datestr_fmt = '%Y-%m-%d'


def get_gaf_file(tax_id=9606):
    root_url = "http://geneontology.org/gene-associations/"
    link_dict = {
        9606: 'goa_human.gaf.gz',
        10090: 'gene_association.mgi.gz'
    }
    url = urlparse.urljoin(root_url, link_dict[tax_id])
    resp = requests.get(url)
    resp.raise_for_status()
    return resp.content


def get_obo_file():
    url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    resp = requests.get(url)
    resp.raise_for_status()
    return resp.content


def check_and_get_file(root_dir, ext, get_func, force_update=False, logger=None):
    """
    Check for a file with a given extension in the supplied dir.
    :param root_dir:
    :param ext:
    :param get_fn: Function handle that, when called with no args, fetches and returns the data
    :param force_update: If True, don't use existing file but instead force a new get call.
    :return: String, giving the filename of the data
    """
    if logger is None:
        logger = log.get_console_logger("check_and_get_file")

    if not os.path.isdir(root_dir):
        logger.warn("Creating %s directory at %s.", ext, root_dir)
        os.makedirs(root_dir)
    flist = glob.glob(os.path.join(root_dir, "*.%s" % ext))
    files_seen = {}
    for full_fn in flist:
        fn = os.path.split(full_fn)[1]
        try:
            d = datetime.datetime.strptime(fn, "%s.%s" % (datestr_fmt, ext))
        except ValueError:
            logger.warn("Failed to parse version of %s file %s. Skipping.", ext, full_fn)
        else:
            files_seen[d] = full_fn
    if force_update or len(files_seen) == 0:
        fn_out = os.path.join(root_dir, "%s.%s" % (datetime.date.today().strftime(datestr_fmt), ext))
        dat = get_func()
        with open(fn_out, 'wb') as fout:
            fout.write(dat)
        logger.info("Downloaded new %s file and saved it at %s", ext, fn_out)
    else:
        latest_date = max(files_seen.keys())
        logger.info("Using existing %s file %s.", ext, files_seen[latest_date])
        fn_out = files_seen[latest_date]

    return fn_out


def parse_gaf(fn):
    if os.path.splitext(fn)[1].lower()[-3:] == '.gz':
        open_func = gzip.open
    else:
        open_func = open
    with open_func(fn, 'rb') as f:
        it = gafiterator(f)
        gene_by_go_term = {}
        go_term_by_gene = {}
        for rec in it:
            gs = rec['DB_Object_Symbol']
            go_id = rec['GO_ID']
            if go_id not in gene_by_go_term:
                gene_by_go_term[go_id] = []
            if gs not in go_term_by_gene:
                go_term_by_gene[gs] = []
            gene_by_go_term[go_id].append(gs)
            go_term_by_gene[gs].append(go_id)
    return gene_by_go_term, go_term_by_gene


def parse_obo(fn):
    f = obo_parser.OBOReader(fn)
    res = list(f)
    # return res
    # go_enrichment.GOEnrichmentStudy()
    # export as a dict
    out = {}
    for r in res:
        if r.id in out:
            print "WARNING: already seen ID %s" % r.id
        out[r.id] = r
    return out


class GOAnalysis(object):

    def __init__(self, tax_id=9606, logger=None, force_update=False):
        self.logger = logger or log.get_console_logger(self.__class__.__name__)
        self.tax_id = tax_id
        if not os.path.isdir(go_dir):
            self.logger.warn("Creating master GO directory at %s.", go_dir)
            os.makedirs(go_dir)
        else:
            self.logger.info("Using existing GO directory at %s.", go_dir)

        # get filenames and parse both GAF and OBO
        self.obo_fn = self.check_and_get_obo(force_update=force_update)
        self.gaf_fn = self.check_and_get_gaf(force_update=force_update)
        self.feat_by_go, self.go_by_feat = parse_gaf(self.gaf_fn)
        self.obo = parse_obo(self.obo_fn)

    def check_and_get_obo(self, force_update=False):
        root_dir = os.path.join(go_dir, "obo")
        return check_and_get_file(root_dir, 'obo', get_obo_file, force_update=force_update, logger=self.logger)

    def check_and_get_gaf(self, force_update=False):
        root_dir = os.path.join(go_dir, "gaf", str(self.tax_id))
        return check_and_get_file(root_dir, 'gz', get_gaf_file, force_update=force_update, logger=self.logger)




