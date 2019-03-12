import os
from settings import OUTPUT_DIR
import requests
import time
import glob
import datetime
from utils import log
from Bio.UniProt.GOA import gafiterator
from goatools import base, obo_parser, associations, go_enrichment
import pandas as pd
import urllib
import gzip
from StringIO import StringIO


DEFAULT_GO_DIR = os.path.join(OUTPUT_DIR, 'go_analysis')
datestr_fmt = '%Y-%m-%d'


def get_gaf_file(retries=3):
    url = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
    i = 0
    dat = None
    while True:
        try:
            dat = urllib.urlopen(url).read()
        except Exception:
            if i < retries:
                time.sleep(1)
            else:
                raise
        else:
            break
        i += 1
    return gzip.GzipFile(fileobj=StringIO(dat)).read()


def get_obo_file(retries=3):
    url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    i = 0
    dat = None
    while i < retries:
        try:
            resp = requests.get(url)
            resp.raise_for_status()
        except Exception:
            if i < retries:
                time.sleep(1)
            else:
                raise
        else:
            break
        i += 1

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


def ens_to_entrez(ens, genetoens_fn=None):
    retries = 3
    if genetoens_fn is None:
        url = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz"
        i = 0
        dat = None
        while True:
            try:
                dat = urllib.urlopen(url).read()
            except Exception:
                if i < retries:
                    time.sleep(1)
                else:
                    raise
            else:
                break
            i += 1
        genetoens_fn = gzip.GzipFile(fileobj=StringIO(dat))

    entrezid_to_ens = pd.read_csv(genetoens_fn, sep='\t', index_col=None, header=0)
    entrezid_to_ens = entrezid_to_ens.loc[entrezid_to_ens['#tax_id'] == 9606, ['GeneID', 'Ensembl_gene_identifier']]
    entrezid_to_ens.columns = ['entrez_id', 'ensembl_id']
    entrezid_to_ens = entrezid_to_ens.loc[~entrezid_to_ens.entrez_id.duplicated()].set_index('entrez_id')

    ens_to_entrezid = entrezid_to_ens.copy()
    ens_to_entrezid.insert(0, 'entrez_id', ens_to_entrezid.index)
    ens_to_entrezid = ens_to_entrezid.loc[~ens_to_entrezid.ensembl_id.duplicated()]
    ens_to_entrezid.set_index('ensembl_id', inplace=True)

    return ens_to_entrezid.reindex(ens).dropna().astype(int)


class GOAnalysis(object):

    def __init__(self, tax_id=9606, logger=None, force_update=False, go_dir=DEFAULT_GO_DIR, bg_genes=None):
        # gene_converter can be used to enable automatic gene conversion
        self.gene_converter = None
        self.logger = logger or log.get_console_logger(self.__class__.__name__)
        self.tax_id = tax_id
        if not os.path.isdir(go_dir):
            self.logger.warn("Creating master GO directory at %s.", go_dir)
            os.makedirs(go_dir)
        else:
            self.logger.info("Using existing GO directory at %s.", go_dir)
        self.base_dir = go_dir

        # get filenames and parse both GAF and OBO
        self.obo_fn = self.check_and_get_obo(force_update=force_update)
        self.gaf_fn = self.check_and_get_gaf(force_update=force_update)
        self.obo = obo_parser.GODag(self.obo_fn)

        self.gaf = associations.read_ncbi_gene2go(self.gaf_fn, taxids=[self.tax_id])
        self.logger.info("{N:,} annotated human genes".format(N=len(self.gaf)))

        self.bg_genes = bg_genes
        if self.bg_genes is not None:
            self.set_bg_genes(bg_genes)


    def check_and_get_obo(self, force_update=False):
        root_dir = os.path.join(self.base_dir, "obo")
        return check_and_get_file(root_dir, 'obo', get_obo_file, force_update=force_update, logger=self.logger)

    def check_and_get_gaf(self, force_update=False):
        root_dir = os.path.join(self.base_dir, "gaf")
        return check_and_get_file(root_dir, 'gaf', get_gaf_file, force_update=force_update, logger=self.logger)

    def set_bg_genes(self, genes):
        if self.gene_converter is not None:
            genes = self._convert(genes)
        self.bg_genes = genes
        self.initialise_goea()

    def initialise_goea(self):
        self.goe_obj = go_enrichment.GOEnrichmentStudy(
            self.bg_genes,
            self.gaf,
            self.obo,
        )

    def run_one(self, genes, tabulate=True):
        if self.gene_converter is not None:
            genes = self._convert(genes)
        res = self.goe_obj.run_study(genes)
        if self.gene_converter is not None:
            # convert the study_items entry back
            for t in res:
                # if statement not required, but may be faster here?
                if len(t.study_items) > 0:
                    t.study_items = self.gene_converter.index[self.gene_converter.isin(t.study_items)].tolist()
                else:
                    t.study_items = []
        if tabulate:
            prt_flds = self.goe_obj.get_prtflds_default(res)
            res = pd.DataFrame(go_enrichment.get_goea_nts_prt(res, prt_flds))
            res.set_index('GO', inplace=True)
        return res

    def set_gene_conversion(self, series):
        """
        Set the gene conversion dataframe. This will be used for all runs.
        For example, human GO inputs expect Entrez gene IDs, but we might wish to convert from Ensembl IDs.
        :param series: pd.Series indexed by the expected input and with values compatible with the default GO IDs.
        :return:
        """
        series = series.squeeze()
        # sanity check: is there any overlap between the GAF and the values?
        intersct = pd.Index(series.values).intersection(self.gaf.keys())
        if len(intersct) == 0:
            raise ValueError("Incompatible gene conversion series supplied: no overlap with currently loaded GAF.")
        self.gene_converter = series

    def _convert(self, genes):
        overlap = self.gene_converter.index.intersection(genes)
        if len(overlap) != len(genes):
            self.logger.warning(
                "Dropped %d of the %d supplied genes because they could not be converted.",
                len(genes) - len(overlap),
                len(genes)
            )
        return self.gene_converter.loc[overlap].values