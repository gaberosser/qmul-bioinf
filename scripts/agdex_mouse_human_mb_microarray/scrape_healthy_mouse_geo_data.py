import requests
from bs4 import BeautifulSoup
from urllib2 import urlopen
import urlparse
import re
from settings import GIT_LFS_DATA_DIR
import os

GSM_LIST = {
    "GSM1321086": "Cer_CT18",
    "GSM1321087": "Cer_CT20",
    "GSM1321088": "Cer_CT22",
    "GSM1321089": "Cer_CT24",
    "GSM1321090": "Cer_CT26",
    "GSM1321091": "Cer_CT28",
    "GSM1321092": "Cer_CT30",
    "GSM1321093": "Cer_CT32",
    "GSM1321094": "Cer_CT34",
    "GSM1321095": "Cer_CT36",
    "GSM1321096": "Cer_CT38",
    "GSM1321097": "Cer_CT40",
    "GSM1321098": "Cer_CT42",
    "GSM1321099": "Cer_CT44",
    "GSM1321100": "Cer_CT46",
    "GSM1321101": "Cer_CT48",
    "GSM1321102": "Cer_CT50",
    "GSM1321103": "Cer_CT52",
    "GSM1321104": "Cer_CT54",
    "GSM1321105": "Cer_CT56",
    "GSM1321106": "Cer_CT58",
    "GSM1321107": "Cer_CT60",
    "GSM1321108": "Cer_CT62",
    "GSM1321109": "Cer_CT64",
}

BASE_URL = 'https://www.ncbi.nlm.nih.gov'
SERIES_URL = urlparse.urljoin(BASE_URL, 'geo/query/acc.cgi')
TABLE_HEADER = ['ID_REF', 'VALUE']
OUT_DIR = os.path.join(GIT_LFS_DATA_DIR, 'microarray_GSE54650')

if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

# retrieve full data for each
res = {}
for acc, sname in GSM_LIST.items():
    surl = SERIES_URL + '?acc=%s' % acc
    gsm_html = urlopen(surl).read()
    soup = BeautifulSoup(gsm_html, "lxml")
    onclick = soup.find('input', value="View full table...")['onclick']
    gsm_url = urlparse.urljoin(
        BASE_URL,
        re.search(r"window.open\('(?P<url>.*)', '_blank'\)", onclick).group('url')
    )
    # retrieve the HTML data for the full table
    this_html = urlopen(gsm_url).read()
    # process to include only data and save to file
    with open(os.path.join(OUT_DIR, sname), 'wb') as f:
        f.write('\t'.join(TABLE_HEADER) + '\n')
        for row in this_html.split('\n'):
            if re.match(r'[0-9]', row):
                f.write(row + '\n')
