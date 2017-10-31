import requests
from bs4 import BeautifulSoup
from urllib2 import urlopen
import urlparse
import re
from settings import GIT_LFS_DATA_DIR
import os
import csv


BASE_URL = 'https://www.ncbi.nlm.nih.gov'
SERIES_URL = urlparse.urljoin(BASE_URL, 'geo/query/acc.cgi?acc=GSE28192')
TABLE_HEADER = ['ID_REF', 'VALUE', 'Avg_NBEADS', 'BEAD_STDERR', 'Detection Pval']
OUT_DIR = os.path.join(GIT_LFS_DATA_DIR, 'microarray_GSE28192')

# get download links
html = urlopen(SERIES_URL).read()
soup = BeautifulSoup(html, "lxml")
regex = re.compile('\/geo\/query\/acc\.cgi\?acc=GSM')

# get sample links
soup_samples = soup.find_all('a', href=regex)
samples = {}

for s in soup_samples:
    sname = s.parent.find_next_sibling('td').text
    surl = urlparse.urljoin(BASE_URL, s['href'])
    samples[sname] = surl

# retrieve full data for each
res = {}
for sname, surl in samples.items():
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
            if row[:4] == 'ILMN':
                f.write(row + '\n')

