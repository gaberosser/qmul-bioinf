import requests
from bs4 import BeautifulSoup
from urllib2 import urlopen
import urlparse
import re
from settings import DATA_DIR
import os
import collections
import csv
import pandas as pd

BASE_URL = 'https://www.ncbi.nlm.nih.gov'
SERIES_URL = urlparse.urljoin(BASE_URL, 'geo/query/acc.cgi?acc=GSE37382')
SAMPLE_URL = urlparse.urljoin(BASE_URL, '/geo/query/acc.cgi?acc=')
TABLE_HEADER = ['ID_REF', 'VALUE']
OUT_DIR = os.path.join(DATA_DIR, 'microarray_GSE37382')

TABLE_START_LINE = 22
TABLE_END_LINE = -12

SAMPLES_ALREADY_LOADED = os.listdir(OUT_DIR)

# get download links
html = urlopen(SERIES_URL).read()
soup = BeautifulSoup(html, "lxml")
regex = re.compile('\/geo\/query\/acc\.cgi\?acc=GSM')

# get sample links
soup_samples = soup.find_all('a', href=regex)

samples = {}

for s in soup_samples:
    sname = s.parent.find_next_sibling('td').text.strip('Primary medulloblastoma sample_')
    surl = urlparse.urljoin(BASE_URL, s['href'])
    samples[sname] = surl

errors = {}
sample_info = collections.defaultdict(dict)
for sname, surl in samples.items():
    if sname in SAMPLES_ALREADY_LOADED:
        print "Already d/l sample %s" % sname
        continue
    try:
        gsm_html = urlopen(surl).read()
        soup = BeautifulSoup(gsm_html, "lxml")
        el = soup.find('td', text=re.compile(r'Characteristics')).find_next_sibling()
        for e in el.contents:
            if 'age' in e:
                t = e.strip('age: ')
                sample_info[sname]['age_years'] = float(t) if t != 'N/A' else t
            if 'gender' in e:
                sample_info[sname]['gender'] = e.strip('gender: ')
            if 'subgroup' in e:
                sample_info[sname]['subgroup'] = e.strip('subgroup: ')

        onclick = soup.find('input', value="View full table...")['onclick']
        gsm_url = urlparse.urljoin(
            BASE_URL,
            re.search(r"window.open\('(?P<url>.*)', '_blank'\)", onclick).group('url')
        )
        # retrieve the HTML data for the full table
        this_html = urlopen(gsm_url).read()
        tab = this_html.split('\n')[TABLE_START_LINE:TABLE_END_LINE]
        # process to include only data and save to file
        with open(os.path.join(OUT_DIR, sname), 'wb') as f:
            f.write('\n'.join(tab))
        print "Success: %s" % sname
    except Exception as exc:
        print "Failed: %s" % sname
        errors[sname] = exc

sample_info = dict(sample_info)
for sname, sinfo in sample_info.items():
    sinfo['name'] = sname

index_fn = os.path.join(OUT_DIR, 'sources.csv')
with open(index_fn, 'ab') as f:
    c = csv.DictWriter(f, fieldnames=['name', 'age_years', 'gender', 'subgroup'])
    c.writeheader()
    c.writerows(sample_info.values())


# combine into a single pandas matrix
res = pd.DataFrame()
sample_files = os.listdir(OUT_DIR)
for sname in sample_files:
    if sname[-4:] == '.csv':
        continue
    s = pd.read_csv(os.path.join(OUT_DIR, sname), header=None, index_col=0, sep='\t').loc[:, 1]
    res[sname] = s
