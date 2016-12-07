import pandas as pd
import re
import os

infile = '/home/gabriel/Downloads/GSE50765_series_matrix.txt'

meta_map = {
    'Sample_title': 'title',
    'Sample_geo_accession': 'accession',
    'Sample_characteristics_subgroup': 'subgroup',
    'Sample_characteristics_age': 'age_years',
    'Sample_characteristics_histology': 'histology',
    'Sample_characteristics_followup': 'followup',
    'Sample_characteristics_death': 'death',
}
# meta = collections.defaultdict(list)
meta = pd.DataFrame()

# raw = []
# gene_ids = []
with open(infile, 'rb') as f:
    while True:
        line = f.readline().strip('\n')
        if len(line) == 0:
            continue
        header = re.match(r'^!(?P<hd>[^\t]*)', line).group('hd')
        if header in meta_map:
            meta[meta_map[header]] = [t.strip('"') for t in line.split('\t')[1:]]
        if line == '!series_matrix_table_begin':
            break

# set the index of meta
meta.set_index('title', inplace=True)
meta.loc[:, 'age_years'] = meta.loc[:, 'age_years'].astype(float)
meta.loc[meta.death == 'N/A', 'death'] = None

# age_str = meta.loc[:, 'age_months']
# intage = age_str.apply(lambda x: 12 * int(re.sub(r'(?P<y>[0-9]*)yrs.*', r'\g<y>', x)) + int(re.sub(r'[0-9]*yrs (?P<m>[0-9]*)mos', r'\g<m>', x)))
# meta.loc[:, 'age_months'] = intage

meta.to_csv("sources.csv", sep=',')
# standardise title formats
# titles = []
# for t in meta.index:
#     tn = re.sub('^[iI][cC](?P<i>[0-9])', 'ICb\g<i>', t)
#     tn = re.sub('^[iI][cC][bB]', 'ICb', tn)
#     tn = re.sub('^[pP][tT]', 'Pt', tn)
#     tn = re.sub(r'(?P<a>Pt|ICb|NT|NCb)-(?P<b>[0-9]*)', '\g<a>\g<b>', tn)
#     tn = re.sub('(?P<a>-i+)', lambda x: x.group('a').upper(), tn)
#     tn = re.sub(r'MB', '', tn)
#     titles.append(tn)

# meta.index = titles
