import requests
import os
import json
import collections
import csv
import pandas as pd
from utils.output import unique_output_dir


base_url = "https://www.encodeproject.org/"
expt_url = "https://www.encodeproject.org/experiments/{eid}/"
expt_ids = [
    # long RNA-Seq
    'ENCSR490SQH',
    'ENCSR950PSB',
    'ENCSR000EYP',
    'ENCSR000COU',
    'ENCSR244ISQ',
    'ENCSR977XUX',
    'ENCSR572EET',
    'ENCSR670WQY',
    'ENCSR043RSE',
    # smRNA-Seq
    'ENCSR828LSC',
    'ENCSR000CRJ',
    'ENCSR911GQI',
    'ENCSR844HLP',
    'ENCSR291IZK',
    # shotgun sequencing?
    'ENCSR547BAN',
    'ENCSR601RDJ',
    'ENCSR832UXG'
]

res = []
meta = pd.DataFrame(columns=['lab', 'title', 'assay', 'description', 'num_files'])

for eid in expt_ids:
    resp = requests.get(expt_url.format(eid=eid))
    the_dat = [t for t in resp.content.splitlines() if len(t) and t[0] == '{']
    if len(the_dat) != 1:
        print "Error: experiment %s, found %d matches for the data block." % (eid, len(the_dat))
        continue
    else:
        the_dat = json.loads(the_dat[0])

    the_files = [t for t in the_dat['files'] if t['file_type'] == 'fastq']
    print "Expt %s, %d files" % (eid, len(the_files))

    meta.loc[eid] = pd.Series({
        'lab': the_dat['lab']['title'],
        'title': the_dat['biosample_summary'],
        'assay': the_dat['assay_title'],
        'description': the_dat['description'],
        'num_files': len(the_files)
    })

    matches = {}
    for f in the_files:
        ttl = f['title']
        lib = f['replicate']['library'].strip('/libraries/').strip('/')
        biol_repl = int(f['replicate']['biological_replicate_number'])
        tech_repl = int(f['replicate']['technical_replicate_number'])
        href = f['href']
        this_row = {
            'experiment': eid,
            'library': lib,
            'biol_repl': biol_repl,
            'tech_repl': tech_repl
        }
        if 'paired_end' in f:
            pe_n = int(f['paired_end'])
            pe_match = f['paired_with'].strip('/files/')
            if pe_match in matches:
                # we've already seen the matching read, so this one completes the row
                this_row = matches.pop(pe_match)
                this_row['read_%d' % pe_n] = href
                if pe_n == 1:
                    this_row['title'] = ttl
                res.append(this_row)
            else:
                this_row['read_%d' % pe_n] = href
                if pe_n == 1:
                    this_row['title'] = ttl
                matches[ttl] = this_row
        else:
            this_row['read_1'] = href
            this_row['title'] = ttl
            res.append(this_row)

# now go back through
# each library that has multiple listings with the same technical repl number and biol repl number is just an
# artificial sequencing run
run_counts = collections.Counter()
for t in res:
    k = (t['library'], t['biol_repl'], t['tech_repl'])
    run_counts[k] += 1
    t['seq_run'] = run_counts[k]
    t['outdir'] = "%s/%d" % (t['experiment'], t['seq_run'])

# finally, generate separate entries for each paired end entry
res_sep = []
for t in res:
    datum1 = dict(t)
    datum1['href'] = datum1.pop('read_1')
    if 'read_2' in t:
        datum1.pop('read_2')
        outfile1 = "%s_biol%d_tech%d_1.fastq.gz" % (
            datum1['title'],
            datum1['biol_repl'],
            datum1['tech_repl'],
        )
        datum1['outfile'] = outfile1

        datum2 = dict(t)
        datum2['href'] = datum2.pop('read_2')
        datum2.pop('read_1')
        outfile2 = "%s_biol%d_tech%d_2.fastq.gz" % (
            datum2['title'],
            datum2['biol_repl'],
            datum2['tech_repl'],
        )
        datum2['outfile'] = outfile2
        res_sep.append(datum2)
    else:
        datum1['outfile'] = "%s_biol%d_tech%d.fastq.gz" % (
            datum1['title'],
            datum1['biol_repl'],
            datum1['tech_repl'],
        )

    res_sep.append(datum1)


# now create a parameters file for Apocrita
outdir = unique_output_dir("encode_link_list", reuse_empty=True)
with open(os.path.join(outdir, "encode_roadmap.params"), 'wb') as f:
    c = csv.writer(f)
    for r in res_sep:
        c.writerow([
            r['experiment'],
            base_url + r['href'].strip('/'),
            r['outdir'],
            r['outfile'],
            ''
        ])  # that last comma ensures no random characters on the end of the filename

# write meta file
meta.to_csv(os.path.join(outdir, 'sources.csv'))