import requests
import json
import collections


def _run_dgidb_interaction_lookup(genes):
    url = 'http://dgidb.org/api/v2/interactions.json'

    resp = requests.get(url, {'genes': ','.join(genes)})
    dat = json.loads(resp.content)
    interactions = dict(
        [(t['searchTerm'], t['interactions']) for t in dat['matchedTerms']]
    )
    ambiguous = collections.defaultdict(list)
    for t in dat['ambiguousTerms']:
        ambiguous[t['searchTerm']].append(t)
    return {
        'interactions': interactions,
        'unmatched': dat['unmatchedTerms'][0].split(', ') if len(dat['unmatchedTerms']) > 0 else [],
        'ambiguous': ambiguous,
    }


def dgidb_lookup_drug_gene_interactions(genes):
    get_genes = ','.join(genes)
    if len(get_genes) > 2000:
        # chunk
        gg = list(genes)
        chunks = []
        while True:
            if len(gg) == 0:
                break
            this_ch = []
            while True:
                try:
                    this_g = gg.pop()
                except IndexError:
                    break

                this_ch.append(this_g)
                the_str = ','.join(this_ch)
                if len(the_str) > 1900:
                    break

            if len(this_ch) > 0:
                chunks.append(this_ch)

    else:
        chunks = [genes]

    res = {}
    for ch in chunks:
        this_res = _run_dgidb_interaction_lookup(ch)
        for k, v in this_res.items():
            if isinstance(v, dict):
                res.setdefault(k, {}).update(v)
            elif isinstance(v, list):
                res.setdefault(k, []).extend(v)
            else:
                raise TypeError("Unable to add chunked results")

    return res
