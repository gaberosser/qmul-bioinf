from bs4 import BeautifulSoup
import requests
import re
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt, patches, collections
from collections import OrderedDict


def get_titles():
    url = 'http://www.nature.com/nrg/series/nextgeneration/index.html?WT.ec_id=NRG-201507'
    resp = requests.get(url)
    soup = BeautifulSoup(resp.content, 'html.parser')
    qset = soup.find_all('h4', attrs = {'class': 'atl'})
    return [t.text.strip().replace(u'\u2013', '-').replace(u'\u2014', '-') for t in qset]


def mpl_circle(coords, *colours, **kwargs):
    if 'radius' in kwargs:
        radius = kwargs['radius']
    else:
        radius = 1

    if 'theta0' in kwargs:
        th0 = kwargs['theta0']
    else:
        th0 = 90
    n = len(colours)
    if n == 0:
        raise ValueError("Must specify at least one colour")
    dth = 360 / float(n)
    res = []
    for i in range(n):
        res.append(patches.Wedge(coords, radius, th0, th0 + dth, edgecolor='none', facecolor=colours[i]))
        th0 += dth
    return collections.PatchCollection(res, match_original=True)


if __name__ == '__main__':
    SHOW_SINGLETONS = False
    categories = [
        ('health', [
            r'clinic.*',
            r'disorder.*',
            r'cancer',
            r'disease',
            r'medicine',
        ]),
        ('cancer', [
            r'cancer',
            r'oncol.*'
        ]),
        ('population history', [
            r'histor.*'
        ]),
        ('epigenetics', [
            r'epigen.*',
            r'histone',
            r'methyl.*',
        ]),
        ('RNA-Seq', [
            r'RNA',
            r'rna.seq.*',
            r'non.coding',
            r'small rna',
            'transcriptom.*'
        ]),
        ('ChIP-Seq', [
            r'chip-seq',
        ]),
        ('prokaryotes', [
            r'prokaryot.*',
            r'bacteri.*',
            r'microb.*'
        ]),
        ('marine', [
            r'marine',
            r'aquatic'
        ]),
        ('ribosome', [r'ribosom.*']),
        ('evolution', [
            r'phylogen.*',
            r'evol.*'
        ]),
        ('protein', [
            r'protein',
            r'proteomi.*'
        ]),
        ('human', [r'human']),
        ('mouse', [r'mouse']),
        ('mammal', [r'human', r'mouse', r'primate', r'chimp.*', r'mammal.*']),
        ('regulation', [r'transcription factor', r'regulat.*']),
        ('variation', [r'polymorphi.*', r'muta.*', r'varia.*']),
        ('methods', [r'single.cell', r'integrat.*', r'replicate', r'reproducib.*', r'annotat.*', r'reference'])
    ]

    # titles = get_titles()
    with open('nature_ngs_titles.txt', 'rb') as f:
        titles = f.read().split('\n')

    ncat = len(categories)
    cat_regex = [(a, re.compile('|'.join(b), flags=re.I)) for a, b in categories]
    cmap = plt.cm.jet
    cat_colours = OrderedDict([
        (cat, cmap(np.linspace(0, 1, ncat)[i])) for i, (cat, _) in enumerate(categories)
    ])
    # make graph
    g = nx.Graph()
    node_colours = {}
    node_groups = {}
    in_graph = {}

    for t in titles:
        node_groups.setdefault(t, set())
        added = False
        for i, (cat, regex) in enumerate(cat_regex):
            if re.search(regex, t):
                added = True
                node_groups[t].add(cat)
                g.add_node(t, groups=node_groups[t])
                if len(in_graph.setdefault(cat, [])) > 0:
                    [g.add_edge(t, x) for x in in_graph[cat]]
                in_graph[cat].append(t)

        if not added and SHOW_SINGLETONS:
            g.add_node(t, groups=node_groups[t])

    for t in titles:
        if (t in g.nodes()) or SHOW_SINGLETONS:
            node_colours[t] = [cat_colours[x] for x in g.node[t]['groups']] or ['gray']

    coords = nx.spring_layout(g)  # key is node ID
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    lc = nx.draw_networkx_edges(g, pos=coords, ax=ax, alpha=0.4)
    kwds = {'radius': 0.02}
    for t, arr in node_colours.items():
        pc = mpl_circle(coords[t], *arr, **kwds)  # PatchCollection
        ax.add_collection(pc)

    ax.set_aspect('equal')
    ax.set_xlim([-.1, 1.5])
    ax.set_ylim([-.1, 1.1])
    ax.set_axis_off()

    # legend
    for cat, col in cat_colours.items():
        p = patches.Circle((999, 999), radius=0., label=cat, edgecolor='none', facecolor=col)
        ax.add_patch(p)

    ax.legend(loc='center right')
    fig.tight_layout()




