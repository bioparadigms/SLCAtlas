#!/usr/bin/env python

import mysql.connector
import sql
import sql.aggregate
import logging
import sys

sys.path.append('./pylib')

import glob
import re
import os.path
import networkx as nx

from historytable import Database, HistoryTable, Table

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

DB = Database()

hmm_hits = sql.Table('hmm_hits')
# we will use print_insert_sql, so we use Table, not sql.Table
hmm_hits_clusters = Table('hmm_hits_clusters', DB)
hmm_hits_cluster_constraints = sql.Table('hmm_hits_cluster_constraints')
uniprot_proteins = HistoryTable('uniprot_proteins', DB)
uniprot_dbrefs = HistoryTable('uniprot_dbrefs', DB)

q = hmm_hits.join(
    uniprot_proteins,
    type_='LEFT',
    condition=(hmm_hits.accession == uniprot_proteins.accession)
).select(
    hmm_hits.accession,
    uniprot_proteins.uniprot_id,
    uniprot_proteins.name,
    uniprot_proteins.symbol,
    uniprot_proteins.tax_id,
    uniprot_proteins.reviewed,
    distinct=True,
    where=(hmm_hits.status > 0))
cur = DB.query(q, dictionary=True)
proteins = cur.fetchall()
cur.close()

logging.info('proteins: {:d}'.format(len(proteins)))

clusters = dict()
uniprot_fingerprints = dict()
uniprot_tax_ids = dict()
sp_accessions = set()

# here we go through each sequence and create a fingerprint for it based on
# uniprot annotations
# we also collect all sequences that belong to a specific annotation
for protein in proteins:
    fp = dict()

    if protein['symbol'] is not None:
        fp['Gene Symbol'] = protein['symbol']
    if protein['reviewed'] == 1: sp_accessions.add(protein['accession'])
    q = uniprot_dbrefs.select(
        where=(uniprot_dbrefs.accession == protein['accession']))
    cur = DB.query(q, dictionary=True)
    dbrefs = cur.fetchall()
    cur.close()

    for dbref in dbrefs:
        if dbref['db'] not in fp:
            fp[dbref['db']] = dbref['xref']
        elif dbref['xref'] != fp[dbref['db']]:
            logging.warning('contradictions/double annotation for sequence: '
                            '{} {} old={} new={}'.format(
                                protein['accession'], dbref['db'],
                                fp[dbref['db']], dbref['xref']))
            if dbref['db'] == 'UniGene':
                logging.warning('ambivalent UniGene association, removing '
                                'UniGene fingerprint from sequence {} '
                                'completely.'.format(protein['accession']))
                fp[dbref['db']] = None
    del_keys = [key for key, value in fp.items() if value is None]
    for key in del_keys: del fp[key]

    uniprot_fingerprints[protein['accession']] = fp
    uniprot_tax_ids[protein['accession']] = protein['tax_id']
    for protein_id in fp.items():
        if (protein_id, protein['tax_id']) not in clusters:
            clusters[(protein_id, protein['tax_id'])] = []
        clusters[(protein_id, protein['tax_id'])].append(protein['accession'])

# we add each sequence as a node to the graph
G = nx.Graph()
for protein in proteins:
    G.add_node(protein['accession'])

# take clustering constraints into account
q = hmm_hits_cluster_constraints.select(
    hmm_hits_cluster_constraints.accession1,
    hmm_hits_cluster_constraints.accession2,
    hmm_hits_cluster_constraints.kind
)
cur = DB.query(q)
never_join = set()
# I want to keep track of those kind==2 constraints that have never been used
# at the end of clustering
never_join_hit = set()
# we use "coloring" to specify nodes connected to other nodes that are involved
# in a split constraint, to make the detection of the violation of split
# constraints easier
# in this system, one node can have several colors
#    acc -> set(2, 6, 4)
# this is because pairs clusters with non-forbidden color pairs can still be
# joined
# but still the color system helps to avoid pairwise lookups of cluster members
# upon a cluster join
node_colors = dict()
# then, we convert split constraints to forbidden color pairs
# this is a dict
#    (color1, color2) => (acc1, acc2), the original constraint that led to this
#    color constraint
forbidden_color_pairs = dict()
# color -> set([node1, node2, ...])
color2nodes = dict()
next_color = 1
for acc1, acc2, kind in cur.fetchall():
    if kind == 1:
        # join elements now
        logging.info('constraints: joining proteins '
                     '{} and {}'.format(acc1, acc2))
        G.add_edge(acc2, acc2)
    elif kind == 2:
        # never join these elements
        if acc2 < acc1: acc1, acc2 = acc2, acc1
        never_join.add((acc1, acc2))
        if acc1 in node_colors: color1 = node_colors[acc1]
        else:
            color1 = next_color
            next_color += 1
            node_colors[acc1] = frozenset([color1])
            color2nodes[color1] = set([acc1])
        if acc2 in node_colors: color2 = node_colors[acc2]
        else:
            color2 = next_color
            next_color += 1
            node_colors[acc2] = frozenset([color2])
            color2nodes[color2] = set([acc2])
        if color2 < color1: color1, color2 = color2, color1
        forbidden_color_pairs[(color1, color2)] = (acc1, acc2)
cur.close()

# here, we connect groups of sequences that share an annotation (xref).
for key, l in clusters.items():
    # split constraints are always going to lead to arbitraty cluster members
    # basically, so we try to introduce some consistency here
    # also this helps us in the loops below
    l.sort()
    logging.info('cluster "{!r}": {!r}'.format(key, l))
    for i in range(len(l)-1):
        dont_join = False
        if (l[i] in node_colors) and (l[i+1] in node_colors):
            # check if any color pairs are forbidden
            logging.info('check color clash for {!r} and {!r}'.format(
                node_colors[l[i]], node_colors[l[i+1]]))
            for c1 in node_colors[l[i]]:
                for c2 in node_colors[l[i+1]]:
                    if c2 < c1: c1, c2 = c2, c1
                    if (c1, c2) in forbidden_color_pairs:
                        # color clash, hit a split constraint
                        logging.info('not joining proteins {} and {} because '
                                     'of split constraint {!r} despite their '
                                     'common annotation '
                                     '{!r}'.format(l[i], l[i+1],
                                                  forbidden_color_pairs[(c1, c2)],
                                                  key))
                        never_join_hit.add(forbidden_color_pairs[(c1, c2)])
                        dont_join = True
        if dont_join: continue
        # calculate final color of all nodes as a "sum of all colors"
        color = (node_colors.get(l[i], frozenset([])) | 
                 node_colors.get(l[i+1], frozenset([])))
        if color:
            # every node of each color will get every color
            for c in color:
                for acc in color2nodes[c]:
                    logging.info('setting color of node {} to {!r}'.format(acc, color))
                    node_colors[acc] = color
                # this is important in case one of the nodes have been
                # uncolored
                node_colors[l[i]] = color
                node_colors[l[i+1]] = color
                color2nodes[c].add(l[i])
                color2nodes[c].add(l[i+1])
        G.add_edge(l[i], l[i+1])

# annotate clusters
# we collect common fingerprints of all members of each cluster
cl_fingerprints = dict()
u2cl = dict()
for c in nx.connected_components(G):
    logging.info('NEW CLUSTER')
    cl = tuple(sorted(c))
    fp = dict()
    for u in cl:
        logging.info('new sequence in cluster: {}, fp: {!r}'.format(u, uniprot_fingerprints[u]))
        for key, value in uniprot_fingerprints[u].items():
            if key in fp:
                if (fp[key] != value) and (key != 'Gene Symbol'):
                    logging.warning('contradicting annotations in cluster: '
                                    'current cluster fp: {!r}, offending acc, '
                                    'db and xref: {}, {}, {}'.format(
                                        fp, u, key, value))
            else:
                logging.info('found xref: {}: {}'.format(key, value))
                fp[key] = value
        u2cl[u] = cl
    cl_fingerprints[cl] = fp

# read BLAST results and connect sequences, if it doesn't contradict with gene
# annotation information.
sequence_lengths = dict()
for blast_out_fn in glob.glob('hits.*.blast.out'):
    with open(blast_out_fn, 'rt', encoding='utf-8') as f:
        for line in f:
            query_def, query_len, hit_def, hit_len, hsp_query_len, hsp_hit_len, identity = line.strip().split("\t")
            query_len = int(query_len)
            hit_len = int(hit_len)
            hsp_query_len = int(hsp_query_len)
            hsp_hit_len = int(hsp_hit_len)
            identity = int(identity)

            query_coverage = float(hsp_query_len)/float(query_len)
            hit_coverage = float(hsp_hit_len)/float(hit_len)
            query_identity_percent = float(identity)/float(query_len)
            hsp_identity = float(identity)/float(hsp_query_len)

            query_uniprot_acc = query_def.split("|")[1]
            hit_uniprot_acc = hit_def.split("|")[1]
            sequence_lengths[query_uniprot_acc] = query_len
            sequence_lengths[hit_uniprot_acc] = hit_len

            # this is now not necessary because we do BLAST searches
            # per-organism
            #### don't connect proteins from different organisms
            ###if (uniprot_tax_ids[query_uniprot_acc] !=
            ###    uniprot_tax_ids[hit_uniprot_acc]): continue

            if hsp_identity > 0.95:
                # connect query_def and hit_def
                if u2cl[query_uniprot_acc] == u2cl[hit_uniprot_acc]: continue
                fp_q = cl_fingerprints[u2cl[query_uniprot_acc]]
                fp_h = cl_fingerprints[u2cl[hit_uniprot_acc]]
                # check if compatible
                dont_join = False
                fp = dict()
                for key in sorted(set(list(fp_q.keys()) + list(fp_h.keys()))):
                    if key == "Gene Symbol": continue
                    if (key in fp_q) and (key in fp_h) and (fp_q[key] != fp_h[key]):
                        logging.warning('!!!!! contradicting join based on BLAST '
                                        'results, perhaps chimera sequence, or '
                                        'very similar proteins? BLAST query: {}, '
                                        'fp={!r}, hit: {}, fp={!r}'.format(
                                            query_def, fp_q, hit_def, fp_h))
                        #print(query_def)
                        #print(fp_q)
                        ##print(u2cl[query_uniprot_acc])
                        #print(hit_def)
                        #print(fp_h)
                        ##print(u2cl[hit_uniprot_acc])
                        dont_join = True
                        break
                    else:
                        fp[key] = fp_q.get(key, fp_h.get(key, None))
                # we bail out of color analysis if we anyway wouldn't join these nodes
                if dont_join: continue
                # check for split constraints
                if (query_uniprot_acc in node_colors) and (hit_uniprot_acc in node_colors):
                    # check if any color pairs are forbidden
                    for c1 in node_colors[query_uniprot_acc]:
                        for c2 in node_colors[hit_uniprot_acc]:
                            if c2 < c1: c1, c2 = c2, c1
                            if (c1, c2) in forbidden_color_pairs:
                                # color clash, hit a split constraint
                                logging.info('not joining proteins {} and {} '
                                             'because of split constraint {!r} '
                                             'despite high-scoring BLAST '
                                             'hit'.format(query_uniprot_acc,
                                                          hit_uniprot_acc,
                                                          forbidden_color_pairs[(c1, c2)])
                                             )
                                never_join_hit.add(forbidden_color_pairs[(c1, c2)])
                                dont_join = True
                # bail out if there is a color clash
                if dont_join: continue
                # calculate final color of all nodes as a "sum of all colors"
                color = (node_colors.get(query_uniprot_acc, frozenset([])) | 
                         node_colors.get(hit_uniprot_acc, frozenset([])))
                if color:
                    # every node of each color will get every color
                    for c in color:
                        for acc in color2nodes[c]: node_colors[acc] = color
                        node_colors[query_uniprot_acc] = color
                        node_colors[hit_uniprot_acc] = color
                        color2nodes[c].add(query_uniprot_acc)
                        color2nodes[c].add(hit_uniprot_acc)
                cl = tuple(sorted(list(u2cl[query_uniprot_acc]) + list(u2cl[hit_uniprot_acc])))
                cl_fingerprints[cl] = fp
                old_cl1 = u2cl[query_uniprot_acc]
                old_cl2 = u2cl[hit_uniprot_acc]
                del cl_fingerprints[old_cl1]
                del cl_fingerprints[old_cl2]
                del old_cl1
                del old_cl2
                for u in cl:
                    u2cl[u] = cl
                G.add_edge(query_uniprot_acc, hit_uniprot_acc)

proteins_by_acc = dict([(p['accession'], p) for p in proteins])

str_protein = lambda p: '{:<20s} ({}) {:<15s} {:>5d} "{}"'.format(
    p['accession'], p['reviewed'], 'NULL' if p['symbol'] is None else
    p['symbol'], uniprot_tax_ids[p['accession']], p['name'])

cluster_id = 1
print('TRUNCATE TABLE {!s};'.format(hmm_hits_clusters))

for c in nx.connected_components(G):
    reps = set()
    sp_ids = []
    for u in c:
        if u in sp_accessions: sp_ids.append(proteins_by_acc[u])
    if len(c) == 1:
        logging.info('***** SINGLE MEMBER CLUSTER')
        # then the only member is the representative
        reps.add(list(c)[0])
    if len(sp_ids) == 0:
        logging.info('***** NON-SWISSPROT CLUSTER')
        # we choose the longest sequence as representative
        cl_rep = max(sorted([(x, sequence_lengths.get(x, 0)) for x in c],
                            key=lambda a: a[1]), key=lambda a: a[1])
        logging.info('cluster representative: {}'.format(
            str_protein(proteins_by_acc[cl_rep[0]])))
        reps.add(cl_rep[0])
    elif len(sp_ids) == 1:
        logging.info('cluster representative: {}'.format(
            str_protein(sp_ids[0])))
        reps.add(sp_ids[0]['accession'])
    elif len(sp_ids) > 1:
        logging.warning('***** MULTI-SWISSPROT CLUSTER (this should not '
                        'normally happen)')
        for u in sp_ids:
            logging.warning(repr(u))
            reps.add(u['accession'])
    for u in c:
        logging.info(str_protein(proteins_by_acc[u]))
        hmm_hits_clusters.print_insert_sql(
            {'accession': u,
             'cluster_id': cluster_id,
             'representative': 1 if u in reps else 0})
    cluster_id += 1

for acc1, acc2 in never_join - never_join_hit:
    logging.info('split constraint never used: {} {}'.format(acc1, acc2))

