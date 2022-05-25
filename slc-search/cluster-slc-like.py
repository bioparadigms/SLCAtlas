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
import numpy as np
import scipy
import scipy.cluster.hierarchy as hac

from historytable import Database, HistoryTable, Table

user = 'cluster-slc-like.py'

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

DB = Database()

hmm_hits = sql.Table('hmm_hits')
# we will use print_insert_sql, so we use Table, not sql.Table
slc_like = HistoryTable('slc_like', DB)
slc_like_cluster_constraints = HistoryTable('slc_like_cluster_constraints', DB)
slc_like_clusters_linkage = Table('slc_like_clusters_linkage', DB)
slc_like_clusters = Table('slc_like_clusters', DB)
slc_families_defining_members = \
        HistoryTable('slc_families_defining_members', DB)
uniprot_proteins = HistoryTable('uniprot_proteins', DB)

q = slc_like.select(
    slc_like.accession,
    distinct=True,
    where=(slc_like.status >= 0),
    order_by=slc_like.record_id.asc)
cur = DB.query(q)
all_accessions = [x[0] for x in cur.fetchall()]
cur.close()
logging.info('read {:d} proteins'.format(len(all_accessions)))

# read uniprot gene symbols, we need these later
gene_symbols = dict()
for acc in all_accessions:
    q = uniprot_proteins.select(
        uniprot_proteins.symbol,
        where=(uniprot_proteins.accession == acc))
    cur = DB.query(q)
    rows = cur.fetchall()
    cur.close()
    if (len(rows) > 0) and (rows[0][0] is not None):
        gene_symbols[acc] = rows[0][0]
logging.info('found gene symbols for {:d} '
             'accessions'.format(len(gene_symbols)))

# read and filter hits, calculate score
# (acc, pfam_id) is unique in hmm_hits (see upload-hmm-hits.py)
slc_hits = dict()
all_pfams = set()
for acc in all_accessions:
    q = hmm_hits.select(
        hmm_hits.pfam_id,
        hmm_hits.bit_score,
        where=((hmm_hits.accession == acc) & (hmm_hits.bit_score > 25)))
    cur = DB.query(q, dictionary=True)
    for hit in cur.fetchall():
        slc_hits[(acc, hit['pfam_id'])] = hit['bit_score'] - 20.0
        all_pfams.add(hit['pfam_id'])
    cur.close()
all_pfams = sorted(all_pfams)
logging.info('read HMM hits of {:d} HMMs'.format(len(all_pfams)))

# turn hits into matrix
hits_matrix = np.zeros((len(all_accessions), len(all_pfams)))
for i, acc in enumerate(all_accessions):
    for j, pfam_id in enumerate(all_pfams):
        hits_matrix[i, j] = slc_hits.get( (acc, pfam_id), 0.0 )

Z = hac.linkage(hits_matrix, method='average', metric='cosine')
logging.info('length of linkage (Z) matrix: {:d}'.format(len(Z)))

##### WRITE Z-matrix here
print('TRUNCATE TABLE {!s};'.format(slc_like_clusters_linkage))
for i in range(len(Z)):
    slc_like_clusters_linkage.print_insert_sql(
        {'branch_id': i+len(all_accessions),
         'cluster_id1': int(Z[i, 0]),
         'cluster_id2': int(Z[i, 1]),
         'distance': Z[i, 2],
         'n_elements': int(Z[i, 3])})
#####

# try: cluster only to distance 0.7, find flat clusters
# here we collect all branches
# branches = {<branch_id>: [<acc1>, <acc2>]}
branches = dict()
# here only the top level ones
top_branches = dict()
for i, acc in enumerate(all_accessions):
    branches[i] = [acc]
    top_branches[i] = branches[i]
#print(branches)
i = 0
# do a complete clustering, but keep top branches that are below threshold
i_below_threshold = 0
while i < len(Z):
    k, l = int(Z[i, 0]), int(Z[i, 1])
    #print('joining {:d} ({!r}) and {:d} ({!r})'.format(k, branches[k], l, branches[l]))
    new_cl = len(all_accessions)+i
    branches[new_cl] = branches[k] + branches[l]
    try:
        assert(len(branches[new_cl]) == int(Z[i, 3]))
    except AssertionError as e:
        logging.error('branches[new_cl] = {!r}'.format(branches[new_cl]))
        logging.error(str(int(Z[i, 3])))
        raise e
    if Z[i, 2] < 0.7:
        i_below_threshold = i
        top_branches[new_cl] = branches[new_cl]
        del top_branches[k]
        del top_branches[l]
    i += 1
logging.info('used first {:d} rows of Z'.format(i_below_threshold))

##### APPLY CLUSTER CONSTRAINTS HERE

# read cluster constraints
q = slc_like_cluster_constraints.select(
    slc_like_cluster_constraints.accession1,
    slc_like_cluster_constraints.accession2,
    slc_like_cluster_constraints.kind,
    slc_like_cluster_constraints.comments)
cur = DB.query(q)
split_cons = []
join_cons = []
# (acc1, acc2, kind) => comments
cons_comments = dict()
for row in cur.fetchall():
    if row[2] == 1: join_cons.append((row[0], row[1]))
    elif row[2] == 2: split_cons.append((row[0], row[1]))
    cons_comments[(row[0], row[1], row[2])] = row[3]
cur.close()
logging.info('read {:d} join and {:d} split constraints'.format(
    len(join_cons), len(split_cons)))

acc2branch_id = dict([ (z, x) for x, y in top_branches.items() for z in y ])

queue = list(top_branches.keys())
# check split violations
used_split_cons = set()
for a1, a2 in split_cons:
    while True:
        br1 = acc2branch_id[a1]
        br2 = acc2branch_id[a2]
        if br1 != br2: break
        # we have a split violation
        branch_id = br1
        # this should not happen
        if branch_id < len(all_accessions): break
        # split cluster to new br1, br2
        br1, br2 = (int(Z[branch_id-len(all_accessions), 0]),
                    int(Z[branch_id-len(all_accessions), 1]))
        logging.info('splitting branch #{:d} into #{:d} and #{:d} due to '
                     'violation of constraint between {} and {} ("{}")'.format(
                         branch_id, br1, br2, a1, a2, 
                         cons_comments[(a1, a2, 2)]))
        used_split_cons.add((a1, a2))
        top_branches[br1] = branches[br1]
        top_branches[br2] = branches[br2]
        del top_branches[branch_id]
        for acc in branches[br1]:
            acc2branch_id[acc] = br1
        for acc in branches[br2]:
            acc2branch_id[acc] = br2

# check for unused split constraints
for a1, a2 in split_cons:
    if (a1, a2) not in used_split_cons:
        logging.warning('unused split constraint between {} and {} ("{}")'.format(
            a1, a2, cons_comments[(a1, a2, 2)]))

# check join constraints
used_join_cons = set()
for a1, a2 in join_cons:
    br1 = acc2branch_id[a1]
    br2 = acc2branch_id[a2]
    if br1 == br2:
        logging.debug('constraint satisfied, both are in branch #{:d}: {} {} '
                      '("{}")'.format(br1, a1, a2, cons_comments[(a1, a2, 1)]))
        continue
    logging.debug('constraint violated: #{:d} and #{:d}, {} {} ("{}")'\
                  .format(br1, br2, a1, a2, cons_comments[(a1, a2, 1)]))
    # we have a join violation
    # find lowest common parent
    b1 = br1
    b2 = br2
    i = min(b1, b2)-len(all_accessions)+1
    while (b1 != b2) and (i < len(Z)):
        Zi0 = int(Z[i, 0])
        Zi1 = int(Z[i, 1])
        if (b1 == Zi0) or (b1 == Zi1):
            b1 = i+len(all_accessions)
            logging.debug('merging #{:d} and #{:d} into #{:d} (b1)'\
                          .format(Zi0, Zi1, b1))
        if (b2 == Zi0) or (b2 == Zi1):
            b2 = i+len(all_accessions)
            logging.debug('merging #{:d} and #{:d} into #{:d} (b2)'\
                          .format(Zi0, Zi1, b2))
        i += 1
    branch_id = b1
    logging.info('joining branches #{:d} and #{:d} to lowest common '
                 'parent #{:d} due to join constraint between {} and '
                 '{} ("{}")'.format(br1, br2, branch_id, a1, a2,
                                   cons_comments[(a1, a2, 1)]))
    used_join_cons.add((a1, a2))
    logging.debug('adding top branch #{:d}'.format(branch_id))
    top_branches[branch_id] = branches[branch_id]
    for acc in top_branches[branch_id]:
        br = acc2branch_id[acc]
        if (br != branch_id) and (br in top_branches):
            logging.debug('removing top branch #{:d}'.format(br))
            del top_branches[br]
        acc2branch_id[acc] = branch_id

for a1, a2 in join_cons:
    if (a1, a2) not in used_join_cons:
        logging.warning('unused join constraint between {} and {} ("{}")'.format(
            a1, a2, cons_comments[(a1, a2, 1)]))

#####


##### BELOW HERE, CLUSTERS ARE NOT CHANGED, JUST FAMILY LABELING IS DONE

# read defining members
q = slc_families_defining_members.select()
cur = DB.query(q, dictionary=True)
defining_members = dict()
for row in cur.fetchall():
    defining_members[row['accession']] = row['family_name']
cur.close()

logging.info('read {:d} defining family members'.format(len(defining_members)))

rx_numbers = re.compile(r'(\d+)')

# get flat clusters
# also do labeling and check labeling conflicts
clusters = []
family2branch_id = dict()
for branch_id, members in top_branches.items():
    # from gene symbols
    family_names1 = set()
    # from defining members
    family_names2 = set()
    for acc in members:
        symbol = gene_symbols.get(acc, acc)
        # let's take into account non-human names as well
        # to pinpoint nomenclature problems
        if symbol.startswith('SLCO') or symbol.startswith('Slco'):
            family_names1.add('SLCO')
        elif symbol.startswith('SLC') or symbol.startswith('Slc'):
            family_names1.add(
                ''.join(rx_numbers.split(symbol)[:2]).upper() )
        if acc in defining_members:
            family_names2.add(defining_members[acc])
    # default cluster/family name
    family_names1 = sorted(family_names1)
    family_name = 'pSLC.{:d}'.format(branch_id)
    if len(family_names1) > 0:
        family_name = family_names1[0]
    if len(family_names1) > 1:
        logging.warning('conflicting family names based on gene symbols: '
                        '{} in branch #{:d}'\
                        .format(' '.join(family_names1), branch_id))
    # defining members take precedence
    family_names2 = sorted(family_names2)
    if len(family_names2) > 0:
        family_name = family_names2[0]
    if len(family_names2) > 1:
        logging.warning('multiple defining members with conflicting family '
                        'names: {} in branch #{:d}'\
                        .format(' '.join(family_names2), branch_id))
    if family_name in family2branch_id:
        family2branch_id[family_name].append(branch_id)
    else:
        family2branch_id[family_name] = [branch_id]
    clusters.append((family_name, branch_id, members))

for family_name, branch_ids in family2branch_id.items():
    if len(branch_ids) > 1:
        logging.warning('multiple branches for family {}: {}'.format(
            family_name, ' '.join(['#{:d}'.format(x) for x in branch_ids])))

#clusters.sort(key=lambda a: len(a[1]), reverse=True)
# sort by branch_id
clusters.sort(key=lambda a: a[1])
print('TRUNCATE TABLE {!s};'.format(slc_like_clusters))
for name, branch_id, members in clusters:
    member_names = ['/'.join([gene_symbols[x], x]) if x in gene_symbols else x
                    for x in members]
    logging.info('#{:d}: {}: {!r}'.format(branch_id, name, member_names))
    for acc in members:
        slc_like_clusters.print_insert_sql(
            {'accession': acc,
             'branch_id': branch_id,
             'family_name': name})

