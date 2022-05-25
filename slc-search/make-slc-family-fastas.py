#!/usr/bin/env python
# vim: ai

import sql
import sql.aggregate
import logging
import sys

sys.path.append('./pylib')

import glob
import os.path
import subprocess
import itertools

from historytable import Database, HistoryTable, Table

import hashlib

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

DB = Database()

slc_like = HistoryTable('slc_like', DB)
slc_like_clusters = Table('slc_like_clusters', DB)
uniprot_proteins = HistoryTable('uniprot_proteins', DB)

q = slc_like.join(
    uniprot_proteins, type_='LEFT',
    condition=(slc_like.accession == uniprot_proteins.accession)
).join(
    slc_like_clusters, type_='LEFT',
    condition=(slc_like.accession == slc_like_clusters.accession)
).select(
    slc_like.accession,
    slc_like_clusters.branch_id,
    slc_like_clusters.family_name,
    uniprot_proteins.tax_id,
    uniprot_proteins.seq_fasta,
    where=(slc_like.status >= 0)
)

cur = DB.query(q, dictionary=True)
all_accessions = cur.fetchall()
cur.close()

all_accessions.sort(key=lambda a: a['family_name'])

# branch_id IS NOT compatible with family_name
# one family can have multiple branches, sometimes it's not possible to merge
# all branches, e.g. SLC38 or SLC51 families.
for family_name, members in itertools.groupby(all_accessions,
                                              key=lambda a: a['family_name']):
    members = sorted(members, key=lambda a: a['accession'])
    # generate key unique to family members
    key = hashlib.sha1(
        ('\n'.join([x['accession'] for x in members])).encode('utf-8')
    ).hexdigest()
    logging.info('{}: ({}) {:d} members'.format(
        members[0]['family_name'],
        ', '.join(['#{:d}'.format(x) for x in set([y['branch_id'] for y in
                                                   members])]),
        len(members)
    ))
    fasta_fn = 'branches-fastas/branch_{}.fasta'.format(key)
    if not os.path.exists(fasta_fn):
        logging.info(f'writing {fasta_fn}')
        with open(fasta_fn, 'wt',
                  encoding='utf-8') as f:
            for m in members:
                f.write(m['seq_fasta'])

