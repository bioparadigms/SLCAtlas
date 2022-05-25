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

from historytable import Database, HistoryTable

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

DB = Database()

hmm_hits = sql.Table('hmm_hits')
uniprot_proteins = HistoryTable('uniprot_proteins', DB)

q = hmm_hits.join(uniprot_proteins, type_='LEFT',
                  condition=(hmm_hits.accession == uniprot_proteins.accession)
                 ).select(
                     hmm_hits.accession,
                     uniprot_proteins.tax_id,
                     uniprot_proteins.seq_fasta,
                     distinct=True,
                     where=(hmm_hits.status > 0))
cur = DB.query(q)
all_accessions = cur.fetchall()
cur.close()

all_accessions.sort(key=lambda a: a[0])

n_tax = dict()
for acc, tax_id, seq_fasta in all_accessions:
    n_tax[tax_id] = n_tax.get(tax_id, 0) + 1

for tax_id, cnt in n_tax.items():
    logging.info('found {:d} sequences for tax_id={:d}'.format(cnt, tax_id))
    with open('hits.{:d}.fasta'.format(tax_id), 'wt', encoding='utf-8') as f:
        for acc, t_id, seq_fasta in all_accessions:
            #print(acc, t_id, tax_id, t_id == tax_id)
            if t_id == tax_id:
                print('writing {} to file {}, len {:d}'.format(acc, f.name,
                                                               len(seq_fasta)))
                f.write(seq_fasta)

