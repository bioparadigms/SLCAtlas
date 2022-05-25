#!/bin/env/python
# vim: ai

import mysql.connector
import sql
import sql.aggregate
import logging
import sys

sys.path.append('../pylib')

from historytable import Database, HistoryTable, Table

logging.basicConfig(level=logging.DEBUG)

DB = Database()

tcdb_sequences = sql.Table('tcdb_sequences')

q = tcdb_sequences.select()
cur = DB.query(q, dictionary=True)
all_selected_sequences = cur.fetchall()
cur.close()

for seq in all_selected_sequences:
    if seq['seq_fasta'] is None:
        logging.warning('{} / {} -> no sequence'.format(seq['tcdb_id'], seq['accession']))
        continue
    else:
        logging.info('{} / {}'.format(seq['tcdb_id'], seq['accession']))
        print(seq['seq_fasta'])

