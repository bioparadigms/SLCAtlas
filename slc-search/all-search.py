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
from Bio import SeqIO
from io import StringIO

from historytable import Database, HistoryTable

user = 'all-search.py'

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

DB = Database()

tcdb_families = HistoryTable('tcdb_families', DB)
tcdb_members = HistoryTable('tcdb_members', DB)
tcdb_sequences = HistoryTable('tcdb_sequences', DB)
pfam_domains = HistoryTable('pfam_domains', DB)

# active families #1
q = tcdb_families.select(
    tcdb_families.family_id3,
    where=(tcdb_families.status > 0))
cur = DB.query(q)
active_tcdb_families = [x[0] for x in cur.fetchall()]
cur.close()

q = tcdb_sequences.select(
    tcdb_sequences.tcdb_id,
    distinct=True,
    where=(tcdb_sequences.status > 0))
cur = DB.query(q)
active_tcdb_sequences = [x[0] for x in cur.fetchall()]
cur.close()

active_tcdb_families2 = ['.'.join(x.split('.')[:3]) for x in
                         active_tcdb_sequences]

if active_tcdb_families != active_tcdb_families2:
    d1 = list(set(active_tcdb_families)-set(active_tcdb_families2))
    if d1:
        logging.warning('deselecting family because they do not contain any '
                        'selected sequences: {!r}'.format(d1))
        for family_id3 in d1:
            active_tcdb_families.remove(family_id3)

    d2 = list(set(active_tcdb_families2)-set(active_tcdb_families))
    if d2:
        logging.warning('selecting family because it has selected sequences: '
                        '{!r}'.format(d2))
        for family_id3 in d2:
            active_tcdb_families.append(family_id3)

active_tcdb_subfamilies = ['.'.join(x.split('.')[:4]) for x in
                           active_tcdb_sequences]

q = pfam_domains.select(
    pfam_domains.dom_name,
    distinct=True,
    where=(pfam_domains.status > 0))
cur = DB.query(q)
active_pfam_domains = [x[0] for x in cur.fetchall()]
cur.close()

# check if all files exist
for family_id3 in active_tcdb_families:
    normal_fns = list(glob.glob('hmms/{}.hmm'.format(family_id3)))
    split_fns = list(glob.glob('hmms/{}-*.hmm'.format(family_id3)))
    if not (normal_fns + split_fns):
        logging.warning('cannot find HMM for family {}'.format(family_id3))
    else:
        # check if sequences are the same
        file_seqs = []
        # pool split parts together
        for fasta_fn in (
            list(glob.glob('hmm-fastas/{}.fasta'.format(family_id3))) +
            list(glob.glob('hmm-fastas/{}-*.fasta'.format(family_id3))) ):
            file_seqs.extend(list(SeqIO.parse(fasta_fn, 'fasta')))
        # get corresponding sequences
        q = tcdb_sequences.select(
            tcdb_sequences.record_id,
            tcdb_sequences.tcdb_id,
            tcdb_sequences.seq_fasta,
            where=((tcdb_sequences.tcdb_id.like('{}.%'.format(family_id3))) &
                   (tcdb_sequences.status > 0)))
        cur = DB.query(q)
        db_seqs = [(x[0], x[1], SeqIO.read(StringIO(x[2]), 'fasta'))
                   for x in cur.fetchall()]
        cur.close()
        file_seqs_set = frozenset([str(x.seq) for x in file_seqs])
        db_seqs_set = frozenset([str(x[2].seq) for x in db_seqs])
        for seq in file_seqs:
            if str(seq.seq) not in db_seqs_set:
                logging.warning('unknown sequence in HMM: {}'.format(seq.id))
        for record_id, tcdb_id, seq in db_seqs:
            if str(seq.seq) not in file_seqs_set:
                logging.warning('sequence not in HMM: {} record_id={:d} '
                                'tcdb_id={}'.format(
                                    seq.id, record_id, tcdb_id))
        logging.info('{}: file={:d}, db={:d}'.format(family_id3,
                                                     len(file_seqs),
                                                     len(db_seqs)))

# there are duplicates somehow?
active_tcdb_subfamilies = sorted(set(active_tcdb_subfamilies))
for family_id4 in active_tcdb_subfamilies:
    #if not (list(glob.glob('hmms/{}.hmm'.format(family_id4))) +
    #        list(glob.glob('hmms/{}-*.hmm'.format(family_id4)))):
    #    logging.warning('cannot find HMM for subfamily {}'.format(family_id4))
    normal_fns = list(glob.glob('hmms/{}.hmm'.format(family_id4)))
    split_fns = list(glob.glob('hmms/{}-*.hmm'.format(family_id4)))
    if not (normal_fns + split_fns):
        logging.warning('cannot find HMM for subfamily {}'.format(family_id4))
    else:
        # check if sequences are the same
        file_seqs = []
        # pool split parts together
        for fasta_fn in (
            list(glob.glob('hmm-fastas/{}.fasta'.format(family_id4))) +
            list(glob.glob('hmm-fastas/{}-*.fasta'.format(family_id4))) ):
            file_seqs.extend(list(SeqIO.parse(fasta_fn, 'fasta')))
        # get corresponding sequences
        q = tcdb_sequences.select(
            tcdb_sequences.record_id,
            tcdb_sequences.tcdb_id,
            tcdb_sequences.seq_fasta,
            where=((tcdb_sequences.tcdb_id.like('{}.%'.format(family_id4))) &
                   (tcdb_sequences.status > 0)))
        cur = DB.query(q)
        db_seqs = [(x[0], x[1], SeqIO.read(StringIO(x[2]), 'fasta'))
                   for x in cur.fetchall()]
        cur.close()
        file_seqs_set = frozenset([str(x.seq) for x in file_seqs])
        db_seqs_set = frozenset([str(x[2].seq) for x in db_seqs])
        for seq in file_seqs:
            if str(seq.seq) not in db_seqs_set:
                logging.warning('unknown sequence in HMM: {}'.format(seq.id))
        for record_id, tcdb_id, seq in db_seqs:
            if str(seq.seq) not in file_seqs_set:
                logging.warning('sequence not in HMM: {} record_id={:d} '
                                'tcdb_id={}'.format(
                                    seq.id, record_id, tcdb_id))
        logging.info('{}: file={:d}, db={:d}'.format(family_id4,
                                                     len(file_seqs),
                                                     len(db_seqs)))

for pfam_name in active_pfam_domains:
    if not (list(glob.glob('hmms/{}.hmm'.format(pfam_name)))):
        logging.warning('cannot find HMM for Pfam domain {}'.format(pfam_name))

active_hmms = []
for hmm_fn in glob.glob('hmms/*.hmm'):
    hmm_name = os.path.splitext(os.path.basename(hmm_fn))[0]
    if '.' in hmm_name:
        # tcdb family
        tcdb_name = hmm_name.split('-')[0]
        if ((tcdb_name not in active_tcdb_families) and 
            (tcdb_name not in active_tcdb_subfamilies)):
            logging.warning('unknown TCDB HMM found: {}'.format(hmm_fn))
            continue
    else:
        if hmm_name not in active_pfam_domains:
            logging.warning('unknown Pfam HMM found: {}'.format(hmm_fn))
            continue
    active_hmms.append((hmm_name, hmm_fn))

# do the actual search
if not os.path.isdir('output'): os.mkdir('output')

active_hmms.sort(key=lambda a: a[0])

for hmm_name, hmm_fn in active_hmms:
    print(hmm_name)
    subprocess.call('hmmsearch '
                    '--domtblout output/{hmm_name}.domtblout.txt '
                    '--tblout output/{hmm_name}.tblout.txt '
                    '-o output/{hmm_name}.out '
                    '{hmm_fn} uniprot-all.fasta'.format(
                        hmm_name=hmm_name, hmm_fn=hmm_fn), shell=True)

