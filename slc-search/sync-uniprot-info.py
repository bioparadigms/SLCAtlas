#!/usr/bin/env python

import sql
import sql.aggregate
import logging
import sys

sys.path.append('./pylib')

import glob
import gzip
import collections
import itertools
import os.path

from historytable import Database, HistoryTable

user = 'sync-uniprot-info.py'

# SYNC infos stored in:
#   uniprot_dbrefs / history / records
#   uniprot_proteins
# BASED ON hmm_hits

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

DB = Database()

hmm_hits = sql.Table('hmm_hits')

slc_like = HistoryTable('slc_like', DB)

uniprot_proteins = HistoryTable('uniprot_proteins', DB)
uniprot_dbrefs = HistoryTable('uniprot_dbrefs', DB)

# get selected hits
q = hmm_hits.select(hmm_hits.accession,
                    distinct=True,
                    where=(hmm_hits.status > 0))
cur = DB.query(q)
all_accessions = set([x[0] for x in cur.fetchall()])
cur.close()

# get slc-like proteins (could have manually added proteins!)
q = slc_like.select(slc_like.accession,
                    distinct=True,
                    where=(slc_like.status > 0))
cur = DB.query(q)
all_slc_like_accessions = set([x[0] for x in cur.fetchall()])
cur.close()

uniprot_fns = sys.argv[1:]

comments1 = []
comments1.append('Updated based on hmm_hits with status > 0.')
comments1.append('Files used:')
for fn in uniprot_fns:
    comments1.append(fn)

comments2 = comments1[:]
comments2[0] = 'Updated based on slc_like with status > 0.'

def fasta_records(f):
    record = None
    for line in f:
        if not line.strip():
            if record: yield record
            record = []
            for line in f:
                if line.startswith(">"): break
        if line.startswith(">"):
            if record: yield record
            record = []
        record.append(line)
    if record: yield record

fasta_store = dict()

# set everything to accepted
for uniprot_fn in uniprot_fns:
    if uniprot_fn.endswith('.gz'):
        f = gzip.open(uniprot_fn, 'rt', encoding='utf-8')
    else:
        f = open(uniprot_fn, 'rt', encoding='utf-8')

    if uniprot_fn.endswith('.fasta') or uniprot_fn.endswith('.fasta.gz'):
        for fasta_record in fasta_records(f):
            accession = fasta_record[0].split()[0].split('|')[1]
            if ((accession in all_accessions) and 
                (accession not in fasta_store)):
                fasta_store[accession] = fasta_record
    elif uniprot_fn.endswith('.txt') or uniprot_fn.endswith('.txt.gz'):
        uniprot_id = None
        reviewed = None
        name = None
        symbol = None
        tax_id = None
        fragment = False
        accessions = set()
        annots = []
        for line in f:
            if line.startswith('ID   '):
                uniprot_id = line.strip().split()[1].rstrip(';')
                reviewed = (line.strip().split()[2] == 'Reviewed;')
            elif line.startswith('AC   '):
                for acc in line.strip().split()[1:]:
                    accessions.add(acc.rstrip(';'))
            elif ((line.startswith('DE   RecName: Full=') or
                   line.startswith('DE   SubName: Full=')) and (name is None)):
                name = line.strip().split('=')[1].rstrip(';')\
                        .split('{')[0].strip()
            elif line.startswith('DE   Flags: Fragment;'):
                fragment = True
            elif line.startswith('GN   Name='):
                symbol = line.strip().split()[1].split('=')[1].rstrip(';')
            elif line.startswith('OX   NCBI_TaxID='):
                tax_id = int(line.strip().split()[1].split('=')[1].rstrip(';'))
            elif line.startswith("DR   HGNC;"):
                annots.append( ("HGNC", line.strip().split(";")[1].strip()) )
            elif line.startswith("DR   GeneID;"):
                annots.append( ("GeneID", line.strip().split(";")[1].strip()) )
            elif line.startswith("DR   UniGene;"):
                annots.append( ("UniGene", line.strip().split(";")[1].strip()) )
            elif line.startswith("DR   FlyBase;"):
                annots.append( ("FlyBase", line.strip().split(";")[1].strip()) )
            elif line.startswith("DR   KEGG;"):
                annots.append( ("KEGG", line.strip().split(";")[1].strip()) )
            elif line.startswith('//'):
                if fragment and (name is not None): name = name + ' (Fragment)'
                for acc in accessions & (all_accessions |
                                         all_slc_like_accessions):
                    # check GN, OX and ID fields...
                    q = uniprot_proteins.select(
                        where=(uniprot_proteins.accession == acc))
                    cur = DB.query(q, dictionary=True)
                    rows = cur.fetchall()
                    cur.close()

                    if acc not in fasta_store:
                        logging.warning('accession {} has no associated fasta '
                                        'sequence, skipping'.format(acc))
                        continue

                    new_row = {
                        'tax_id': tax_id,
                        'accession': acc,
                        'uniprot_id': uniprot_id,
                        'name': name,
                        'symbol': symbol,
                        'seq_fasta': ''.join(fasta_store[acc]),
                        'reviewed': int(reviewed)}

                    if acc in all_accessions: comments = comments1
                    if acc in all_slc_like_accessions: comments = comments2
                    if len(rows) == 0:
                        uniprot_proteins.insert_row(new_row, manual=0,
                                                    user=user,
                                                    comments=comments)
                    else:
                        if len(rows) > 1:
                            logging.warning('duplicate records found for '
                                            '{}'.format(acc))
                        for old_row in rows:
                            uniprot_proteins.update_row(
                                old_row, new_row, manual=0, user=user,
                                comments=comments)

                    # remove duplicates
                    annots = list(dict(zip(annots, itertools.repeat(1))).keys())

                    q = uniprot_dbrefs.select(
                        uniprot_dbrefs.db,
                        uniprot_dbrefs.xref,
                        where=(uniprot_dbrefs.accession == acc))
                    cur = DB.query(q)
                    old_annots = [tuple(x) for x in cur.fetchall()]
                    if old_annots != annots:
                        uniprot_dbrefs.update_one_to_many(('db', 'xref'),
                                                          {'accession': acc},
                                                          annots, manual=0,
                                                          user=user,
                                                          comments=comments)

                uniprot_id = None
                reviewed = None
                name = None
                symbol = None
                tax_id = None
                fragment = False
                accessions = set()
                annots = []
    f.close()

