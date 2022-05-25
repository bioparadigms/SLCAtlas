#!/usr/bin/env python
# vim: ai

import glob
import sys
import gzip
import logging

import sql
import sql.aggregate

sys.path.append('./pylib')

from historytable import Database, HistoryTable, Table

user = 'upload-hmm-hits.py'

#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

DB = Database()

hmm_hits = Table('hmm_hits', DB)

# we upload HMM hits and mark them as selected according to criteria
# (bit_score)

hits = dict()
pfam_names = dict()
for fn in glob.glob("output/*.domtblout.txt"):
    for line in open(fn):
        if line.startswith("#"): continue
        l = line.split(None, 22)
        uniprot_acc = l[0].split('|')[1]
        gene_description = l[22].strip()
        gene_name = None
        for s in gene_description.split():
            if s.startswith("GN="):
                gene_name = s.split("=", 1)[1]
                break
        pfam_id = l[4].split('.')[0]
        pfam_name = l[3]
        if pfam_id == "-": pfam_id = pfam_name
        if pfam_id not in pfam_names: pfam_names[pfam_id] = pfam_name
        hmm_start, hmm_end = int(l[15]), int(l[16])
        start, end = int(l[19]), int(l[20])
        hmm_coverage = float(int(l[16])-int(l[15]))/float(l[5])
        coverage = float(int(l[20])-int(l[19]))/float(l[2])
        hit_evalue = float(l[6])
        # the i-Evalue from domtblout
        dom_ievalue = float(l[12])
        bit_score = float(l[7])
        key = (pfam_id, uniprot_acc)
        #if (hit_evalue > 0.001) or (dom_evalue > 0.001) or (hmm_coverage < 0.40):
        #    #if (hit_score < 20.0):
        #    logging.info("skipping hit because of low score: {}".format(line))
        #    continue
        if key not in hits:
            logging.info('found hit for gene {}: {}'.format(uniprot_acc,
                                                            line.strip()))
            hits[key] = (gene_name, hmm_coverage, coverage, hit_evalue, dom_ievalue, bit_score,
                         hmm_start, hmm_end, start, end)
        else:
            if (bit_score > hits[key][5]):
                logging.info('updating hit for {} due to higher score: '
                             '{}'.format(uniprot_acc, line.strip()))
                hits[key] = (gene_name, hmm_coverage, coverage, hit_evalue, dom_ievalue, bit_score,
                             hmm_start, hmm_end, start, end)
            elif (bit_score == hits[key][5]) and \
                 ((hmm_coverage > hits[key][1]) or (coverage > hits[key][2])) and \
                 (dom_ievalue < hits[key][4]):
                #sys.stderr.write('updating hit for {} due to higher coverage: {}'\
                #                 .format(gene_name, line))
                hits[key] = (gene_name, hmm_coverage, coverage, hit_evalue, dom_ievalue, bit_score, 
                             hmm_start, hmm_end, start, end)

# make a dict:
#   uniprot1 => (pfam1, pfam2, pfam3...), uniprot2 => (pfam2, pfam4, pfam5, ..)
hits_by_gene = dict( (y, [x[0] for x in hits.keys() if x[1] == y]) for y in
                    [z[1] for z in hits.keys()] )

comments = []
comments.append('Based on HMM search on set of UniProt sequences, unfiltered, '
                'and marked as selected for bit score > 50.')

#q = hmm_hits.select(sql.aggregate.Max(hmm_hits.record_id))
#cur = DB.query(q)
#record_id = cur.fetchone()[0]
#if record_id is None: record_id = 1
#else: record_id += 1
#cur.close()

print('TRUNCATE TABLE {!s};'.format(hmm_hits))
for key, hit in hits.items():
    pfam_id, uniprot_acc = key
    pfam_name = pfam_names[pfam_id]
    status = 0
    if hit[5] > 50.0:
        # select hits with bit score > 50
        status = 1
    hmm_hits.print_insert_sql(
        {'accession': uniprot_acc,
         'pfam_id': pfam_id,
         'pfam_name': pfam_name,
         'hit_evalue': hit[3],
         'dom_ievalue': hit[4],
         'bit_score': hit[5],
         'hmm_coverage': hit[1],
         'hmm_start': int(hit[6]),
         'hmm_end': int(hit[7]),
         'seq_coverage': hit[2],
         'seq_start': int(hit[8]),
         'seq_end': int(hit[9]),
         'status': status},
        comments=comments)

    """
    print("INSERT INTO `hmm_hits` ("
          "`record_id`, "
          "`accession`, `pfam_id`, `pfam_name`, "
          "`hit_evalue`, `dom_ievalue`, `bit_score`, "
          "`hmm_coverage`, `hmm_start`, `hmm_end`, "
          "`seq_coverage`, `seq_start`, `seq_end`, "
          "`status`, "
          "`user`, `comments`) "
          "VALUES ({:d}, "
          "{}, {}, {}, "
          "{}, {}, {}, "
          "{}, {}, {}, "
          "{}, {}, {}, "
          "{}, "
          "{}, {} );".format(
              record_id,
              *[Database.format(x) for x in (
                  uniprot_acc, pfam_id, pfam_name,
                  hit[3], hit[4], hit[5],
                  hit[1], int(hit[6]), int(hit[7]),
                  hit[2], int(hit[8]), int(hit[9]),
                  status )],
              Database.format(user), Database.format('\n'.join(comments))))
    record_id += 1
    """
