#!/usr/bin/env python
# vim: ai

import sql
import sql.aggregate
import logging
import sys

sys.path.append('../pylib')

import glob
import os.path
import subprocess
import itertools

from historytable import Database, HistoryTable, Table

from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
import dendropy as dpy
import io

import hashlib

logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(level=logging.INFO)

DB = Database()

hmm_hits = sql.Table('hmm_hits')
uniprot_proteins = HistoryTable('uniprot_proteins', DB)

##### START

fasta_fn = sys.argv[1]

aln_fn = fasta_fn.replace('.fasta', '.aln.fasta')
if not os.path.exists(aln_fn):
    subprocess.run(['clustalo', '-i', fasta_fn, '--iter=5', '-o', aln_fn])
# convert to phylip
msa = AlignIO.read(open(aln_fn, 'rt', encoding='utf-8'), format='fasta')
# rename
for i in range(len(msa)):
    msa[i].id = msa[i].id.split('|')[1]

# check for duplicate sequences and remove them
seqs = dict()
for seq in msa:
    seq_key = hashlib.sha256(str(seq.seq).encode('latin1')).hexdigest()
    logging.debug('read sequence {} with hash {}'.format(seq.id, seq_key))
    if seq_key in seqs:
        logging.warning('duplicate sequences: {} and '
                        '{}'.format(seqs[seq_key][0].id, seq.id))
        seqs[seq_key].append(seq)
    else:
        seqs[seq_key] = [seq]

final_seqs = []
for seq_key, seq_list in seqs.items():
    logging.debug('keeping sequence {} with hash {}'.format(seq_list[0].id, seq_key))
    final_seqs.append(seq_list[0])
msa_old = msa
msa = MultipleSeqAlignment(final_seqs)


phylip_fn = aln_fn.replace('.fasta', '.phylip')
if not os.path.exists(phylip_fn):
    logging.info('file {} does not exist, creating.'.format(phylip_fn))
    AlignIO.write(msa, open(phylip_fn, 'wt', encoding='utf-8'), format='phylip')

sms_dir = phylip_fn.replace('.phylip', '.sms')
if not os.path.exists('{}/{}_phyml_tree.txt'.format(sms_dir, phylip_fn)):
    logging.info('file {} does not exist, creating.'.format(
        '{}/{}_phyml_tree.txt'.format(sms_dir, phylip_fn)
    ))
    subprocess.run('/opt/sms-1.8.4/sms.sh -i {0} -d aa -o {1} -t -s SPR -r 10 -b -4 >{1}.out'.format(phylip_fn, sms_dir), shell=True)

# rename in tree
tax2tag = {
        9606: '_Homo',
        10090: '_Mus',
        10116: '_Rattus',
        7227: '_Drosophila',
        9031: '_Gallus',
        6239: '_Caenorhabditis',
        7955: '_Danio'
}
with open('{}/{}_phyml_tree.txt'.format(sms_dir, phylip_fn), 'rt',
          encoding='utf-8') as f:
    tree_str = f.read()

# get slc names
descriptions = [x for x in open(fasta_fn, 'rt', encoding='utf-8') if
                x.startswith('>')]
slc_names = dict()
for i in range(len(msa_old)):
    uid = msa_old[i].id
    q = uniprot_proteins.select(where=(uniprot_proteins.accession == uid))
    cur = DB.query(q, dictionary=True)
    p = cur.fetchone()
    cur.close()

    if p is not None:
        if p['symbol'] is not None: slc_name = p['symbol'] + tax2tag.get(p['tax_id'], '_other')
        else: slc_name = p['accession'] + tax2tag.get(p['tax_id'], '_other')
    else:
        # this is just in case we include custom proteins in our alignment,
        # e.g. as an outgroup
        tax_name = 'other'
        slc_name = uid
        for desc in descriptions:
            if uid in desc:
                line_parts = desc.split()
                for s in line_parts:
                    if s.startswith('GN='):
                        slc_name = s.split('=')[1]
                    elif s.startswith('OS='):
                        tax_name = s.split('=')[1]
                break
        slc_name = slc_name + '_' + tax_name
    # colons prevent dendropy from reading the tree normally
    slc_name = slc_name.replace(':', '-').replace('\\', '-')
    slc_names[uid] = slc_name+'_'+uid

    tree_str = tree_str.replace(uid, slc_name+'_'+uid)

# put back removed duplicate sequences
for seq_key, seq_list in seqs.items():
    if len(seq_list) > 1:
        slc_name1 = slc_names[seq_list[0].id]
        subtree_str = '(' + ','.join([
            '{}:0.0'.format(slc_names[x.id]) for x in seq_list
        ]) + ')1.0'
        tree_str = tree_str.replace(slc_name1, subtree_str)

renamed_tree_fn = '{}/{}_phyml_tree.renamed.txt'.format(sms_dir, phylip_fn)
if not os.path.exists(renamed_tree_fn):
    with open(renamed_tree_fn, 'wt',
              encoding='utf-8') as f:
        f.write(tree_str)

tree = dpy.Tree.get(path=renamed_tree_fn, schema='newick')
taxa = set([x.lstrip('_') for x in tax2tag.values()])
outgroup_node = None
for n in tree.nodes():
    if (n.taxon is not None) and (not any(x in taxa for x in
                                          n.taxon.label.split())):
        outgroup_node = n
        break
if outgroup_node is not None:
    print('found outgroup:', outgroup_node)
    tree.to_outgroup_position(outgroup_node, update_bipartitions=True)
tree2 = Phylo.read(io.StringIO(tree.as_string(schema='newick')), 'newick')

Phylo.draw_ascii(tree2)

