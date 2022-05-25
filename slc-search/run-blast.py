#!/usr/bin/env python
# vim: ai

from bs4 import BeautifulSoup
import re
import glob
import os.path
import subprocess


# run blast for each organism
for hits_fn in glob.glob('hits.*.fasta'):
    hits_fn_base = os.path.splitext(hits_fn)[0]
    subprocess.call('/opt/blast+/bin/makeblastdb '
                    '-dbtype prot '
                    '-in {0}.fasta'.format(hits_fn_base), shell=True)
    subprocess.call('/opt/blast+/bin/blastp '
                    '-db {0}.fasta '
                    '-query {0}.fasta '
                    '-out {0}.blast.xml '
                    '-outfmt 5 '
                    '-num_threads 8'.format(hits_fn_base), shell=True)

    # preprocess output

    rx_sp = re.compile(r'[|>]sp\|(\S+)\b')

    buf = open('{}.blast.xml'.format(hits_fn_base)).read()
    rx_iteration = re.compile(r'<Iteration>.*?</Iteration>', re.DOTALL)

    with open('{}.blast.out'.format(hits_fn_base), 'wt', encoding='utf-8') as f:
        for mo in rx_iteration.finditer(buf):
            soup = BeautifulSoup(mo.group(0), "lxml")
            for iteration in soup.find_all("iteration"):
                query_def = iteration.find("iteration_query-def").text
                query_name = query_def
                query_len = int(iteration.find("iteration_query-len").text)
                #print "query_len={:d}".format(query_len)
                for h in iteration.find_all("hit"):
                    hit_num = int(h.find("hit_num").text)
                    hit_id = h.find("hit_id").text
                    hit_def = h.find("hit_def").text
                    hit_accession = h.find("hit_accession").text
                    hit_len = int(h.find("hit_len").text)
                    if hit_def == query_def: continue
                    hsp_query_len = 0
                    hsp_hit_len = 0
                    identity = 0
                    #for hsp in h.find_all("hsp"):
                    hsp = h.find("hsp")
                    hsp_query_len += int(h.find("hsp_query-to").text)-int(h.find("hsp_query-from").text)+1
                    hsp_hit_len += int(h.find("hsp_hit-to").text)-int(h.find("hsp_hit-from").text)+1
                    identity += int(h.find("hsp_identity").text)

                    query_coverage = float(hsp_query_len)/float(query_len)
                    hit_coverage = float(hsp_hit_len)/float(hit_len)
                    query_identity_percent = float(identity)/float(query_len)
                    hsp_identity = float(identity)/float(hsp_query_len)

                    print("\t".join([str(x) for x in [
                        query_def, query_len, hit_def, hit_len, hsp_query_len,
                        hsp_hit_len, identity]]), file=f)

