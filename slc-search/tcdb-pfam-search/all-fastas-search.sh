#!/bin/bash

hmmsearch --cpu 8 --domtblout all-fastas-search.domtblout --tblout all-fastas-search.tblout -o all-fastas-search.out ../../pfam-domains/pfam-20190214/Pfam-A.hmm.gz all-fastas.fasta 2>all-fastas-search.err

