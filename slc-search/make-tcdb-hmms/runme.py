#!/usr/bin/env python

import os
import os.path
import sys
import glob
import subprocess
from Bio import AlignIO

for fn in glob.iglob("input/*.aln"):
    name = os.path.splitext(os.path.basename(fn))[0]
    sys.stderr.write("{}\n".format(name))
    msa = AlignIO.read(open(fn), "clustal")
    f = open("output/{}.fasta".format(name), "w")
    AlignIO.write(msa, f, "fasta")
    f.close()
    subprocess.call("hmmbuild -n {0} output/{0}.hmm output/{0}.fasta >output/{0}.hmmbuild.out".format(name), shell=True)

