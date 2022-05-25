#!/usr/bin/env python
# vim: ai

import os.path
import glob
import itertools
import sys
import subprocess

n = int(sys.argv[1])
N = int(sys.argv[2])

for infn in itertools.islice( sorted(glob.iglob("input/*.fasta")), n, None, N ):
    outfn = infn.replace(".fasta", "").replace("input", "output")
    if not os.path.exists("{}.hhblits".format(outfn)):
        subprocess.call(
            "/home/ubelix/dbmr/gg19u692/hhsuite-3.3.0-SSE2/bin/hhblits "
            "-cpu 1 -v 2 -i {0} -e 1e-3 "
            "-d /home/ubelix/dbmr/gg19u692/dbs/UniRef30_2020_06/UniRef30_2020_06 "
            "-o {1}.hhblits -oa3m {1}.a3m -n 3 -mact 0.35 >{1}.stdout 2>{1}.stderr".format(infn, outfn), shell=True)
    if not os.path.exists("{}.hhm".format(outfn)):
        subprocess.call(
            "/home/ubelix/dbmr/gg19u692/hhsuite-3.3.0-SSE2/bin/hhmake "
            "-i {1}.a3m -o {1}.hhm -v 2 >>{1}.stdout 2>>{1}.stderr".format(infn, outfn), shell=True)

#            "HHLIB=/home/ubelix/dbmr/gg19u692/hhsuite-3.3.0-SSE2/lib/hh "
