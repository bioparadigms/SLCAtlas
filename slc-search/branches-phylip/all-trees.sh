#!/bin/bash

NJOBS=4

# conda activate pyweb

for fn in branch_*[0-9a-f].fasta ; do
	if [ $(jobs | wc -l) -ge $NJOBS ]; then wait -n ; fi
	echo $fn
	NSEQ=$(grep -c '^>' $fn)
	if [ $NSEQ -le 1 ]; then
		echo "single sequence in fasta file, ignoring"
	elif [ $NSEQ -le 2 ]; then
		echo "two sequences in fasta file, only clustalo"
		aln_fn=${fn/.fasta/.aln.fasta}
		if [ ! -s $aln_fn ]; then 
			echo "clustalo -i $fn --iter=5 -o $aln_fn &"
			clustalo -i $fn --iter=5 -o $aln_fn &
		fi
	else
		#python single-tree.py $fn >${fn/.fasta/.aln.tree.txt} &
		tree_fn="${fn/.fasta/.aln.sms}/${fn/.fasta/.aln.phylip_phyml_tree.renamed.txt}"
		if [ ! -s "$tree_fn" ]; then 
			echo "python -u single-tree2.alt.py $fn 2>${fn/.fasta/.single-tree2.err} | tee ${fn/.fasta/.single-tree2.out} &"
			python -u single-tree2.alt.py $fn 2>${fn/.fasta/.single-tree2.err} | tee ${fn/.fasta/.single-tree2.out} &
		fi
	fi
done

wait

