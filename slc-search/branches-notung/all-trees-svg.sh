#!/bin/bash

# conda activate pyweb

for fn in branch_*[0-9a-f] ; do
	echo $fn
	tree_fn=$fn/${fn}.tree.rooting.0.rearrange.0.reconciled
	if [ ! -e $tree_fn ]; then
		tree_fn=$fn/${fn}.tree.rooting.0.rearrange.0
		if [ ! -e $tree_fn ]; then
			tree_fn=$fn/${fn}.tree.rooting.0
			if [ ! -e $tree_fn ]; then
				echo "suitable tree not found: $fn"
				continue
			fi
		fi
	fi
	python -u make-tree-svg.py $tree_fn 2>$fn/make-tree-svg.err >$fn/make-tree-svg.out
	# this is done in the make-slc-family-trees.py script
	#inkscape -z -A $fn/${fn}.pdf $fn/make-tree-svg.out
done

