#!/bin/bash

for fastafn in ../slc-fastas/*/*.fasta ; do
	acc=$(basename $fastafn)
	acc=${acc%.fasta}
	if [ ! -e "output/${acc}.hhr" ]; then
		echo "$acc"
	fi
done
	
