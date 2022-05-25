#!/bin/bash

# n=0
for i in ../branches-phylip/branch_*.aln.sms ; do
    branch=$(basename $i | sed -e 's/\.aln\.sms$//')
    echo $branch
    mkdir -p $branch
    cd $branch
    if [ ! -s "${branch}.tree" ]; then
        sed -e 's/_\(Danio\|Homo\|Gallus\|Rattus\|Mus\|Caenorhabditis\|Drosophila\)_\([0-9A-Z]\+\)/_\2_\1/g' \
            -e 's/zgc:/zgc./g' \
            -e 's/si:/si./g' \
            -e 's/NEST:/NEST./g' \
            -e 's/im:/im./g' \
            -e 's/anon-EST:/anon-EST./g' \
            -e 's/BcDNA:/BcDNA./g' \
            -e 's/EG:/EG./g' \
            -e 's/zmp:/zmp./g' \
            -e 's/(2)/-2-/g' \
            ../$i/${branch}.aln.phylip_phyml_tree.renamed.txt \
            >${branch}.tree
    fi
    if [ ! -s "${branch}.tree.rooting.ntglog" ]; then
        /opt/jre1.8.0_92/bin/java -jar ~/Notung-2.9.1.5/Notung-2.9.1.5.jar \
            -g ${branch}.tree \
            -s ../7species.stree \
            --root \
            --prune \
            --speciestag postfix \
            --maxtrees 20 --allopt \
            --log \
            --parsable \
            --treestats \
            --rootscores 
    fi
    for j in ${branch}.tree.rooting.{0..19} ; do
        if [ ! -e $j ]; then break ; fi
        if [ ! -s "$j.rearrange.ntglog" ]; then
            /opt/jre1.8.0_92/bin/java -jar ~/Notung-2.9.1.5/Notung-2.9.1.5.jar \
                -g $j \
                --rearrange \
                --prune \
                --speciestag postfix \
                --maxtrees 20 \
                --threshold 0.9 \
                --log \
                --parsable \
                --treestats
        fi
        for k in $j.rearrange.{0..19} ; do
            if [ ! -e $k ]; then break ; fi
            if [ ! -s "$k.reconciled.ntglog" ]; then
                /opt/jre1.8.0_92/bin/java -jar ~/Notung-2.9.1.5/Notung-2.9.1.5.jar \
                    -g $k \
                    --reconcile \
                    --prune \
                    --speciestag postfix \
                    --maxtrees 20 \
                    --log \
                    --parsable \
                    --treestats \
                    --homologtabletabs
            fi
        done
    done
    cd ..
    # n=$((n+1))
    # if [ $n -ge 5 ]; then break ; fi
done

