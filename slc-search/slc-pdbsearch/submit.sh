#!/bin/bash
# vim: ai

HHSUITE=$HOME/hhsuite-3.3.0-SSE2
PDB70=$HOME/dbs/pdb70_from_mmcif_210804/pdb70

JOBNAME=slc-pdbsearch

INPUTDIR=$PWD/input
OUTPUTDIR=$PWD/output

# required for SGE
mkdir -p $OUTPUTDIR

FILES=()
for fn in $INPUTDIR/*.hhm ; do
	name=$(basename $fn)
	outfn="output/${name%.hhm}.hhr"
	if [ ! -e "$outfn" ]; then
		FILES+=("$fn")
	fi
done
#FILES=($INPUTDIR/*.hhm)
N=${#FILES[@]}

echo "$N"
#echo "${FILES[@]}"

sbatch --array=1-$N << _EOF_
#!/bin/bash

#SBATCH --mail-user=gergely.gyimesi@dbmr.unibe.ch
#SBATCH --mail-type=begin,fail,end
#SBATCH --output=$OUTPUTDIR/slurm-%j.out
#SBATCH --error=$OUTPUTDIR/slurm-%j.err
#SBATCH --requeue
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=$JOBNAME

FILES=(${FILES[*]})

FASTA=\${FILES[\$((\$SLURM_ARRAY_TASK_ID-1))]}
FASTA_BASE=\`basename \$FASTA\`
PROTNAME=\${FASTA_BASE/.hhm/}

# check if input dir is reachable
i=1 ; while [ ! -d $INPUTDIR -a "\$i" -le 24 ]; do echo "waiting for $INPUTDIR" ; sleep 600 ; i=\$((i+1)) ; done
if [ \$i -gt 24 ]; then echo "giving up." ; exit -1 ; fi

# check if output dir is reachable
i=1 ; while [ ! -d $OUTPUTDIR -a "\$i" -le 24 ]; do mkdir -p $OUTPUTDIR ; sleep 600 ; i=\$((i+1)) ; done
if [ \$i -gt 24 ]; then echo "cannot reach $OUTPUTDIR, giving up." ; exit -1 ; fi

cd $OUTPUTDIR

$HHSUITE/bin/hhsearch -cpu 1 -v 2 -i $INPUTDIR/\${PROTNAME}.hhm -d $PDB70 -o \${PROTNAME}.hhr -p 20 -P 20 -Z 100 -z 1 -b 1 -B 100 -seq 1 -aliw 80 -local -ssm 2 -norealign -sc 1 -dbstrlen 10000

_EOF_

