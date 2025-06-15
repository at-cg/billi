#!/bin/bash
#PBS -N Test-Annotate
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log-test
#PBS -e op.err-test
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

set -e

TOOL=TEST
CHR=chrM
DIR_PATH=/home/daanish/projects/Pangene/Data
# DATA="$DIR_PATH/$TOOL/$TOOL.gfa"
DATA_VG="$DIR_PATH/$TOOL/hprc/$CHR.vg"
DATA_PATH="$DIR_PATH/$TOOL/hprc/$CHR-path.txt"
DATA_IDX="$DIR_PATH/$TOOL/hprc/$CHR.xg"
DATA_VCF="$DIR_PATH/$TOOL/hprc/$CHR.vcf"
FILT_VCF="$DIR_PATH/$TOOL/hprc/FILT-$CHR.vcf"

BUBBLE_PATH="$DIR_PATH/$TOOL/hprc/bubble.txt"

# /usr/bin/time -v gfatools bubble $DATA > $BUBBLE_PATH

# vg convert -g $DATA > $DATA_VG

# vg paths -L -x $DATA_VG > $DATA_PATH

# vg index $DATA_VG -x $DATA_IDX

## vg find -x $DATA_IDX -p "chr7" | vg view -a - | grep '"type":"Variant"' > "$DIR_PATH/log.txt"

vg deconstruct $DATA_IDX -P "HG04199#2#JBHDTI010000083.1#0[0-15992]" -a > $DATA_VCF
gzip $DATA_VCF # https://github.com/pangenome/vcfbub/issues/4

# vg deconstruct $DATA -P chr15 -a > $DATA_VCF
# -H "#" -P S288C -e -a -t 8

vcfbub -l 0 -r 10000000 -i "$DATA_VCF.gz" > $FILT_VCF
