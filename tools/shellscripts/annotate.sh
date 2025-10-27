#!/bin/bash
#PBS -N Annotate
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log
#PBS -e op.err
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

set -e

TOOL=HG38
DIR_PATH="../../../Billi_data/Data"
DATA="$DIR_PATH/$TOOL/$TOOL.gfa"
DATA_VG="$DIR_PATH/$TOOL/hprc/$TOOL.vg"
DATA_PATH="$DIR_PATH/$TOOL/hprc/$TOOL-path.txt"
DATA_IDX="$DIR_PATH/$TOOL/hprc/$TOOL.xg"
DATA_VCF="$DIR_PATH/$TOOL/hprc/$TOOL.vcf"

# BUBBLE_PATH="$DIR_PATH/$TOOL/hprc/bubble.txt"

# /usr/bin/time -v gfatools bubble $DATA > $BUBBLE_PATH

# vg convert -g $DATA > $DATA_VG

# vg paths -L -x $DATA_VG > $DATA_PATH

# vg index $DATA_VG -x $DATA_IDX

vg deconstruct $DATA -p x -n > $DATA_VCF

# vg deconstruct $DATA_IDX -e -p "GRCh38#0#chrUn_KI270389v1" > $DATA_VCF