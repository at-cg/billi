#!/bin/bash
#PBS -N Gunzip
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

DIR='../../../panbubble_data/Data/CHM13/v2'

cd $DIR

gunzip /home/daanish/projects/panbubble_data/Data/CHM13/v2/hprc-v2.0-mc-chm13.gfa.gz

cd ../../HG38/v2

gunzip /home/daanish/projects/panbubble_data/Data/HG38/v2/hprc-v2.0-mc-grch38.gfa.gz
