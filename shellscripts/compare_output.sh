#!/bin/bash
#PBS -N Compare
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

EXE=/home/daanish/projects/panbubble_scratch/src/src_new/comparison/compare_with_vg.py

python3 $EXE
