#!/bin/bash
#PBS -N Convert
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

DIR='../../../Billi_data/Data/CHM13/chromosomes'

vg index "$DIR/chr1.vg" -x "$DIR/chr1.xg"
vg convert -f "$DIR/chr1.xg" > "$DIR/chr1.gfa"
