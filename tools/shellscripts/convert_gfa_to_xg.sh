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

FILE='/home/daanish/projects/Billi/test/data/gfa_files/EC1.gfa'

vg convert -g $FILE > 'EC1.vg'
vg index -x 'EC1.xg' 'EC1.vg'