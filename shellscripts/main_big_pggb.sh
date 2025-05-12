#!/bin/bash
#PBS -N Main-PGGB
#PBS -l nodes=1:ppn=48
#PBS -q largemem
#PBS -l walltime=192:00:00
#PBS -o op.log_big_pggb1
#PBS -e op.err_big_pggb1
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

set -e

cd ../src
make clean
make
ulimit -s unlimited

## HPRC

# PGGB
/usr/bin/time -v ./main /scratch/projects/daanish/data/Bubbles/Data/PGGB/hprc-v1.0-pggb.gfa /home/daanish/projects/Pangene/results/Billi/PGGB
