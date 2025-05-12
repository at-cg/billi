#!/bin/bash
#PBS -N Fast-PGGB
#PBS -l nodes=1:ppn=48
#PBS -q largemem
#PBS -l walltime=192:00:00
#PBS -o log_pggb
#PBS -e err_pggb
 
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
/usr/bin/time -v ./main /home/daanish/projects/Pangene/Data/PGGB/hprc-v1.0-pggb.gfa /home/daanish/projects/Pangene/results/Billi/PGGB
