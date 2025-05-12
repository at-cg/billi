#!/bin/bash
#PBS -N Main-PGGB-D
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log_big_pggb2
#PBS -e op.err_big_pggb2
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

set -e

cd ../src
make -f makefile_debug clean
make -f makefile_debug
ulimit -s unlimited

## HPRC

# PGGB
/usr/bin/time -v ./main_debug /home/daanish/projects/Pangene/Data/PGGB/hprc-v1.0-pggb.gfa /home/daanish/projects/Pangene/results/Billi/PGGB

