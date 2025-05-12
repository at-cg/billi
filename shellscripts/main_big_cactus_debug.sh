#!/bin/bash
#PBS -N Main-Cactus-D
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o log_debug_cactus
#PBS -e err_debug_cactus
 
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

# CHM13
/usr/bin/time -v ./main_debug /home/daanish/projects/Pangene/Data/CHM13/hprc-v1.1-mc-chm13.gfa /home/daanish/projects/Pangene/results/Billi/CHM13

