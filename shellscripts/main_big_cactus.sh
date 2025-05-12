#!/bin/bash
#PBS -N Main-CACTUS
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log_big_cactus
#PBS -e op.err_big_cactus
 
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

# HG38
/usr/bin/time -v ./main /home/daanish/projects/Pangene/Data/HG38/hprc-v1.1-mc-grch38.gfa /home/daanish/projects/Pangene/results/Billi/HG38

# CHM13
/usr/bin/time -v ./main /home/daanish/projects/Pangene/Data/CHM13/hprc-v1.1-mc-chm13.gfa /home/daanish/projects/Pangene/results/Billi/CHM13
