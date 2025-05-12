#!/bin/bash
#PBS -N Main-CACTUS
#PBS -l nodes=1:ppn=48
#PBS -q largemem
#PBS -l walltime=192:00:00
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
/usr/bin/time -v ./main /scratch/projects/daanish/data/Bubbles/Data/HPRC/GRCh38/hprc-v1.1-mc-grch38.gfa /home/daanish/projects/Pangene/results/Billi/HG38

# CHM13
/usr/bin/time -v ./main /scratch/projects/daanish/data/Bubbles/Data/HPRC/CHM13/hprc-v1.1-mc-chm13.gfa /home/daanish/projects/Pangene/results/Billi/CHM13
