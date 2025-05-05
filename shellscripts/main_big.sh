#!/bin/bash
#PBS -N Main-big
#PBS -l nodes=1:ppn=48
#PBS -q largemem
#PBS -l walltime=24:00:00
#PBS -o op.log_big
#PBS -e op.err_big
 
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

# PGGB
/usr/bin/time -v ./main /scratch/projects/daanish/data/Bubbles/Data/PGGB/hprc-v1.0-pggb.gfa /home/daanish/projects/Pangene/results/Billi/PGGB
