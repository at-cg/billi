#!/bin/bash
#PBS -N printBubble
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log_print
#PBS -e op.err_print
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

# set -e

cd ../src
make -f makefile_print clean
make -f makefile_print
ulimit -s unlimited

# E-coli
/usr/bin/time -v ./printBubble /scratch/projects/daanish/data/Bubbles/Data/Ecoli/Ecoli.gfa /home/daanish/projects/Pangene/results/Billi/Ecoli
