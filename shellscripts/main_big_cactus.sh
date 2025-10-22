#!/bin/bash
#PBS -N OPT-CACTUS
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o log_cactus
#PBS -e err_cactus
 
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

DEPTH=1000000000

## HPRC

# HG38
/usr/bin/time -v ./main decompose -i /home/daanish/projects/panbubble_data/Data/HG38/v1/HG38.gfa -d $DEPTH -c true -r true -p true -o /home/daanish/projects/Pangene/results/panbubble/HG38

# CHM13-V1
/usr/bin/time -v ./main decompose -i /home/daanish/projects/panbubble_data/Data/CHM13/v1/CHM13.gfa -d $DEPTH -c true -r true -p true -o /home/daanish/projects/Pangene/results/panbubble/CHM13

