#!/bin/bash
#PBS -N Main-small
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log_small
#PBS -e op.err_small
 
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

# C4-90
/usr/bin/time -v ./main decompose -i /scratch/projects/daanish/data/Bubbles/Data/C4/C4-90.gfa -d $DEPTH -c -r -p -o /home/daanish/projects/panbubble_data/results/small/panbubble/C4-90

# MHC
/usr/bin/time -v ./main decompose -i /scratch/projects/daanish/data/Bubbles/Data/MHC/MHC.gfa -d $DEPTH -c -r -p -o /home/daanish/projects/panbubble_data/results/small/panbubble/MHC

# E-coli
/usr/bin/time -v ./main decompose -i /scratch/projects/daanish/data/Bubbles/Data/Ecoli/Ecoli.gfa -d $DEPTH -c -r -p -o /home/daanish/projects/panbubble_data/results/small/panbubble/Ecoli

# Human100
/usr/bin/time -v ./main decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/human100.gfa -d $DEPTH -c -r -p -o /home/daanish/projects/panbubble_data/results/small/panbubble/human100

# Human100p10
/usr/bin/time -v ./main decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/human100p10.gfa -d $DEPTH -c -r -p -o /home/daanish/projects/panbubble_data/results/small/panbubble/human100p10

# Mtb152m-p0a1
/usr/bin/time -v ./main decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152m-p0a1.gfa -d $DEPTH -c -r -p -o /home/daanish/projects/panbubble_data/results/small/panbubble/Mtb152m-p0a1

# Mtb152m-p1a2
/usr/bin/time -v ./main decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152m-p1a2.gfa -d $DEPTH -c -r -p -o /home/daanish/projects/panbubble_data/results/small/panbubble/Mtb152m-p1a2

# Mtb152p
/usr/bin/time -v ./main decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152p.gfa -d $DEPTH -c -r -p -o /home/daanish/projects/panbubble_data/results/small/panbubble/Mtb152p

