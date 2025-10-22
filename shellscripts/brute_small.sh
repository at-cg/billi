#!/bin/bash
#PBS -N Brute-small
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log_brute_small
#PBS -e op.err_brute_small
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

set -e

cd ../src
make -f makefile_brute clean
make -f makefile_brute
ulimit -s unlimited


# C4-90
/usr/bin/time -v ./main_brute decompose -i /scratch/projects/daanish/data/Bubbles/Data/C4/C4-90.gfa -o /home/daanish/projects/panbubble/test/data/results_brute/C4-90 -c 

# MHC
/usr/bin/time -v ./main_brute decompose -i /scratch/projects/daanish/data/Bubbles/Data/MHC/MHC.gfa -o /home/daanish/projects/panbubble/test/data/results_brute/MHC -c 

# E-coli
/usr/bin/time -v ./main_brute decompose -i /scratch/projects/daanish/data/Bubbles/Data/Ecoli/Ecoli.gfa -o /home/daanish/projects/panbubble/test/data/results_brute/Ecoli -c 

# Human100
/usr/bin/time -v ./main_brute decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/human100.gfa -o /home/daanish/projects/panbubble/test/data/results_brute/human100 -c

# Human100p10
/usr/bin/time -v ./main_brute decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/human100p10.gfa -o /home/daanish/projects/panbubble/test/data/results_brute/human100p10 -c 

# Mtb152m-p0a1
/usr/bin/time -v ./main_brute decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152m-p0a1.gfa -o /home/daanish/projects/panbubble/test/data/results_brute/Mtb152m-p0a1 -c 

# Mtb152m-p1a2
/usr/bin/time -v ./main_brute decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152m-p1a2.gfa -o /home/daanish/projects/panbubble/test/data/results_brute/Mtb152m-p1a2 -c

# Mtb152p
/usr/bin/time -v ./main_brute decompose -i /scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152p.gfa -o /home/daanish/projects/panbubble/test/data/results_brute/Mtb152p -c

