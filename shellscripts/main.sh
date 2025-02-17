#!/bin/bash
#PBS -N Testing
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log
#PBS -e op.err
 
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
# time ./main /scratch/projects/daanish/data/Bubbles/Data/HPRC/GRCh38/hprc-v1.1-mc-grch38.gfa /home/daanish/projects/Pangene/results/Billi/HG38

# CHM13
# time ./main /scratch/projects/daanish/data/Bubbles/Data/HPRC/CHM13/hprc-v1.1-mc-chm13.gfa /home/daanish/projects/Pangene/results/Billi/CHM13

# PGGB
# time ./main /scratch/projects/daanish/data/Bubbles/Data/PGGB/hprc-v1.0-pggb.gfa /home/daanish/projects/Pangene/results/Billi/PGGB

## Others

# C4-90
# time ./main /scratch/projects/daanish/data/Bubbles/Data/C4/C4-90.gfa /home/daanish/projects/Pangene/results/Billi/C4-90

# MHC
# time ./main /scratch/projects/daanish/data/Bubbles/Data/MHC/MHC.gfa /home/daanish/projects/Pangene/results/Billi/MHC

# E-coli
# time ./main /scratch/projects/daanish/data/Bubbles/Data/Ecoli/Ecoli.gfa /home/daanish/projects/Pangene/results/Billi/Ecoli

# Human100
# time ./main /scratch/projects/daanish/data/Bubbles/Data/Human/human100.gfa /home/daanish/projects/Pangene/results/Billi/human100

# Human100p10
# time ./main /scratch/projects/daanish/data/Bubbles/Data/Human/human100p10.gfa /home/daanish/projects/Pangene/results/Billi/human100p10

# Mtb152m-p0a1
# time ./main /scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152m-p0a1.gfa /home/daanish/projects/Pangene/results/Billi/Mtb152m-p0a1

# Mtb152m-p1a2
# time ./main /scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152m-p1a2.gfa /home/daanish/projects/Pangene/results/Billi/Mtb152m-p1a2

# Mtb152p
# time ./main /scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152p.gfa /home/daanish/projects/Pangene/results/Billi/Mtb152p

## Testing 

# time ./main ../test/data/gfa_files/worst_case\(100000\).gfa ../test/data/results/worst_case\(100000\)
# time ./main ../test/data/gfa_files/EC16.gfa ../test/data/results/EC16


# DIR=../test/data/
# DDIR=$DIR/gfa_files
# RDIR=$DIR/results

# if [ ! -d "$RDIR" ]; then
#     mkdir -p $RDIR
# fi

# for FILE in "$DDIR"/*; do
#     if [ -f "$FILE" ]; then
#         FILENAME=$(basename "$FILE")
#         # echo "$FILENAME"
#         IFS='.' read -ra PARTS <<< "$FILENAME"
#         TAG="${PARTS[0]}"
#         # if [[ "$TAG" = "EC1" ]] || [[ "$TAG" = "EC2" ]]; then
#         #     continue
#         # fi
#         echo "$TAG"
#         ./main $DDIR/$FILENAME $RDIR/$TAG
#     fi
# done
