#!/bin/bash
#PBS -N Benchmark
#PBS -l nodes=1:ppn=48
#PBS -q largemem
#PBS -l walltime=192:00:00
#PBS -o benchmark.log
#PBS -e benchmark.err
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

set -e

DATA_DIR="../../test/data/gfa_files"
OUT_DIR="../../../panbubble_data/results/benchmark"

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

LOG=log-benchmark.txt

echo "Starting Driver-Benchmark Script" >> $LOG

# DATA=("/scratch/projects/daanish/data/Bubbles/Data/C4/C4-90-modified.gfa"
#       "/scratch/projects/daanish/data/Bubbles/Data/MHC/MHC-modified.gfa"
#       "/scratch/projects/daanish/data/Bubbles/Data/Human/human100.gfa"
#       "/scratch/projects/daanish/data/Bubbles/Data/Human/human100p10.gfa"
#       "/scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152m-p0a1.gfa"
#       "/scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152m-p1a2.gfa"
#       "/scratch/projects/daanish/data/Bubbles/Data/Human/Mtb152p.gfa"
#       "/home/daanish/projects/panbubble_data/Data/CHM13/v1/CHM13.gfa"
#       "/home/daanish/projects/panbubble_data/Data/HG38/v1/HG38.gfa"
#       "/home/daanish/projects/panbubble_data/Data/CHM13/v2/hprc-v2.0-mc-chm13.gfa"
#       "/home/daanish/projects/panbubble_data/Data/HG38/v2/hprc-v2.0-mc-grch38.gfa")
    #   ""
    #   )
DATA=("/home/daanish/projects/panbubble_data/Data/CHM13/v1/CHM13.gfa"
      "/home/daanish/projects/panbubble_data/Data/HG38/v1/HG38.gfa")

# TAG=("Test")
# TAG=("C4-90" "MHC" "human100" "human100p10" "Mtb152m-p0a1" "Mtb152m-p1a2" "Mtb152p" "CHM13_v1" "HG38_v1" "CHM13_v2" "HG38_v2")
TAG=("CHM13_v1" "HG38_v1")
    #     "" ""
    #     "" ""
    # )
     
MAKE_FILE=../../src
make -C $MAKE_FILE clean
make -C $MAKE_FILE
ulimit -s unlimited

SHELL_EXE1=tool_scripts/panbubble.sh
SHELL_EXE2=tool_scripts/vg.sh
SHELL_EXE3=tool_scripts/pangene.sh

chmod +x $SHELL_EXE1 $SHELL_EXE2 $SHELL_EXE3

sz=${#DATA[@]}

## panbubble
for ((i=0; i<=$((sz-1)); i++)); do
    
    sh $SHELL_EXE1 ${DATA[$i]} $OUT_DIR ${TAG[$i]} $LOG

done

## vg
# for ((i=0; i<=$((sz-1)); i++)); do
    
#     sh $SHELL_EXE2 ${DATA[$i]} $OUT_DIR ${TAG[$i]} $LOG

# done

## pangene
# for ((i=0; i<=$((sz-1)); i++)); do
    
#     sh $SHELL_EXE3 ${DATA[$i]} $OUT_DIR ${TAG[$i]} $LOG

# done