#!/bin/bash
#PBS -N Benchmark
#PBS -l nodes=1:ppn=48
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -o benchmark_final1.log
#PBS -e benchmark_final1.err
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

set -e

OUT_DIR="../../../panbubble_data/results/benchmark_final"

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

LOG=log-benchmark.txt

echo "Starting Driver-Benchmark Script" >> $LOG

DATA=("/scratch/projects/daanish/data/Bubbles/Data/Ecoli/Ecoli.gfa"
      "/scratch/projects/daanish/data/Bubbles/Data/Human/human100p10.gfa"
      "/scratch/projects/daanish/data/Bubbles/Data/C4/C4-90-modified.gfa"
      "/scratch/projects/daanish/data/Bubbles/Data/MHC/MHC-modified.gfa"
      "/home/daanish/projects/panbubble_data/Data/CHM13/v2/chromosomes/chrY.gfa"
      "/home/daanish/projects/panbubble_data/Data/CHM13/v2/chromosomes/chrX.gfa")

TAG=("Ecoli" "human100p10" "C4-90" "MHC" "CHM13_v2_Y" "CHM13_v2_X")
     
MAKE_FILE=../../src
make -C $MAKE_FILE clean
make -C $MAKE_FILE
ulimit -s unlimited

SHELL_EXE1=tool_scripts/panbubble.sh
SHELL_EXE2=tool_scripts/vg.sh
SHELL_EXE3=tool_scripts/pangene.sh

chmod +x $SHELL_EXE1 $SHELL_EXE2 $SHELL_EXE3

sz=${#DATA[@]}

start=100000000
end=100000000

## panbubble
# for ((i=0; i<=$((sz-1)); i++)); do
#     if (( i > 3 )); then 
#         start=1
#         end=10
#     fi

#     for ((j=$start; j<=$end; j++)); do
#         sh $SHELL_EXE1 ${DATA[$i]} $OUT_DIR ${TAG[$i]} $j $LOG
#     done
# done

## vg
for ((i=0; i<=$((sz-1)); i++)); do
    
    sh $SHELL_EXE2 ${DATA[$i]} $OUT_DIR ${TAG[$i]} $LOG

done

## pangene
# for ((i=0; i<=$((sz-1)); i++)); do
    
#     sh $SHELL_EXE3 ${DATA[$i]} $OUT_DIR ${TAG[$i]} $LOG

# done