#!/bin/bash
#PBS -N Driver-Small
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

DATA_DIR="/scratch/projects/daanish/data/Bubbles/Data"
OUT_DIR="../../../panbubble_data/results/small"
LOG=log-small.txt

echo "Starting Driver-Small Script" > $LOG

SHELL_EXE=run.sh
chmod +x $SHELL_EXE

MAKE_FILE=../../src
make -C $MAKE_FILE clean
make -C $MAKE_FILE
ulimit -s unlimited

# C4-90
sh $SHELL_EXE "$DATA_DIR/C4/C4-90-modified.gfa" $OUT_DIR "C4-90" $LOG

# # MHC
sh $SHELL_EXE "$DATA_DIR/MHC/MHC-modified.gfa" $OUT_DIR "MHC" $LOG

# # E-coli
sh $SHELL_EXE "$DATA_DIR/Ecoli/Ecoli.gfa" $OUT_DIR "Ecoli" $LOG

# # Human100
sh $SHELL_EXE "$DATA_DIR/Human/human100.gfa" $OUT_DIR "human100" $LOG

# # Human100p10
sh $SHELL_EXE "$DATA_DIR/Human/human100p10.gfa" $OUT_DIR "human100p10" $LOG

# # Mtb152m-p0a1
sh $SHELL_EXE "$DATA_DIR/Human/Mtb152m-p0a1.gfa" $OUT_DIR "Mtb152m-p0a1" $LOG

# # Mtb152m-p1a2
sh $SHELL_EXE "$DATA_DIR/Human/Mtb152m-p1a2.gfa" $OUT_DIR "Mtb152m-p1a2" $LOG

# # Mtb152p
sh $SHELL_EXE "$DATA_DIR/Human/Mtb152p.gfa" $OUT_DIR "Mtb152p" $LOG
