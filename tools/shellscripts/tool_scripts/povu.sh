#!/bin/bash
#PBS -N Run_povu
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

if [ $# -ne 4 ]; then
    echo "full_file_path, out_dir, file_name, log"
    exit 1
fi

FULL_PATH=$1
OUT_DIR=$2
FILE_NAME=$3
LOG=$4

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

# povu
OUT_DIR3="$OUT_DIR/povu/$FILE_NAME"
if [ ! -d "$OUT_DIR3" ]; then
    mkdir -p $OUT_DIR3
fi

povu decompose -h -s -i $FULL_PATH -o $OUT_DIR3

echo "Done running povu for input $FILE_NAME" >> $LOG
