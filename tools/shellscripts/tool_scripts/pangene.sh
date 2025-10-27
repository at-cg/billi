#!/bin/bash
#PBS -N Run_pangene
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

echo "$(pwd)"

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

# pangene
OUT_DIR5="$OUT_DIR/pangene"
if [ ! -f "$OUT_DIR5/$FILE_NAME.txt" ]; then
    mkdir -p $OUT_DIR5
    touch "$OUT_DIR5/$FILE_NAME.txt"
fi

cd ../../../panbubble_data/pangene
/usr/bin/time -v k8 pangene.js call $FULL_PATH > "../results/benchmark/pangene/$FILE_NAME.txt"

echo "Done running pangene for input $FILE_NAME" >> $LOG
