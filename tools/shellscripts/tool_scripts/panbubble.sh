#!/bin/bash
#PBS -N Run_panbubble
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

## panbubble
OUT_DIR4="$OUT_DIR/panbubble_correct/$FILE_NAME"
if [ ! -d "$OUT_DIR4" ]; then
    mkdir -p $OUT_DIR4
fi

MAXDEPTH=1000000000
EXE=../../src/main_brute
/usr/bin/time -v $EXE decompose -i $FULL_PATH -d $MAXDEPTH -c -r -p -o $OUT_DIR4

echo "Done running panbubble for input $FILE_NAME" >> $LOG
