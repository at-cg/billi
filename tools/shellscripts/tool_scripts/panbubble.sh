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

if [ $# -ne 5 ]; then
    echo "full_file_path, out_dir, file_name, offset, log"
    exit 1
fi

FULL_PATH=$1
OUT_DIR=$2
FILE_NAME=$3
OFFSET=$4
LOG=$5

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p $OUT_DIR
fi

## panbubble
OUT_DIR4="$OUT_DIR/panbubble_debug/$FILE_NAME-$OFFSET"
if [ ! -d "$OUT_DIR4" ]; then
    mkdir -p $OUT_DIR4
fi

EXE=../../src/main1
/usr/bin/time -v $EXE decompose -i $FULL_PATH -f $OFFSET -c -r -o $OUT_DIR4

echo "Done running panbubble for input $FILE_NAME-$OFFSET" >> $LOG
