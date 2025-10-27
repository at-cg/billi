#!/bin/bash
#PBS -N Run_bubblegun
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

## BubbleGun
OUT_DIR1=$OUT_DIR/BubbleGun
if [ ! -d "$OUT_DIR1" ]; then
    mkdir -p $OUT_DIR1
fi

# Compacting the graph
if [ ! -d "$(dirname $FULL_PATH)/BubbleGun" ]; then
    mkdir -p "$(dirname $FULL_PATH)/BubbleGun"
fi

COMPACT_PATH="$(dirname $FULL_PATH)/BubbleGun/${FILE_NAME}_compact.gfa"

BubbleGun -g $FULL_PATH compact $COMPACT_PATH

BubbleGun -g $FULL_PATH bchains --bubble_json "$OUT_DIR1/$FILE_NAME.json"
BubbleGun -g $COMPACT_PATH bchains --bubble_json "$OUT_DIR1/${FILE_NAME}_compact.json"

rm -rf "BubbleGun.${FILE_NAME}.log"

echo "Done running BubbleGun for input $FILE_NAME" >> $LOG