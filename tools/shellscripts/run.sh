#!/bin/bash
#PBS -N Run
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
# OUT_DIR1=$OUT_DIR/BubbleGun
# if [ ! -d "$OUT_DIR1" ]; then
#     mkdir -p $OUT_DIR1
# fi

# BubbleGun -g $FULL_PATH bchains --bubble_json "$OUT_DIR1/$FILE_NAME.json"
# rm -rf "BubbleGun.${FILE_NAME}.log"

# echo "Done running BubbleGun for input $FILE_NAME" >> $LOG


## vg
OUT_DIR2="$OUT_DIR/vg"
if [ ! -d "$OUT_DIR2" ]; then
    mkdir -p $OUT_DIR2
fi

vg convert -g $FULL_PATH > "$FILE_NAME.vg"
vg index -x "$FILE_NAME.xg" "$FILE_NAME.vg"

## Choose the right method and output before running
vg snarls -A cactus -a "$FILE_NAME.xg" > "$OUT_DIR2/$FILE_NAME.snarls"

vg view -j -R "$OUT_DIR2/$FILE_NAME.snarls" > "$OUT_DIR2/${FILE_NAME}_snarls.json"

rm -rf "$FILE_NAME.vg" "$FILE_NAME.xg" "$OUT_DIR2/$FILE_NAME.snarls"

echo "Done running vg for input $FILE_NAME" >> $LOG


## povu
# OUT_DIR3="$OUT_DIR/povu/$FILE_NAME"
# if [ ! -d "$OUT_DIR3" ]; then
#     mkdir -p $OUT_DIR3
# fi

# povu deconstruct -v 2 -h -t 48 -i $FULL_PATH -o $OUT_DIR3

# echo "Done running povu for input $FILE_NAME" >> $LOG


## Billi
OUT_DIR4="$OUT_DIR/Billi/$FILE_NAME"
if [ ! -d "$OUT_DIR4" ]; then
    mkdir -p $OUT_DIR4
fi

EXE=../../src/main
$EXE $FULL_PATH $OUT_DIR4

echo "Done running Billi for input $FILE_NAME" >> $LOG


## pangene
# OUT_DIR5="$OUT_DIR/pangene"
# if [ ! -f "$OUT_DIR5/$FILE_NAME.txt" ]; then
#     mkdir -p $OUT_DIR5
#     touch "$OUT_DIR5/$FILE_NAME.txt"
# fi

# cd ../../../Billi_data/pangene
# k8 pangene.js call $FULL_PATH > "../results/test/pangene/$FILE_NAME.txt"

# echo "Done running pangene for input $FILE_NAME" >> $LOG
