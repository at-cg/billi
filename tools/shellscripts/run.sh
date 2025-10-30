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

## Compacting the graph
# if [ ! -d "$(dirname $FULL_PATH)/BubbleGun" ]; then
#     mkdir -p "$(dirname $FULL_PATH)/BubbleGun"
# fi

# COMPACT_PATH="$(dirname $FULL_PATH)/BubbleGun/${FILE_NAME}_compact.gfa"

# BubbleGun -g $FULL_PATH compact $COMPACT_PATH

# BubbleGun -g $FULL_PATH bchains --bubble_json "$OUT_DIR1/$FILE_NAME.json"
# BubbleGun -g $COMPACT_PATH bchains --bubble_json "$OUT_DIR1/${FILE_NAME}_compact.json"

# rm -rf "BubbleGun.${FILE_NAME}.log"

# echo "Done running BubbleGun for input $FILE_NAME" >> $LOG


## vg
# OUT_DIR2="$OUT_DIR/vg"
# if [ ! -d "$OUT_DIR2" ]; then
#     mkdir -p $OUT_DIR2
# fi

# TYPE=('' '_mod')
# BTYPE=('' '_ultrabubble')
# MAXNODES=1000000000

# vg convert -g $FULL_PATH > "$FILE_NAME.vg"
# vg mod -n "$FILE_NAME.vg" > "${FILE_NAME}_mod.vg"

# for TY in "${TYPE[@]}"; do

#     vg index -x "${FILE_NAME}$TY.xg" "${FILE_NAME}$TY.vg"

#     ## Choose the right method and output before running (-A cactus)
#     for BTY in "${BTYPE[@]}"; do
#         if [[ "$BTY" == "_ultrabubble" ]]; then
#             vg snarls -o -m $MAXNODES "${FILE_NAME}$TY.xg" > "$OUT_DIR2/${FILE_NAME}${TY}${BTY}.snarls"
#         else
#             vg snarls -a -m $MAXNODES "${FILE_NAME}$TY.xg" > "$OUT_DIR2/${FILE_NAME}${TY}${BTY}.snarls"
#         fi

#         vg view -j -R "$OUT_DIR2/${FILE_NAME}${TY}${BTY}.snarls" > "$OUT_DIR2/${FILE_NAME}${TY}${BTY}_snarls.json"

#         rm -rf "$OUT_DIR2/${FILE_NAME}${TY}${BTY}.snarls"

#         echo "Done running vg for input ${FILE_NAME}${TY}${BTY}" >> $LOG
#     done

#     rm -rf "${FILE_NAME}$TY.vg" "${FILE_NAME}$TY.xg"

# done


## povu
# OUT_DIR3="$OUT_DIR/povu/$FILE_NAME"
# if [ ! -d "$OUT_DIR3" ]; then
#     mkdir -p $OUT_DIR3
# fi

# povu decompose -h -s -i $FULL_PATH -o $OUT_DIR3

# echo "Done running povu for input $FILE_NAME" >> $LOG


## panbubble
OUT_DIR4="$OUT_DIR/panbubble_stress/$FILE_NAME"
if [ ! -d "$OUT_DIR4" ]; then
    mkdir -p $OUT_DIR4
fi

EXE=../../src/main
$EXE decompose -i $FULL_PATH -c -r -f 1 -o $OUT_DIR4

echo "Done running panbubble for input $FILE_NAME" >> $LOG


## pangene
# OUT_DIR5="$OUT_DIR/pangene"
# if [ ! -f "$OUT_DIR5/$FILE_NAME.txt" ]; then
#     mkdir -p $OUT_DIR5
#     touch "$OUT_DIR5/$FILE_NAME.txt"
# fi

# cd ../../../panbubble_data/pangene
# k8 pangene.js call $FULL_PATH > "../results/small/pangene/$FILE_NAME.txt"

# echo "Done running pangene for input $FILE_NAME" >> $LOG
