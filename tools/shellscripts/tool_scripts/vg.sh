#!/bin/bash
#PBS -N Run_vg
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

## vg
OUT_DIR2="$OUT_DIR/vg"
if [ ! -d "$OUT_DIR2" ]; then
    mkdir -p $OUT_DIR2
fi

# TYPE=('' '_mod')
# BTYPE=('' '_ultrabubble')
TYPE=('')
BTYPE=('')
MAXNODES=1000000000

vg convert -g $FULL_PATH > "$FILE_NAME.vg"
vg mod -n "$FILE_NAME.vg" > "${FILE_NAME}_mod.vg"

for TY in "${TYPE[@]}"; do

    vg index -x "${FILE_NAME}$TY.xg" "${FILE_NAME}$TY.vg"

    ## Choose the right method and output before running (-A cactus)
    for BTY in "${BTYPE[@]}"; do
        if [[ "$BTY" == "_ultrabubble" ]]; then
            /usr/bin/time -v vg snarls -t 48 -o -m $MAXNODES "${FILE_NAME}$TY.xg" > "$OUT_DIR2/${FILE_NAME}${TY}${BTY}.snarls"
        else
            /usr/bin/time -v vg snarls -t 48 -a -m $MAXNODES "${FILE_NAME}$TY.xg" > "$OUT_DIR2/${FILE_NAME}${TY}${BTY}.snarls"
        fi

        vg view -j -R "$OUT_DIR2/${FILE_NAME}${TY}${BTY}.snarls" > "$OUT_DIR2/${FILE_NAME}${TY}${BTY}_snarls.json"

        rm -rf "$OUT_DIR2/${FILE_NAME}${TY}${BTY}.snarls"

        echo "Done running vg for input ${FILE_NAME}${TY}${BTY}" >> $LOG
    done

    rm -rf "${FILE_NAME}$TY.vg" "${FILE_NAME}$TY.xg"

done
