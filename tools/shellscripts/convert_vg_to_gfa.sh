#!/bin/bash
#PBS -N Convert_HG38
#PBS -l nodes=1:ppn=48
#PBS -q long
#PBS -l walltime=72:00:00
#PBS -o op.log
#PBS -e op.err
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

set -e

# REFERENCE=("CHM13" "HG38")
REFERENCE=("HG38")

for REF in "${REFERENCE[@]}"; do

    DIR="../../../panbubble_data/Data/$REF/v2/chromosomes"

    for FILE in "$DIR"/*.vg; do
        if [ -f "$FILE" ]; then

            base="${FILE%.*}"

            vg index $FILE -x "$base.xg"
            vg convert -f "$base.xg" > "$base.gfa"

            rm -rf "$base.xg"

        fi
    done

done