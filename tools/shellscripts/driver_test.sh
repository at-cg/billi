#!/bin/bash
#PBS -N Driver-Test
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

DATA_DIR="../../test/data/gfa_files"
OUT_DIR="../../../Billi_data/results/test"
LOG=log-test.txt

echo "Starting Driver-Test Script" > $LOG

SHELL_EXE=run.sh
chmod +x $SHELL_EXE

MAKE_FILE=../../src
make -C $MAKE_FILE clean
make -C $MAKE_FILE
ulimit -s unlimited

for FILE in "$DATA_DIR"/*; do
    if [ -f "$FILE" ]; then
        REAL_PATH=$(realpath $FILE)
        # echo "$REAL_PATH" >&2

        FILENAME=$(basename "$FILE")
        IFS='.' read -ra PARTS <<< "$FILENAME"
        TAG="${PARTS[0]}"
        # if [[ "$TAG" = "worst_case(100000)" ]]; then
        #     continue
        # fi
        if [[ "$TAG" != "EC19" ]]; then
            continue
        fi
        # echo "$TAG"

        sh $SHELL_EXE $REAL_PATH $OUT_DIR $TAG $LOG
    fi
done
