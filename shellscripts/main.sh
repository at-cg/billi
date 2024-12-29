#!/bin/bash
#PBS -N Billi
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

cd ../src
make clean
make

DIR=../test/data/
DDIR=$DIR/gfa_files
RDIR=$DIR/results

if [ ! -d "$RDIR" ]; then
    mkdir -p $RDIR
fi

for FILE in "$DDIR"/*; do
    if [ -f "$FILE" ]; then
        FILENAME=$(basename "$FILE")
        # echo "$FILENAME"
        IFS='.' read -ra PARTS <<< "$FILENAME"
        TAG="${PARTS[0]}"
        # echo "$TAG"
        ./main $DDIR/$FILENAME $RDIR/$TAG.txt  
    fi
done
