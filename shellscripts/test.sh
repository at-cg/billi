#!/bin/bash
#PBS -N Test-Small
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

# set -e

cd ../src
make clean
make
ulimit -s unlimited

# /usr/bin/time -v ./main ../test/data/gfa_files/worst_case\(100000\).gfa ../test/data/results/worst_case\(100000\)
# ./main ../test/data/gfa_files/t1-1.gfa ../test/data/results/t1-1

DIR=../test/data
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
        # if [[ "$TAG" != "EC7" ]]; then
        #     continue
        # fi
        echo "$TAG"
        ./main decompose -i $DDIR/$FILENAME -d 10 -c true -r true -p true -o $RDIR/$TAG
        # ./main decompose -i $DDIR/$FILENAME -o $RDIR/$TAG
    fi
done
