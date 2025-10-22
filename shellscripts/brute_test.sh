#!/bin/bash
#PBS -N Brute
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o op.log_brute
#PBS -e op.err_brute
 
if [ -n "$PBS_JOBID" ] || [ -n "$PBS_O_WORKDIR" ]; then
    cd "$PBS_O_WORKDIR"
else
    cd "$(dirname "$0")" || exit 1
fi

# set -e

cd ../src
make -f makefile_brute clean
make -f makefile_brute
ulimit -s unlimited

# /usr/bin/time -v ./main_brute ../test/data/gfa_files/worst_case\(100000\).gfa ../test/data/results/worst_case\(100000\)
# ./main_brute ../test/data/gfa_files/EC5.gfa ../test/data/results/EC5


DIR=../test/data/
DDIR=$DIR/gfa_files
RDIR=$DIR/results_brute

if [ ! -d "$RDIR" ]; then
    mkdir -p $RDIR
fi

for FILE in "$DDIR"/*; do
    if [ -f "$FILE" ]; then
        FILENAME=$(basename "$FILE")
        # echo "$FILENAME"
        IFS='.' read -ra PARTS <<< "$FILENAME"
        TAG="${PARTS[0]}"
        # if [[ "$TAG" != "EC15" ]]; then
        #     continue
        # fi
    
        ./main_brute decompose -i $DDIR/$FILENAME -o $RDIR/$TAG -c

        if ! diff -q "$DIR/results/$TAG/panbubble.txt" "$DIR/results_brute/$TAG/panbubble.txt" > /dev/null; then
            echo "$TAG"
        fi

    fi
done
