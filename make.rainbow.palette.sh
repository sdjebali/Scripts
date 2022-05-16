#!/bin/bash
set -Eexo pipefail

# make.rainbow.palette.sh

# example
# cd ~/fragencode/tools/multi/Renv
# module load compiler/gcc-9.3.0
# module load system/R-4.0.4_gcc-9.3.0
# pgm=~/fragencode/tools/multi/Scripts/make.rainbow.palette.sh
# time $pgm 45 ~/fragencode/workspace/sdjebali/geneswitch/pipelines/rnaseq/debugging/may12th2022/1fcf2f0e13487ddabfd460f79886aa
# real	0m0.458s
# has produced the file ~/fragencode/workspace/sdjebali/geneswitch/pipelines/rnaseq/debugging/may12th2022/1fcf2f0e13487ddabfd460f79886aa/rainbow.45.txt 
# #FF0000
# #FF2200
# 45 (1 fields)


if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo "Usage: make.rainbow.palette.sh <integer> <outdir>" >&2
    echo "" >&2
    echo "takes as input: " >&2
    echo "- an integer N" >&2
    echo "- an output directory" >&2
    echo "" >&2
    echo "produces as output in the output directory:" >&2
    echo "- a rainbow palette called rainbow.N.txt with this number N of colors" >&2
    echo "" >&2
    exit 1
fi

echo '
p=rainbow("'$1'")
write.table(p, file = "'$2'/rainbow.'$1'.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
' | R --vanilla
