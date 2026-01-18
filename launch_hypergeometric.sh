#!/bin/bash
set -Eexo pipefail

# when we have two sets from a given universe and want to know if their intersection is more than expected by chance

# example
# we have
# - 1829 genes in 1st set
# - 1630 genes in 2nd set
# - 299 in the intersection
# - 10672 in the universe
# then we would have to ask the following to know whether it is more than expected by chance
# phyper(299-1,1829,10672+299-1829,1630,lower.tail=FALSE)
# [1] 0.02789313  *** exactly equivalent to exact fisher test one-sided

# usage:
# launch_hypergeometric.sh nA nB nA&&B nuniv
# exactly equivalent to fisher exact test one-sided (see explanations from Andrea in stat.tests.txt)

# example
# srun --x11 --mem=8G --pty bash
# cd ~/fragencode/tools/multi/Renv
# module load statistics/R/4.3.1
# launch_hypergeometric.sh 1829 1630 299 10672 
# + '[' '!' -n 1829 ']'
# + '[' '!' -n 1630 ']'
# + '[' '!' -n 299 ']'
# + '[' '!' -n 10672 ']'
# + R --vanilla --slave
# + echo 'phyper(299-1,1829,10672+299-1829,1630,lower.tail=FALSE)'
# [1] 0.02789313

if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ]
then
echo "usage: launch_hypergeometric.sh nA nB nA&&B nuniv"
else
echo 'phyper('$3'-1,'$1','$4'+'$3'-'$1','$2',lower.tail=FALSE)' | R --vanilla --slave
fi
