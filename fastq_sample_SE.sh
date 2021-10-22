#!/bin/bash
# set -Eexo pipefail

# usage
# fastq_sample_SE.sh file.fastq.gz pcent
# made for single end data
# for paired end data use fastq_sample.sh


# example
#########
# dir=~/fragencode/data/reads/srnaseq/sus_scrofa
# pgm=~/fragencode/tools/multi/Scripts/fastq_sample_SE.sh
# cd $dir
# time $pgm Susscrofa_FR17MAG201511203_Liver_smallRNA.fastq.gz 1 2> ESusscrofa_FR17MAG201511203_Liver_smallRNA_sample.err


if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo "Usage: fastq_sample_SE.sh file.fastq.gz pcent" >&2
    echo "" >&2
    exit -1
else
    R=$1
    pcent=$2
fi

basetmp=`basename $R`
base=${basetmp%.fastq.gz}

zcat $R | awk 'NR%4==1{n++; print n}' > $base.count.txt
nb=`wc -l $base.count.txt | awk '{print $1}'`
shuf $base.count.txt | head -n $((nb*$pcent/100)) > $base.sample_count.txt
zcat $R | awk -v fileRef=$base.sample_count.txt 'BEGIN{while (getline < fileRef >0) {ok[$1]=1}} NR%4==1{n++; if(ok[n]==1){found=1; print}} (NR%4==2||NR%4==3)&&found==1{print} NR%4==0&&found==1{print; found=0}' | gzip > subset_$pcent\_$base.fastq.gz

rm $base.count.txt
