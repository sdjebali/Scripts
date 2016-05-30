#!/bin/bash

# usage
# fastq_sample.sh R1.fastq.gz R2.fastq.gz pcent

if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo "Usage: fastq_sample.sh R1.fastq.gz R2.fastq.gz pcent" >&2
    echo "" >&2
    exit -1
else
    R1=$1
    R2=$2
    pcent=$3
fi

basetmp1=`basename $R1`
base1=${basetmp1%.fastq.gz}

basetmp2=`basename $R2`
base2=${basetmp2%.fastq.gz}

zcat $R1 | awk 'NR%4==1{n++; print n}' > count.txt
nb=`awk '{n++}END{print n}' count.txt`
shuf count.txt | head -n $((nb*$pcent/100)) > sample_count.txt
zcat $R1 | awk -v fileRef=sample_count.txt 'BEGIN{while (getline < fileRef >0) {ok[$1]=1}} NR%4==1{n++; if(ok[n]==1){found=1; print}} (NR%4==2||NR%4==3)&&found==1{print} NR%4==0&&found==1{print; found=0}' | gzip > subset_$pcent\_$base1.fastq.gz
zcat $R2 | awk -v fileRef=sample_count.txt 'BEGIN{while (getline < fileRef >0) {ok[$1]=1}} NR%4==1{n++; if(ok[n]==1){found=1; print}} (NR%4==2||NR%4==3)&&found==1{print} NR%4==0&&found==1{print; found=0}' | gzip > subset_$pcent\_$base2.fastq.gz
