#!/bin/bash

# usage
# fastq_sample.sh R1.fastq.gz R2.fastq.gz pcent
# takes 4 minutes on two files of 9G each (117 876 320 read pairs), and resulting files are 300kb each

# example
#########
# cd ~/Downloads/rnaseq_course/SIB_2020
# pgm=~/Downloads/Scripts/Bash/fastq_sample.sh
# time $pgm ENCFF000EWJ.fastq.gz ENCFF000EWX.fastq.gz 1 2> ENCFF000EWJ.fastq_sample.err
# Improved on June 15th 2020 so that intermediate files are named after the first fastq.gz file and uses wc -l instead of an awk script to count rows in count file

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

zcat $R1 | awk 'NR%4==1{n++; print n}' > $base1.count.txt
nb=`wc -l $base1.count.txt | awk '{print $1}'`
shuf $base1.count.txt | head -n $((nb*$pcent/100)) > $base1.sample_count.txt
zcat $R1 | awk -v fileRef=$base1.sample_count.txt 'BEGIN{while (getline < fileRef >0) {ok[$1]=1}} NR%4==1{n++; if(ok[n]==1){found=1; print}} (NR%4==2||NR%4==3)&&found==1{print} NR%4==0&&found==1{print; found=0}' | gzip > subset_$pcent\_$base1.fastq.gz
zcat $R2 | awk -v fileRef=$base1.sample_count.txt 'BEGIN{while (getline < fileRef >0) {ok[$1]=1}} NR%4==1{n++; if(ok[n]==1){found=1; print}} (NR%4==2||NR%4==3)&&found==1{print} NR%4==0&&found==1{print; found=0}' | gzip > subset_$pcent\_$base2.fastq.gz
