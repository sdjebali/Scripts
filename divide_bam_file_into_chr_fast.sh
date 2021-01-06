#!/bin/bash
set -Eexo pipefail 

# divide_bam_file_into_chr_fast.sh
# !!! improved version of divide_bam_file_into_chr.sh !!!

# this script takes as input:
#############################
# - a bam file
# - an optional tsv file without header that has chr names we want in the first column
# and produces as output in the current working directory:
########################################################## 
# - as many bam files as there are chromosomes in the second file (or as many as present in the bam file otherwise)

# Notes
#######
# - needs samtools to be available (tested with version 1.10 on genologin slurm cluster)
# - could be improved by adding a number of threads for samtools commands

# Example of usage
##################
# cd ~/regenet/tools/abc.model
# module load bioinfo/samtools-1.10
# pgm=~/fragencode/tools/multi/Scripts/divide_bam_file_into_chr_fast.sh
# bamfile=~/regenet/results/cage/homo_sapiens/hg38/ENCFF928ULT.bam
# chrfile=~/fragencode/data/species/homo_sapiens/GRCh38/chrom.1-22XY.ucsc.len
# time $pgm $bamfile $chrfile 2> divide_bam_file_into_chr.err
# real	4m55.787s   *** instead of 7m13 for the non fast version but that produces sorted and indexed bam file
# the two versions have their merit and this one is a bit but not tremendously faster at the end

# second input file
# chr1	248956422
# 24 (2 fields)

# Check the compulsory input does exist
#######################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo "Usage: divide_bam_file_into_chr_fast.sh file.bam <chr.tsv>" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a bam file" >&2
    echo "- an optional tsv file without header that has chr names we are interested in in 1st column" >&2
    echo "" >&2
    echo "produces as output in the current working directory:" >&2
    echo "-  as many bam files as there are chromosomes in chr.tsv if provided or present in the bam file otherwise," >&2
    echo "   files that will be named after the input bam file base name" >&2
    echo "" >&2
    echo "Note: needs samtools to be available" >&2
    exit 1
fi

# Assign variables
##################
path="`dirname \"$0\"`" # relative path to exec
rootDir="`( cd \"$path\" && pwd )`" # absolute path to exec
bamfile=$1
base=`basename ${bamfile%.bam}`

# 1. Get the chromosomes present in the bam file and assign the wanted chromosome file
######################################################################################
echo "I am getting the chromosomes present in the bam file and I am assigning the wanted chromosome file variable" >&2 
samtools view -H $bamfile | awk '$1=="@SQ"{split($2,a,":"); print a[2]}' | sort -V | uniq > bamfile.chr.txt
if [ ! -n "$2" ]
then
    chrfile=bamfile.chr.txt
else
    chrfile=$2
fi
echo "done" >&2
# ERCC-00002
# 292 (1 fields)

# 2. Get the chromosomes that are both wanted and present in the bam file
##########################################################################
echo "I am getting the chromosomes that are both wanted and present in the bam file" >&2 
cat bamfile.chr.txt $chrfile | awk '{print $1}' | sort -V | uniq -c | awk '$1==2{print $2}' > common.chr.txt
echo "done" >&2 
# chr1
# 25 (1 fields)

# 3. Make the individual chromosome bam files by reading the input bam file once only
######################################################################################
echo "I am making the individual chromosome bam files by reading the input bam file once only" >&2 
samtools view -H $bamfile > header.txt
samtools view $bamfile | awk -v fileRef=common.chr.txt -v base=$base 'BEGIN{OFS="\t"; while (getline < fileRef >0){ok[$1]=1}} ok[$3]==1{print > base"."$3".sam"}'
cat common.chr.txt | while read c
do
    cat header.txt $base.$c.sam | samtools view -b > $base.$c.bam
done
echo "done" >&2

# 6. Removes unuseful files
###########################
echo "I am removing unuseful files" >&2 
rm bamfile.chr.txt
rm header.txt
cat common.chr.txt | while read c
do
    rm $base.$c.sam
done
rm common.chr.txt
echo "done" >&2 
