#!/bin/bash
set -Eexo pipefail 

# divide_bam_file_into_chr.sh
# !!! there is a new script called divide_bam_file_into_chr_fast.sh that takes at least twice less time !!!

# this script takes as input:
#############################
# - a bam file
# - an optional tsv file without header that has chr names we want in the first column
# and produces as output
########################
# - in the same directory as the bam file a sorted bam file and its bai index by samtools
# - in the current working directory:
#   * as many bam files as there are chromosomes in the second file (or as many as present in the bam file otherwise)

# Notes
#######
# - needs samtools to be available (tested with version 1.10 on genologin slurm cluster)
# - could be improved by adding a number of threads for samtools commands

# Example of usage
##################
# cd ~/regenet/tools/abc.model
# module load bioinfo/samtools-1.10
# pgm=~/fragencode/tools/multi/Scripts/divide_bam_file_into_chr.sh
# bamfile=~/regenet/results/cage/homo_sapiens/hg38/ENCFF928ULT.bam
# chrfile=~/fragencode/data/species/homo_sapiens/GRCh38/chrom.1-22XY.ucsc.len
# time $pgm $bamfile $chrfile 2> divide_bam_file_into_chr.err
# real	7m13.224s 

# second input file
# chr1	248956422
# chr2	242193529
# 24 (2 fields)

# Check the compulsory input does exist
#######################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo "Usage: divide_bam_file_into_chr.sh file.bam <chr.tsv>" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a bam file" >&2
    echo "- an optional tsv file without header that has chr names we are interested in in 1st column" >&2
    echo "" >&2
    echo "produces as output:" >&2
    echo "- in the directory of the bam file a sorted version of it and its samtools bai index file" >&2
    echo "- in the current working directory as many bam files as there are chromosomes in chr.tsv if provided," >&2
    echo "  or present in the bam file otherwise, files that will be named after the input bam file base name" >&2
    echo "" >&2
    echo "Note: needs samtools to be available" >&2
    exit 1
fi

# Assign variables
##################
path="`dirname \"$0\"`" # relative path to exec
rootDir="`( cd \"$path\" && pwd )`" # absolute path to exec
bamfile=$1
basetmp=${bamfile%.bam}
base=`basename $basetmp`

# 1. Sort the bam file
######################
echo "I am sorting the bam file" >&2 
samtools sort $bamfile > $base.sorted.bam
echo "done" >&2 

# 2. Index the sorted bam file
##############################
echo "I am sorting the bam file" >&2 
samtools index $base.sorted.bam > $base.sorted.bam.bai
echo "done" >&2 

# 3. Get the chromosomes present in the bam file and assign the wanted chromosome file
######################################################################################
echo "I am getting the chromosomes present in the bam file and assigning the wanted chromosome file" >&2 
samtools view -H $base.sorted.bam | awk '$1=="@SQ"{split($2,a,":"); print a[2]}' | sort -V | uniq > bamfile.chr.txt
if [ ! -n "$2" ]
then
    chrfile=bamfile.chr.txt
else
    chrfile=$2
fi
echo "done" >&2
# ERCC-00002
# 292 (1 fields)

# 4. Get the chromosomes that are both wanted and present in the bam file
##########################################################################
echo "I am getting the chromosomes that are both wanted and present in the bam file" >&2 
cat bamfile.chr.txt $chrfile | awk '{print $1}' | sort -V | uniq -c | awk '$1==2{print $2}' > common.chr.txt
echo "done" >&2 
# chr1
# 25 (1 fields)

# 5. Make the individual chromosome bam files
#############################################
echo "I am making the individual chromosome bam files" >&2 
cat common.chr.txt | while read c
do
    samtools view $base.sorted.bam $c -b > $base.sorted.$c.bam
done
echo "done" >&2
# apparently leaves the complete header in each chr indiv bam file

# 6. Removes unuseful files
###########################
echo "I am removing unuseful files" >&2 
rm bamfile.chr.txt
rm common.chr.txt
rm $base.sorted.bam $base.sorted.bam.bai
echo "done" >&2 
