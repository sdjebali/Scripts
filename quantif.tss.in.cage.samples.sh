#!/bin/bash
set -Eexo pipefail

# quantif.tss.in.cage.samples.sh
#################################
# This is to quantify a tss bed file (extended by a certain window on each side beforehand) in multiple cage samples
# for which we have bam files and associated sample (bioreplicate) identifier in a tsv file
# This script takes as input:
#############################
# - a bed6 file of tss extended as we want and with a tss id in column 4 and strand in column 6
# - a tsv file with header that has the cage library id in 1st column and the cage bam file in 2nd column
# - an optional number of threads to use (1 if nothing specified)
# This script outputs in the current working directory:
#######################################################
# - a matrix tsv file with header that has the tss as rows and the cage samples as columns and contains
#   the number of cage reads strandedly overlapping (by at least 1 bp) the tss extended that we have
# Note: this script needs samtools and bedtools to run (tested with samtools-1.10 and bedtools2-2.29.0)

# Example
#########
# cd ~/regenet/workspace/sdjebali/tss/cage
# pgm=~/fragencode/tools/multi/Scripts/quantif.tss.in.cage.samples.sh
# echo "
#    cd ~/regenet/workspace/sdjebali/tss/cage
#    module load bioinfo/samtools-1.10
#    module load bioinfo/bedtools2-2.29.0
#    $pgm gencode.v19.annotation_capped_sites_nr_ext200bpeachside_nochrM_sorted.bed lid.cage.bam.tsv 2 > quantif.tss.in.cage.samples.out 2> quantif.tss.in.cage.samples.err
# " | awk 'BEGIN{print "\#\!\/bin\/sh"} {print}' > launch.tss.cage.sh
# sbatch --mem=8G --cpus-per-task=2 -J tsscage --mail-user=sarah.djebali@inserm.fr --mail-type=END,FAIL --workdir=$PWD --export=ALL -p workq launch.tss.cage.sh

# output
# tssid	tssinfo	ENCFF928ULT	ENCFF002MHE
# chr10:100011579:100011980:+	ENSG00000230928.1:ENST00000433374.1,:antisense:antisense,:not_low	1	0
# 178815 (4 fields)  *** a bit less rows than unique $1":"$2":"$3":"$6 from initial tss file, to understand at one point


if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo "Usage: quantif.tss.in.cage.samples.sh tssext.bed cage.lid.bam.tsv <threads> > quantif.tss.in.cage.samples.out 2> quantif.tss.in.cage.samples.err" >&2
    echo "where:" >&2
    echo "- tssext.bed is a bed6+ file that has the TSS extended on each side as we want for the cage intersection" >&2
    echo "- cage.lid.bam.tsv is a tsv file with header that has at least 2 columns, the first one being the labExpId of the directional cage experiment" >&2
    echo "  and the second being the corresponding cage bam file" >&2
    echo "- threads is an optional number of threads to use for samtools (default 1)" >&2
    echo "This script produces as output in the current working directory:" >&2
    echo "- a tsv file with header with number of cage reads overlapping each tss window in each sample (1bp stranded overlap)" >&2
    echo "" >&2
    exit 1
fi

# Sets variables
################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
tss=$1
cage=$2
if [ -n "$1" ]
then
    thr=$3
else
    thr=1
fi
basetmp=`basename $tss`
base=${basetmp%%.bed}

# Programs
##########
SUBBAM=$rootDir/sam.subselect.awk


# 1. Make sorted cage bam files that only include the chromosomes that are in the tss file and excluding chrM (for sorting reasons)
###################################################################################################################################
echo "I am making sorted cage bam files that only include the chromosomes that are in the tss file and excluding chrM (for sorting reasons)" >&2
awk 'NR>=2{print $1, $2}' $cage | while read lid bam
do
   samtools view -@ $thr -h $bam | awk -v fileRef=$tss -f $SUBBAM | samtools sort -@ $thr | samtools view -@ $thr -b -o $lid.chrtss.sort.bam
done
echo "done" >&2
					
# 2. Intersect the tss windows with the reads from all the cage bam files and reformat into tss expression matrix
#################################################################################################################
echo "I am intersecting the tss windows with the reads from all the cage bam files and reformatting the results into the wanted tss cage expression matrix" >&2
namelist=`awk 'NR>=2{s=(s)($1)(" ")} END{print s}' $cage`
bamlist=`awk 'NR>=2{s=(s)($1".chrtss.sort.bam")(" ")} END{print s}' $cage`
intersectBed -a $tss -b $bamlist -names $namelist -sorted -s -C -nobuf | awk -v fileRef=$cage 'BEGIN{OFS="\t"; while (getline < fileRef >0){n++; if(n>=2){samp[n-1]=$1}}} {tss[$1":"$2":"$3":"$6]=1; info[$1":"$2":"$3":"$6]=$4; nb[$1":"$2":"$3":"$6,$7]=$8} END{OFS="\t"; for(t in tss){s=t"\t"info[t]"\t"; for(i=1; i<=(n-2); i++){s=(s)(nb[t,samp[i]])("\t")} print (s)(nb[t,samp[i]])}}' | sort -k1,1 | awk -v fileRef=$cage 'BEGIN{OFS="\t"; while (getline < fileRef >0){n++; if(n>=2){samp[n-1]=$1}} s="tssid\ttssinfo\t"; for(i=1; i<=(n-2); i++){s=(s)(samp[i])("\t")} print (s)(samp[i])} {print}' > $base.with.cage.read.count.tsv
echo "done" >&2

# 3. Remove intermediate files
##############################
# awk 'NR>=2{print $1, $2}' $cage | while read lid bam
# do
#     rm $lid.chrtss.sort.bam
# done
