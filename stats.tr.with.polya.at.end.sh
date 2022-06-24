#!/bin/bash
set -Eexo pipefail

# stats.tr.with.polya.at.end.sh
# script that provides number and % of transcripts with polyA sites at their ends
# usually run after find.polya.in.genome.sh that scans a genome to obtain a bed file of polya sites
# - inputs:
#   * polya site positions in genome file (bed format), for example obtained by find.polya.in.genome.sh
#   * gene annotation (gtf format with at least exon rows and with transcript_id in the 9th field)
#   * flank integer (default 50)
# - output:
#   * report about nb and % of transcripts with a motif at its end +- flank bp
# - dependences:
#   * bedsort, bedtools

# example
#########
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/analysis/annotation.quality/tts
# module load bioinfo/bedtools2-2.29.0
# pgm=~/fragencode/tools/multi/Scripts/stats.tr.with.polya.at.end.sh
# polya=~/fragencode/data/species/gallus_gallus/GRCg6a/gallus_gallus.polyAhexamers.bed
# annot=~/fragencode/data/species/gallus_gallus/GRCg6a.102/gallus_gallus.gtf
# time $pgm $polya $annot 2> stats.tr.with.polya.at.end.err
# real	0m19.812s

# output = ~/fragencode/data/species/gallus_gallus/GRCg6a.102/gallus_gallus.tts-50.polyA.stats.tsv
# /home/sdjebali/fragencode/data/species/gallus_gallus/GRCg6a.102/gallus_gallus.gtf	16368	39288	41.66%	14485	33720	42.96%


# Check inputs are provided
###########################
if [ ! -n "$1" ] || [ ! -n "$2" ] 
then
    echo "" >&2
    echo "Usage: stats.tr.with.polya.at.end.sh polya.bed annot.gtf <flank>" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a bed file of polyA site positions in a genome" >&2
    echo "- an annotation file in gtf format (with at least exon rows and with transcript_id in the 9th field)" >&2
    echo "- an optional integer representing the flank (# bp from transcripts'ends where polyA sites will be looked for)" >&2
    echo "" >&2
    echo "produces as output in the same directory as the gff file:" >&2
    echo "- a tsv file with the number and % of transcripts with polyA sites at their ends (+-flank) located at the same place as the gtf file" >&2
    echo "" >&2
    echo "Note: Needs bedtools in your path (tested with v2-2.29.0 on the genotoul slurm cluster)" >&2
    exit 1
fi

# Variable assignment
#####################
polya=$1
annot=$2
re='^[0-9]+$'
# if a third parameter is provided and it is a number then we assign it to flank
if [ -n "$3" ] && [[ $3 =~ $re ]]
then
    flank=$3
else
    flank=50
fi
basetmp=${annot%.gtf}
base=${basetmp%.gff}
all=$base.tts-$flank.all.bed
distinct=$base.tts-$flank.distinct.bed
out=$base.tts-$flank.polyA.stats.tsv

# Program assignment
####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
percent=$rootDir/percentage.awk

# Start of the script
#####################
# 1. Makes the sorted bed file of annotated transcript ends (tts) flanked by flank bp on each side 
##################################################################################################
echo "I am making the sorted bed file of annotated transcript ends (tts) flanked by flank bp on each side" >&2
sed 's/"/ /g' $annot | awk '!/^#/ && $3=="exon"{for(i=7;i<NF;i++){if ($i=="transcript_id"){id=$(i+1); chr[id]=$1; str[id]=$7; if((gbeg[id]=="")||($4<gbeg[id])){gbeg[id]=$4} if((gend[id]=="")||($5>gend[id])){gend[id]=$5}}}} END{for(t in chr){print chr[t], gbeg[t]-1, gend[t], t, ".",str[t]}}' OFS="\t" | sort -k1,1 -k2,2n -k3,3n | uniq | awk -v flank=$flank '$6=="+"{$2=$3-flank}$6=="-"{$3=$2+flank}1' OFS="\t" | sort -k1,1 -k2,2n -k3,3n > $all
# sed 's/"/ /g' $annot | awk '!/^#/ && $3=="transcript"{id=".";for(i=7;i<NF;i++){if ($i=="transcript_id"){id=$(i+1)}} print $1,$4-1,$5,id,".",$7}' OFS="\t" | sort -k1,1 -k2,2n -k3,3n | uniq | awk -v flank=$flank '$6=="+"{$2=$3-flank}$6=="-"{$3=$2+flank}1' OFS="\t" | sort -k1,1 -k2,2n -k3,3n > $all
echo "done" >&2

# 2. Makes the sorted bed file of distinct annotated transcript ends (tts) flanked by flank bp on each side 
###########################################################################################################
echo "I am making the sorted bed file of distinct annotated transcript ends (tts) flanked by flank bp on each side" >&2
cat $all | awk '{$4="."}1' OFS="\t" | sort -k1,1 -k2,2n -k3,3n | uniq > $distinct
echo "done" >&2

# 3. Computes the number of total and distinct annotated transcript ends (tts) flanked by flank bp on each side
###############################################################################################################
#    as well as the subset of them that overlap a polya site
############################################################
echo "I am computing the number of total and distinct annotated transcript ends (tss) flanked by flank bp on each side" >&2
echo "as well as the subet of them that overlap a polyA site" >&2
Na=`cat $all | wc -l | awk '{print $1}'`
Nd=`cat $distinct | wc -l | awk '{print $1}'`
na=`bedtools intersect -u -s -a $all -b $polya | wc -l | awk '{print $1}'`
nd=`bedtools intersect -u -s -a $distinct -b $polya | wc -l | awk '{print $1}'`
if (( "$Na" > 0 ))
then
    pa=`echo $na $Na | $percent`
else
    pa="NA"
fi
if (( "$Nd" > 0 ))
then
    pd=`echo $nd $Nd | $percent`
else
    pd="NA"
fi
echo "done" >&2

# 4. Puts everything into a statistic file in tsv format
########################################################
echo "I am putting all the computed numbers into a statistic file in tsv format" >&2
echo $annot $pa $pd | sed 's/ /\t/g' > $out
echo "done" >&2

# 5. Cleans
###########
echo "I am cleaning" >&2
rm $all
rm $distinct
echo "done" >&2
