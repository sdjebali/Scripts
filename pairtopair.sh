#!/bin/bash

# Make the script fail if there is anything wrong in any pipe
set -Eexo pipefail

##########################
# script pairtopair.sh   #
##########################

# This script takes as input two bedpe files and outputs a bedpe file corresponding to a subset of the 1st bedpe file
# provided as input. This subset corresponds to relations of the 1st file that overlap at least one relation from the
# 2nd file on each side and in the same order, meaning those rows from the 1st file where the 1st segment overlaps
# the 1st segment of a relation from the 2nd file and where the 2nd segment overlaps the 2nd segment of the same relation
# (note that bedtools pairtopair is supposed to do that by default but does not seem to always work)

# Notes:
# - The input files do not have to be sorted
# - Only work on the cis (same chr) relations of the two files for the moment, unless I am doing it differently


# Example
# srun --x11 --mem=16G --time=10:00:00 --pty bash
# cd ~/work/contact_evolution/hs.pchic.compare
# module load bioinfo/bedtools/2.31.1
# pgm=~/fragencode/tools/multi/Scripts/pairtopair.sh
# file1=/work/project/bridge/workspace/sdjebali/cne_regulatory_evolution/results/hs.pchic.compare/jung.pchic.sc2.allrelations.hg38.bedpe
# file2=/work/project/bridge/workspace/sdjebali/cne_regulatory_evolution/results/hs.pchic.compare/laverre.allrelations.bedpe
# time $pgm $file1 $file2 > pairtopair.out 2> pairtopair.err

# input $file1 is like this
# chr1	980140	996796	chr1	781054	790044
# 664418 (6 fields)
# input $file2 is like this
# chr1	912789	915238	chr1	938702	940711
# 910180 (6 fields)
# output jung.pchic.sc2.allrelations.hg38.over.2ndfile.bedpe is like this
# chr1	980140	996796	chr1	926972	939708	chr1_968261_992014:chr1_938702_940711,chr1_968261_992014:chr1_915239_938701,
# 188146 (7 fields)


# Check the inputs exist
########################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: pairtopair.sh file1.bedpe file2.bedpe >&2
    echo "" >&2
    echo "takes as input two bedpe files in which only the 6 first fields will be considered" >&2
    echo "and produces as output a bedpe file that is a subset of the first input bedpe file" >&2
    echo "that includes relations that are overlapped by relations of the 2nd file on each side" >&2
    echo "and with segments of the relation ordered the same way (unlike bedtools pairtopair)" >&2
    echo "Note: The input files do not have to be sorted" >&2
    echo "Note: requires bedtools to be installed" >&2
    exit 1
else
    file1=$1
    file2=$2
fi

# Variable assignment
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
base=`basename ${file1%.bedpe}`

# 1. Make unique sorted bed files of 1st and 2nd segments from the 1st and 2nd files respectively
#################################################################################################
echo "I am making unique sorted bed files of 1st and 2nd segments from 1st and 2nd bedpe input files respectively" >&2
time for f in $file1 $file2
do
    awk -v f=$f 'BEGIN{OFS="\t"; split(f,a,".bedpe")} {seen1[$1"_"$2"_"$3]++; if(seen1[$1"_"$2"_"$3]==1){print $1, $2, $3 > a[1]".1stseg.bed"} seen2[$4"_"$5"_"$6]++; if(seen2[$4"_"$5"_"$6]==1){print $4, $5, $6 > a[1]".2ndseg.bed"}}' $f
    for elt in 1stseg 2ndseg
    do
	sort -V -k1,1 -k2,2n -k3,3n ${f%.bedpe}.$elt.bed > ${f%.bedpe}.$elt.sorted.bed
    done
done
echo "done" >&2
# chr1	980140	996796
# 16774 (3 fields)
# chr1	629952	634968
# 282715 (3 fields)
# chr1	912789	915238
# 19389 (3 fields)
# chr1	83360	84049
# 308359 (3 fields) *** real	0m12.354s

# 2. Intersect the 1st segments of the two input files and the 2nd segments of the two input files
##################################################################################################
echo "I am intersecting the 1st segments together and the 2nd segments together" >&2
for elt in 1stseg 2ndseg
do
    intersectBed -a ${file1%.bedpe}.$elt.sorted.bed -b ${file2%.bedpe}.$elt.sorted.bed -wao | awk '$NF!=0{overlist[$1"_"$2"_"$3]=(overlist[$1"_"$2"_"$3])($4"_"$5"_"$6)(",")} END{OFS="\t"; for(e in overlist){split(e,a,"_"); print a[1], a[2], a[3], overlist[e]}}' > ${file1%.bedpe}.$elt.sorted.over.2ndfile.$elt.bed
done
echo "done" >&2
# chr15	39770954	39786474	chr15_39778011_39785052,
# 13667 (4 fields)
# chrX	40568833	40574572	chrX_40566921_40570741,chrX_40570742_40571950,
# 182900 (4 fields)

# 3. Report in the 1st input file and just for the relations with intersections on both sides (segments)
########################################################################################################
#    the list of 1st and the list of 2nd segments from the 2nd input file that overlap with the 1st and 
#######################################################################################################
#    2nd segments of the 1st file respectively
##############################################
echo "For the relations of the 1st input file that have overlap for both segments with segments from the 2nd input file" >&2
echo "I am reporting the list of the 2nd input file 1st segments and the list of the 2nd input file 2nd segments" >&2
echo "that overlap with them" >&2
awk -v fileRef1=${file1%.bedpe}.1stseg.sorted.over.2ndfile.1stseg.bed -v fileRef2=${file1%.bedpe}.2ndseg.sorted.over.2ndfile.2ndseg.bed 'BEGIN{OFS="\t"; while (getline < fileRef1 >0){overlist1[$1"_"$2"_"$3]=$4} while (getline < fileRef2 >0){overlist2[$1"_"$2"_"$3]=$4}} overlist1[$1"_"$2"_"$3]!=""&&overlist2[$4"_"$5"_"$6]!=""{print $0, overlist1[$1"_"$2"_"$3], overlist2[$4"_"$5"_"$6]}' $file1 > $base.withoverlist.1stseg.2ndseg.2ndfile.bedpe
echo "done" >&2
# chr1	980140	996796	chr1	926972	939708	chr1_968261_992014,	chr1_938702_940711,chr1_915239_938701,
# 417686 (8 fields) *** about 2/3 of the rows here

# 4. For each such 1st input file relation, look in the 2nd input file whether there is one relation
####################################################################################################
#    made of a segment from the 1st list and a segment from the 2nd list and report those rows as output
########################################################################################################
#    with the information of the relations from the 2nd input file that harbour such segments
##############################################################################################
echo "For each such 1st input file relation, look in the 2nd input file whethere there is one relation" >&2
echo "made of a segment from the 1st list and a segment from the 2nd list and report those rows as output" >&2
echo "with the information of the relations from the 2nd input file that harbour such segments" >&2
awk -v fileRef=$file2 'BEGIN{OFS="\t"; while (getline < fileRef >0){ok2[$1"_"$2"_"$3":"$4"_"$5"_"$6]=1}} {s=""; split($7,a,","); split($8,b,","); k=1; while(a[k]!=""){l=1; while(b[l]!=""){if(ok2[a[k]":"b[l]]==1){s=(s)(a[k]":"b[l])(",")} l++} k++} if(s!=""){print $1, $2, $3, $4, $5, $6, s}}' $base.withoverlist.1stseg.2ndseg.2ndfile.bedpe > $base.over.2ndfile.bedpe
echo "done" >&2
# chr1	980140	996796	chr1	926972	939708	chr1_968261_992014:chr1_938702_940711,chr1_968261_992014:chr1_915239_938701,
# 188146 (7 fields)  *** -1/3 of the input and seems fine

# 5. Clean
##########
echo "I am cleaning" >&2
for f in $file1 $file2
do
    for elt in 1stseg 2ndseg
    do
	rm ${f%.bedpe}.$elt.sorted.bed
    done
done
rm ${file1%.bedpe}.1stseg.sorted.over.2ndfile.1stseg.bed ${file1%.bedpe}.2ndseg.sorted.over.2ndfile.2ndseg.bed 
rm $base.withoverlist.1stseg.2ndseg.2ndfile.bedpe
echo "done" >&2
