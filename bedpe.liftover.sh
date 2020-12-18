#!/bin/bash
set -Eexo pipefail

# bedpe.liftover.sh
###################
# This is to lift a bedpe file over from one genome assembly to another genome assembly (from the same species)
# This script takes as input:
#############################
# - a gzipped bedpe file 
# - a chain file that allows to go from the genome assembly of the bedpe file to the target genome
# - the id of the genome assembly in which the bedpe file is
# - the id of the genome target
# This script produces as output:
#################################
# - the bedpe file in the coordinates of the target genome assembly
# - some stat files indicating whether we lost a lot of elements and nucleotides in this process at different levels
# This scripts need liftOver from Kent tree to be installed

# Example
#########
# cd /work2/project/regenet/workspace/sdjebali/egprediction/from3D/predictions/capturehic/jung.ren.2019/HCmerge/score2/pall
# pgm=/work2/project/fragencode/tools/multi/Sccripts/bedpe.liftover.sh
# map=/work2/project/fragencode/data/species/homo_sapiens/liftover/hg19ToHg38.over.chain.gz
# module load bioinfo/kentUtils-v370
# /usr/bin/time -v $pgm HCmerge.pall.score2.gninfo.elt2info.bedpe.gz $map hg19 GRCh38 > bedpe.liftover.out 2> bedpe.liftover.err
# User time (seconds): 2.78
# 	System time (seconds): 0.14
# 	Percent of CPU this job got: 119%
# 	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.46
# 	Average shared text size (kbytes): 0
# 	Average unshared data size (kbytes): 0
# 	Average stack size (kbytes): 0
# 	Average total size (kbytes): 0
# 	Maximum resident set size (kbytes): 13764
# 	Average resident set size (kbytes): 0
# 	Major (requiring I/O) page faults: 0
# 	Minor (reclaiming a frame) page faults: 14994
# 	Voluntary context switches: 9702
# 	Involuntary context switches: 33
# 	Swaps: 0
# 	File system inputs: 0
# 	File system outputs: 13584
# 	Socket messages sent: 0
# 	Socket messages received: 0
# 	Signals delivered: 0
# 	Page size (bytes): 4096
# 	Exit status: 0
# Elements with smaller, bigger and equal length after liftover
# 55728	381	0.683678	286	0.513207	55061	98.8031
# Absolute difference between target and initial element length for smaller lengths and bigger lengths non equal to 0
# 55728	258807	4.64411	5977	0.107253
# Number of initial and lifted over pairs and percent of the second over the first
# 64514	64438	99.8822


# Note: was initially done for HCmerge promoter capture hic data from Jung et al, 2019 but without taking strand into account


if [ ! -n "$1" ] && [ ! -n "$2" ] && [ ! -n "$3" ] && [ ! -n "$4" ]
then
    echo "" >&2
    echo "Usage: bedpe.liftover.sh file.bedpe.gz map.file init.assembly.id target.assembly.id > stats.out 2> liftover.err" >&2
    echo "where:" >&2
    echo "- file.bedpe.gz is the gzipped bedpe file one wants to liftover from init.assembly.id to target.assembly.id" >&2
    echo "- map.file is the chain file that allows to go from init.assembly.id to target.assembly.id" >&2
    echo "This script produces as output:" >&2
    echo "- the bedpe file lifted over to target.assembly.id" >&2
    echo "- some summary statistics about the quality of the liftover at different levels on the standard output" >&2
    echo "" >&2
    exit 1
fi

path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

bedpe=$1
map=$2
assemb1=$3
assemb2=$4

base=`basename ${bedpe%.bedpe.gz}`



# 1) Make a bed file of unique elements from the bedpe input file but remembering the element coordinates
##########################################################################################################
# chr1	1183838	1209657	chr1	1491937	1496017	chr1:1183838:1209657,chr1:1491937:1496017	2.18497234006452	.	.	po	22511	4859	6	297230	2.0085	2.9873	1.0369	1.9571	1	chr1:1189288:1209265:-:ENSG00000160087.16:UBE2J2,	0	NA	elt2A	282280	NA	ENSG00000160075.10,	NA	NA	NA	NA	NA	ENSG00000160075.10,
# 64514 (33 fields)
echo "I am making a bed file of unique elements from the bedpe input file but remembering the element coordinates" >&2
zcat $bedpe | awk 'BEGIN{OFS="\t"} {str1=(($9=="+"||$9=="-") ? $9 : "."); str2=(($10=="+"||$10=="-") ? $10 : "."); print $1, $2, $3, $1":"$2":"$3":"str1, ".", str1; print $4, $5, $6, $4":"$5":"$6":"str2, ".", str2}' | sort -V -k1,1 -k2,2n -k3,3n | uniq > $base.elt1.elt2.named.uniq.bed
# chr1	1183838	1209657	chr1:1183838:1209657:.	.	.
# 55768 (6 fields)
echo "done" >&2

# 2) Lift this bed file over to the target genome assembly and make some stats about the liftover quality
#########################################################################################################
# liftOver oldFile map.chain newFile unMapped
echo "I am lifting this bed file over to the target genome assembly and make some stats about the liftover quality" >&2
liftOver $base.elt1.elt2.named.uniq.bed $map $base.elt1.elt2.named.uniq.$assemb2.bed $base.elt1.elt2.named.uniq.unmapped.to.$assemb2.bed 
# chr1	1248458	1274277	chr1:1183838:1209657:.	.	.
# 55728 (6 fields)
# how many segments with a smaller length? the same length? a bigger length?
echo "Elements with smaller, bigger and equal length after liftover"
awk '{tot++; split($4,a,":"); sz1=$3-$2; sz2=a[3]-a[2]; if(sz2<sz1){less++}else{if(sz2>sz1){more++}else{same++}}} END{OFS="\t"; print tot, less, less/tot*100, more, more/tot*100, same, same/tot*100}' $base.elt1.elt2.named.uniq.$assemb2.bed
# 55728	381	0.683678	286	0.513207	55061	98.8031
# absolute difference between lengths when non equal to 0?
echo "Absolute difference between target and initial element length for smaller lengths and bigger lengths non equal to 0"
awk '{split($4,a,":"); sz1=$3-$2; sz2=a[3]-a[2]; if(sz2<sz1){tot++; less+=(sz1-sz2)}else{tot++; if(sz2>sz1){more+=(sz2-sz1)}}} END{OFS="\t"; print tot, less, less/tot, more, more/tot}' $base.elt1.elt2.named.uniq.$assemb2.bed 
# 55728	258807	4.64411	5977	0.107253 
echo "done" >&2

# 3) Make a bedpe file from the lifted over bed file going over the initial bedpe file
######################################################################################
# chr1	1183838	1209657	chr1	1491937	1496017	chr1:1183838:1209657,chr1:1491937:1496017	2.18497234006452	.	.	po	22511	4859	6	297230	2.0085	2.9873	1.0369	1.9571	1	chr1:1189288:1209265:-:ENSG00000160087.16:UBE2J2,	0	NA	elt2A	282280	NA	ENSG00000160075.10,	NA	NA	NA	NA	NA	ENSG00000160075.10,
# 64514 (33 fields) 
echo "I am making a bedpe file from the lifted over bed file and deriving some statistics about the liftover" >&2
zcat $bedpe | awk -v fileRef=$base.elt1.elt2.named.uniq.$assemb2.bed 'BEGIN{OFS="\t"; while (getline < fileRef >0){coord[$4]=$1":"$2":"$3":"$6}} {c1=coord[$1":"$2":"$3":"$9]; c2=coord[$4":"$5":"$6":"$10]; if((c1!="")&&(c2!="")){split(c1,a,":"); split(c2,b,":"); s=$8"\t"((a[4]=="+"||a[4]=="-") ? a[4] : ".")"\t"((b[4]=="+"||b[4]=="-") ? b[4] : ".")"\t"; for(i=11; i<=NF; i++){s=(s)($i)("\t")} s=(s)($1":"$2":"$3":"$9","$4":"$5":"$6":"$10); print a[1], a[2], a[3], b[1], b[2], b[3], c1","c2, s}}' > $base.$assemb2.bedpe
# chr1	1248458	1274277	chr1	1556557	1560637	chr1:1248458:1274277:.,chr1:1556557:1560637:.	2.18497234006452	.	.	po	22511	4859	6	297230	2.0085	2.9873	1.0369	1.9571	1	chr1:1189288:1209265:-:ENSG00000160087.16:UBE2J2,	0	NA	elt2A	282280	NA	ENSG00000160075.10,	NA	NA	NA	NA	NA	ENSG00000160075.10,	chr1:1183838:1209657:.,chr1:1491937:1496017:.
# 64438 (34 fields)
# how many element pairs are lifted over and what is the proportion over the initial nb of pairs ?
echo "Number of initial and lifted over pairs and percent of the second over the first"
lifted=`wc -l $base.$assemb2.bedpe | awk '{print $1}'`
init=`zcat $bedpe | wc -l | awk '{print $1}'`
echo $init $lifted | awk 'BEGIN{OFS="\t"} {print $1, $2, $2/$1*100}' 
echo "done" >&2

