#!/bin/bash
set -Eexo pipefail

###################################
# genepairs2commonenhancers.sh    #
###################################
# Given:
########
# - a tsv file with header that has gene pairs as rows and gene ids in the two first columns
# - a tsv file with header that has genes as rows with their gene id in column 1, the chromosome
#   (ensembl convention) and the position of their tss in the 3rd and 4th column respectivetly
# - a tsv file with header that has genes as rows with their id in column 1 and then for each
#   cell type the semi column separated list of the enhancers (chr:beg:end, chr in ucsc convention)
#   that are contacted by the gene in each cell type (empty column when no contacted enhancers)

# Provides:
###########
# - a tsv file with header named after the 1st input that has gene pairs on the same chr and closer
#   than 4Mb with their two gene ids first, then their tss positions, then their tss to tss distance
#   (<=4Mb) and finally the cardinal of the union of their enhancers that are between 25kb and 2Mb
#   from their two tss, the cardinal of the intersection of their enhancers that are between 25kb and 2Mb
#   from their two tss and finally the corresponding jaccard percentage
# - on the standard output the format of the different intermediate and final files produced as well
#   as the R summary of the distribution of the jaccard percentage for the cis less than 4Mb gene pairs


# Example
#########
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/work/duplicons/common_enhancers
# module load statistics/R/4.3.1
# genepairs=~/work/duplicons/paper/supplementary_datasets_and_tables/supplementary_datasets7.tsv 
# genetsscoord=~/work/duplicons/paper/supplementary_datasets_and_tables/supplementary_datasets10.tsv
# contactedenhancers=~/work/duplicons/paper/supplementary_datasets_and_tables/supplementary_datasets6.tsv
# pgm=~/fragencode/tools/multi/Scripts/genepairs2commonenhancers.sh
# time $pgm $genepairs $genetsscoord $contactedenhancers > realpairs.common.out 2> realpairs.common.err
# real	0m15.867s


# 1st input file = gene pairs
#############################
# ID1	ID2	LastCommonAncestor	DuplicationAgeClass	DuplicationAge	GenomicLocalization	GenomicDistance	Ohnolog	DeltaExpression	DeltaComplexity	ExpressionDivergence	ContactDivergence
# ENSMUSG00000000056	ENSMUSG00000002280	Gnathostomata	very_old	462	trans	NA	yes	0.0880228204617022	1.181818181818184.63893233215434	4.38864335382784
# 1405 (12 fields)
# 16 (13 fields)

# 2nd input file = gene tss coordinates
#######################################
# ID	Biotype	Chr	TSSPosition	MeanTPM
# ENSMUSG00000064372	Mt_tRNA	MT	15422	NA
# 56306 (5 fields)

# 3rd input file = gene contacted enhancer coordinates
######################################################
# gene	EpiSC	ESC	ESC_18	ESC_NKO	ESC_wild	ESd_starved	ESd_TPO	FLC	preB_aged	preB_young	TSC	preadip_D0	preadip_D2	preadip_4H
# ENSMUSG00000000001	chr3:108101507:108101657;chr3:...
# 25547 (15 fields)  

# output file = supplementary_datasets7.closer4Mb.gn1.gn2.tss1.tss2.dist.contactedenhancers.union.inter.nb.jaccardpcent.tsv
###########################################################################################################################
# gnid1	gnid2	tsspos1	tsspos2	tssdist	nunion	ninter	jaccard.pcent
# ENSMUSG00000001225	ENSMUSG00000020651	12:31390871	12:31559969	169098	12	0	NA
# 330 (8 fields)

# output file = realpairs.common.out
####################################
# gnid1	gnid2	tsspos1	tsspos2	tssdist
# ENSMUSG00000001225	ENSMUSG00000020651	12:31390871	12:31559969	169098
# 330 (5 fields)
# gnid1	gnid2	tsspos1	tsspos2	tssdist	nrenhdistok1	nrenhdistok2
# ENSMUSG00000001225	ENSMUSG00000020651	12:31390871	12:31559969	169098	31700960:31701110,31703360:31703510,31701360:31701510,31121860:31122010,31721560:31721710,31719260:31719410,31720560:31720710,31720960:31721410,31719560:31719710,31927460:31927610,	31261260:31261410,31268560:31268710,
# 330 (7 fields)
# gnid1	gnid2	tsspos1	tsspos2	tssdist	nrenhdistok1	nrenhdistok2	enhunionlist	enhinterlist
# ENSMUSG00000001225	ENSMUSG00000020651	12:31390871	12:31559969	169098	31700960:31701110,31703360:31703510,31701360:31701510,31121860:31122010,31721560:31721710,31719260:31719410,31720560:31720710,31720960:31721410,31719560:31719710,31927460:31927610,	31261260:31261410,31268560:31268710,	31700960:31701110,31703360:31703510,31701360:31701510,31121860:31122010,31721560:31721710,31719260:31719410,31720560:31720710,31720960:31721410,31719560:31719710,31927460:31927610,31261260:31261410,31268560:31268710,	NA
# 330 (9 fields)
# gnid1	gnid2	tsspos1	tsspos2	tssdist	nunion	ninter	jaccard.pcent
# ENSMUSG00000001225	ENSMUSG00000020651	12:31390871	12:31559969	169098	12	0	NA
# 330 (8 fields)
# Read 329 items
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#   0.7246  12.5000  25.7915  30.7282  44.8026 100.0000       75 
# argmax 38

# Check the arguments
#####################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo Usage: genepairs2commonenhancers.sh gnpairs.tsv tsscoord.tsv contactedenh.tsv >&2
    echo "takes as input:" >&2
    echo "- gnpairs.tsv is a tsv file with header that has the gene pairs of interest as rows with their ids in columns 1 and 2" >&2
    echo "- tsscoord.tsv is a tsv file with header that has genes as rows with their id in column 1, the gene tss chromosome" >&2
    echo "  (ensembl convention) in column 3 and the gene tss position in column 4" >&2
    echo "- contactedenh.tsv is a tsv file with header that has genes as rows with their ids in column 1," >&2
    echo "  and then for each cell type the semi colon separated list of enhancers contacted as" >&2
    echo "  chr(ucsc convention):beg:end and an emtpy column if no enhancers contacted" >&2
    echo "and outputs in a file that is called after the without tsv basename of the gene pair file and in the current working directory:" >&2
    echo "- a tsv file with header that only has gene pairs that are on the same chromosome and less" >&2
    echo "  distant than 4Mb with their ids in the 2 first columns, their tss positions in the 2 next columns" >&2
    echo "  their tss to tss distance next and finally the cardinal of the union of contacted enhancers by the" >&2
    echo "  two genes that also are at a distance between 25kb and 2Mb from each of the two genes," >&2
    echo "  then the cardinal of the intersection of those and the Jaccard percentage corresponding to those" >&2 
    echo "!!! Requires R in your path !!!" >&2
    echo "" >&2
    exit 1
fi

# Assigns variables
###################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
gnpair=$1
tsscoord=$2
contactenh=$3
base=`basename ${gnpair%.tsv}`

# Programs
##########
check=$rootDir/check_simple.sh
addenh=$rootDir/add_nrenhdistoklist_to_gene_pairs.awk
unioninter=$rootDir/add_union_inter_enhlists.awk
stats=$rootDir/stats.sh


# 1. Make a tsv file with header for the gene pairs on the same chr and less than 4Mb apart,
###########################################################################################
#    that have gnid1, gnid2, tss1, tss2, dist
#############################################
echo "I am making a tsv file with header for the gene pairs on the same chr and less than 4Mb apart" >&2
awk -v fileRef=$tsscoord 'BEGIN{OFS="\t"; while (getline < fileRef >0){tsschr[$1]=$3; tsspos[$1]=$4}} NR==1{print "gnid1", "gnid2", "tsspos1", "tsspos2", "tssdist"} NR>=2{c1=tsschr[$1]; c2=tsschr[$2]; p1=tsspos[$1]; p2=tsspos[$2]; d=(c1==c2 ? abs(p1-p2) : "NA"); if(c1==c2&&d<=4000000){print $1, $2, c1":"p1, c2":"p2, d}} function abs(x){return (x>=0 ? x : (-1*x))}' $gnpair > $base.closer4Mb.gn1.gn2.tss1.tss2.dist.tsv
$check $base.closer4Mb.gn1.gn2.tss1.tss2.dist.tsv
echo done >&2

# 2. Add the list of non redundant enhancers of gene1 that are at a [25kb; 2Mb] distance from the tss of gene1 and the tss of gene2
###################################################################################################################################
#    and the same for the enhancers of gene2
############################################
echo "I am adding the list of non redundant enhancers of gene1 that are between 25kb and 2Mb from its tss and the tss of gene2" >&2
echo "and the list of non redundant enhancers of gene2 that are between 25kb and 2Mb from its tss and the tss of gene1" >&2
awk -v fileRef1=$tsscoord -v fileRef2=$contactenh -f $addenh $base.closer4Mb.gn1.gn2.tss1.tss2.dist.tsv > $base.closer4Mb.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.tsv
$check $base.closer4Mb.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.tsv
echo done >&2

# 3. Add the nr enh list made of the union of both lists and then the intersection list
#######################################################################################
echo "I am adding the enhancer list made of the union of both lists and the the intersection list" >&2
awk -f $unioninter $base.closer4Mb.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.tsv > $base.closer4Mb.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.union.inter.tsv
$check $base.closer4Mb.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.union.inter.tsv
echo done >&2

# 4. Add for each gene pair the cardinal of the union and of the inter lists as well as the jaccard percentage
##############################################################################################################
#    made by dividing the 2nd by the 1st and multiplying by 100, and remove the intermediate lists
##################################################################################################
echo "I am adding the cardinal of the union and of the inter lists as well as the corresponding Jaccard percentage" >&2
awk 'NR==1{OFS="\t"; print $1, $2, $3, $4, $5, "nunion", "ninter", "jaccard.pcent"}  NR>=2{nu=split($8,a,","); ni=split($9,b,","); if($9=="NA"){inter=0; union=nu-1; jac="NA"}else{if($8=="NA"){union=0; jac=0}else{ union=nu-1; inter=ni-1; jac=inter/union*100}} print $1, $2, $3, $4, $5, union, inter, jac}' $base.closer4Mb.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.union.inter.tsv > $base.closer4Mb.gn1.gn2.tss1.tss2.dist.contactedenhancers.union.inter.nb.jaccardpcent.tsv
$check $base.closer4Mb.gn1.gn2.tss1.tss2.dist.contactedenhancers.union.inter.nb.jaccardpcent.tsv
echo done >&2

# 5. Computes the distribution of the jaccard percentage
########################################################
echo "I am computing the distribution of the Jaccard percentage for the input read pairs closer than 4Mb" >&2
awk 'NR>=2{print $NF}' $base.closer4Mb.gn1.gn2.tss1.tss2.dist.contactedenhancers.union.inter.nb.jaccardpcent.tsv > tmp$$
$stats tmp$$
rm tmp$$
echo done >&2

# 6. Remove intermediate files
##############################
echo "I am removing intermediate files" >&2
rm $base.closer4Mb.gn1.gn2.tss1.tss2.dist.tsv
rm $base.closer4Mb.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.tsv
echo done >&2

