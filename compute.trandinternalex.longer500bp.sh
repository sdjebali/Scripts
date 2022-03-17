#!/bin/bash
set -Eexo pipefail

# compute.trandinternalex.longer500bp.sh
# this script takes as inpute a gff2/gtf file with at least exon rows and that has gene_id and transcript_id
# in the 9th field, and compute the number and % of transcripts with at least 3 exons that have at least
# one internal exon longer than 500 bp as well as the number and % of total and distinct internal exons
# longer than 500 bp


# example of usage
##################
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/analysis/annotation.quality/internalex500
# pgm=~/fragencode/tools/multi/Scripts/compute.trandinternalex.longer500bp.sh
# annot=~/fragencode/data/species/sus_scrofa/Sscrofa11.1.102/sus_scrofa.gtf
# time $pgm $annot 2> compute.trandinternalex.longer500bp.err > compute.trandinternalex.longer500bp.out
# nbtr3ex	withinternalexlonger500bp.nb	withinternalexlonger500bp.pcent	nbinternalex	longer500bp.nb1	longer500bp.pcent1	nbdistinternalex	longer500bp.nb2	longer500bp.pcent2
# 48599	13752	28.2969	525374	17282	3.28947	202258	10318	5.10141
# real	2m34.127s

# Check input file is provided otherwise exit with error
########################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: compute.trandinternalex.longer500bp.sh annot.gff >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a gene annotation file in gtf or gff version 2 format that has at lease exon rows and" >&2
    echo "  with gene_id and transcript_id in the 9th field (anywhere)" >&2
    echo "" >&2
    echo "produces as output:" >&2
    echo "- a tsv file with header that has the number and pcent of transcripts with at least" >&2
    echo "  3 exons that have at least one internal exon longer than 500 bp as well as the number" >&2
    echo "  and percent of total and distinct internal exons that are longer than 500 bp" >&2
    echo "" >&2
    exit 1
fi

# Variable assignment
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
annot=$1
basetmp=`basename ${annot%.gtf}`
base=${basetmp%.gff}

# Programs
##########
MAKEOK=$rootDir/make_gff_ok.awk
GFF2GFF=$rootDir/gff2gff.awk
COMPUTE=$rootDir/compute.internalexlonger500bp.awk



# 1. Make a proper exon gff file with gene_id and transcript_id as the two first keys in the 9th field
######################################################################################################
#    and sorted by transcript id and then exon coord
####################################################
echo "I am making an exon gff file with gene_id and transcript_id as the two first keys in the 9th field and sorted by tr and coord" >&2
awk '$3=="exon"' $annot | awk -f $MAKEOK | awk -f $GFF2GFF | sort -k12,12 -k4,4n -k5,5n > $base\_exons_sorted_by_tr.gff
echo "done" >&2

# 2. Make a file only including transcripts with at least 3 exons
##################################################################
echo "I am making a file only including transcripts with at least 3 exons" >&2
awk '{nbex[$12]++; row[$12,nbex[$12]]=$0} END{for(t in nbex){if(nbex[t]>=3){for(i=1; i<=nbex[t]; i++){print row[t,i]}}}}' $base\_exons_sorted_by_tr.gff | awk -f $MAKEOK | sort -k12,12 -k4,4n -k5,5n > $base\_trmorethan3ex_exons_sorted_by_tr.gff
echo "done" >&2

# 3. Make a file with internal exons of transcripts with at least 3 exons
#########################################################################
echo "I am making a file with internal exons of transcripts with at least 3 exons" >&2
awk '{nbex[$12]++; row[$12,nbex[$12]]=$0} END{for(t in nbex){for(i=2; i<=(nbex[t]-1); i++){print row[t,i]}}}' $base\_trmorethan3ex_exons_sorted_by_tr.gff | awk -f $MAKEOK | sort -k12,12 -k4,4n -k5,5n > $base\_trmorethan3ex_internalexons_sorted_by_tr.gff
echo "done" >&2

# 4. Compute the number and % of transcripts with at least 3 exons that have internal exons more than 500bp
###########################################################################################################
#    as well as the number and % of total and distinct internal exons longer than 500bp
########################################################################################
echo "I am computing the number and percent of transcripts with at least 3 exons that have internal exons longer than 500bp" >&2
echo "as well as the number and percent of all and distinct internal exons longer than 500bp" >&2
awk -f $COMPUTE $base\_trmorethan3ex_internalexons_sorted_by_tr.gff
echo "done" >&2

# 5. Delete intermediate files
##############################
echo "I am removing intermediate files" >&2
rm $base\_exons_sorted_by_tr.gff
rm $base\_trmorethan3ex_exons_sorted_by_tr.gff
rm $base\_trmorethan3ex_internalexons_sorted_by_tr.gff
echo "done" >&2
