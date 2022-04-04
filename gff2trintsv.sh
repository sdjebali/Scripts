#!/bin/bash
set -Eexo pipefail

# gff2trintsv.sh
################
# - takes as input:
###################
#   * a gtf or gff2 annotation file with at least exon rows and with gene_id and transcript_id in the 9th field (anywhere)
# - produces as output:
#######################
#   * a tsv file without header with as many rows as transcripts in the input file and with the following columns
# - trid
# - gnid
# - chr
# - strand
# - nb introns
# - in case nb introns >=1 the comma separated list of beg:end of the introns
# - in case nb introns ==1 the beg:end of the only exons
# this allows in particular to know if a tr merging procedure creates transcripts
# that are completely new compared to what it took as input

# example:
##########
# srun --mem=8G --pty bash
# annot=/work/project/fragencode/workspace/ckurylo/geneswitch/TAGADA/test_results/allbig_stringtie_vs_tmerge/annotation/novel_tmerge.gtf
# pgm=/work/project/fragencode/tools/multi/Scripts/gff2trintsv.sh
# time $pgm $annot 2> gff2trintsv.err > gff2trintsv.out 
# ENSSSCT00000029661	LOC_000000019937	X	+	5	68397940:68403937,68404130:68406503,68406595:68433861,68433945:68434496,68434633:68444177,
# ENSSSCT00000051990	LOC_000000020165	3	-	10	17093964:17094203,17094353:17094672,17094755:17095029,17095204:17096339,17096535:17096665,17096772:17096886,17096975:17097136,17097280:17097381,17097435:17097558,17097660:17098123,
# 71834 (6 fields) ***  real	2m8.279s


# Check the input file is provided
##################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo "Usage: gff2trintsv.sh annot.gff" >&2
    echo "" >&2
    echo "- where:" >&2
    echo "  * annot.gff is a gff2 file with at least exon rows and with gene_id and transcript_id in the 9th field" >&2
    echo "" >&2
    echo "- this script produces as output in the directory where the input file was, a tsv file without header" >&2
    echo "  that has as many rows as transcripts in the input file and with the following columns" >&2
    echo "  trid, gnid, chr, strand, nb introns, a string corresponding to the following:" >&2
    echo "  * in case nb introns >=1 the comma separated list of beg:end of the introns" >&2
    echo "  * in case nb introns ==1 the beg:end of the only exon" >&2
    echo "" >&2
    exit 1
fi

# Variable assigments
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
annot=$1
basetmp=${annot%.gff}
base=${basetmp%.gtf} 

# Programs
##########
MAKEOK=$rootDir/make_gff_ok.awk
GFF2GFF=$rootDir/gff2gff.awk
INTRONS=$rootDir/make_introns.awk


# Start of the process
######################
# 1. Make a gff ok file
#######################
echo "I am making a gff ok file from the input (gene_id and transcript_id as first two keys in 9th field)" >&2
awk '$3=="exon"' $annot | awk -f $MAKEOK | awk -f $GFF2GFF | sort -k12,12 -k4,4n -k5,5n > $base\_exons_sorted_by_tr.gff
echo "done" >&2

# 2. Make the intron file sorted by transcript id and then start and end for the spliced transcripts
####################################################################################################
echo "I am making the intron file sorted by transcript id and then start and end for the spliced transcripts" >&2
awk -v fldgn=10 -v fldtr=12 -f $INTRONS $base\_exons_sorted_by_tr.gff | sort -k12,12 -k4,4n -k5,5n | awk -f $GFF2GFF > $base\_introns_sorted_by_tr.gff
echo "done" >&2

# 3. Make part of the output file corresponding to spliced transcripts of the input
###################################################################################
echo "I am making the part of the wanted output file corresponding to splice transcripts" >&2
awk '{split($10,a,"\""); split($12,b,"\""); gn[b[2]]=a[2]; chr[b[2]]=$1; str[b[2]]=$7; nbint[b[2]]++; intlist[b[2]]=(intlist[b[2]])($4":"$5)(",")} END{OFS="\t"; for(t in gn){print t, gn[t], chr[t], str[t], nbint[t], intlist[t]}}' $base\_introns_sorted_by_tr.gff > $base.splicedtr.tr.gn.id.chr.str.introns.nb.list.tsv
echo "done" >&2

# 4. Make part of the output file corresponding to monoexonic transcripts of the input
######################################################################################
echo "I am making the part of the wanted output file corresponding to monoexonic transcripts" >&2
awk '{nbex[$12]++; if(nbex[$12]==1){split($10,a,"\""); gnid[$12]=a[2]; chr[$12]=$1; str[$12]=$7; coord[$12]=$4":"$5}} END{OFS="\t"; for(t in nbex){if(nbex[t]==1){split(t,a,"\""); print a[2], gnid[t], chr[t], str[t], 0, coord[t]}}}' $base\_exons_sorted_by_tr.gff > $base.monoextr.tr.gn.id.chr.str.intrnb.coord.tsv
echo "done" >&2

# 5. Make the final output file by concatenating both
#####################################################
echo "I am making the final output file by concatenating the spliced tr and the monoex tr parts" >&2
cat $base.splicedtr.tr.gn.id.chr.str.introns.nb.list.tsv $base.monoextr.tr.gn.id.chr.str.intrnb.coord.tsv > $base.alltr.tr.gn.id.chr.str.intrnb.intrlistorexcoord.tsv
echo "done" >&2

# 6. Clean intermediate files
#############################
echo "I am cleaning" >&2
rm $base\_exons_sorted_by_tr.gff
rm $base\_introns_sorted_by_tr.gff
rm $base.splicedtr.tr.gn.id.chr.str.introns.nb.list.tsv
rm $base.monoextr.tr.gn.id.chr.str.intrnb.coord.tsv
echo "done" >&2
