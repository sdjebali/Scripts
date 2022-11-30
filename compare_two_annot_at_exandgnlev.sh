#!/bin/bash
set -Eexo pipefail

# compare_two_annot_at_exandgnlev.sh
# takes as input
# - two gtf files that have at least exon rows and with gene_id and transcript_id in the 9th field
# provides as output
# - the total number of distinct exons of each annot and the common nb of distinct exons between the two annot
# - the total number of genes of each annot and common nb of genes between the two annot, considering two genes as
#   common if they have at least one exonic bp common on the same strand

# !!! it needs bedtools in the path !!!
# !!! all intermediate files are indexed by process id so as to be able to run this script on the same data in parallel !!!

# example
#########
# srun --x11 --mem=8G --pty bash
# source /work/project/fragencode/tools/modules.genologin.load.2022-03-14.sh
# pgm=~/fragencode/tools/multi/Scripts/compare_two_annot_at_exandgnlev.sh
# sp=gallus_gallus
# outdir=~/fragencode/workspace/sdjebali/geneswitch/analysis/gtfs/annotation.quality/references/v102.v108/$sp
# ensv102=`ls ~/fragencode/data/species/$sp/*.102/$sp.gtf`
# ensv108=`ls ~/fragencode/data/species/$sp/*.108/$sp.gtf`
# mkdir -p $outdir
# cd $outdir   
# time $pgm $ensv102 $ensv108 2> $outdir/compare_ens_v102_v108_atexandgnlev.err > $outdir/compare_ens_v102_v108_atexandgnlev.out
# real	2m19.200s

# output $outdir/compare_ens_v102_v108_atexandgnlev.out
# nrex1   nrex2   nrexcommon.nb  nrexcommon.pcent1  nrexcommon.pcent2  gn1    gn2    gn1common.nb  gn1common.pcent  gn2common.nb  gn2common.pcent
# 222593  314417  146662         65.888             46.6457            24356  30862  20814         85.4574          19961         64.6782

# Check the arguments
#####################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: compare_two_annot_at_exandgnlev.sh annot1.gtf annot2.gtf >&2
    echo "where:" >&2
    echo "- annot1.gtf and annot2.gtf are two annotation files with at least exon rows and with gene_id and transcript_id in the 9th field" >&2
    echo "and outputs statistics about number of distinct exons in each annotation and in common and number of genes in each annot and in common" >&2
    echo "(as defined by simple stranded exonic overlap)" >&2
    echo "!!! Requires bedtools to be in your path !!!" >&2
    echo "" >&2
    exit 1
fi

# Assigns variables
###################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
an1=$1
an2=$2
an1base=${an1%.gtf}
an2base=${an2%.gtf}

# Programs
##########
MAKEOK=$rootDir/make_gff_ok.awk
GFF2GFF=$rootDir/gff2gff.awk
INTER=intersectBed

# 1. Compute the distinct exons of each annot
#############################################
echo "I am computing the number of distinct exons of each annotation" >&2
awk '$3=="exon"{print $1":"$4":"$5":"$7}' $an1 | sort | uniq > $an1base.$$.txt
awk '$3=="exon"{print $1":"$4":"$5":"$7}' $an2 | sort | uniq > $an2base.$$.txt
ex1=`wc -l $an1base.$$.txt | awk '{print $1}'`
ex2=`wc -l $an2base.$$.txt | awk '{print $1}'`
echo done >&2


# 2. Compute the distinct exons in common between the two annot
###############################################################
echo "I am computing the number of distinct exons that are common between the two annotations" >&2
excom=`cat $an1base.$$.txt $an2base.$$.txt | sort | uniq -c | awk '$1==2' | wc -l | awk '{print $1}'`
ex1pcent=`echo $ex1 $excom | awk '{print $2/$1*100}'`
ex2pcent=`echo $ex2 $excom | awk '{print $2/$1*100}'`
echo done >&2

# 3. Compute the genes of each annot
####################################
echo "I am computing the number of genes of each annot" >&2
awk '$3=="exon"' $an1 | awk -f $MAKEOK | awk -f $GFF2GFF | sort -k12,12 -k4,4n -k5,5n > $an1base.exons.ok.$$.gff
awk '$3=="exon"' $an2 | awk -f $MAKEOK | awk -f $GFF2GFF | sort -k12,12 -k4,4n -k5,5n > $an2base.exons.ok.$$.gff
gn1=`awk '{print $10}' $an1base.exons.ok.$$.gff | sort | uniq | wc -l | awk '{print $1}'`
gn2=`awk '{print $10}' $an2base.exons.ok.$$.gff | sort | uniq | wc -l | awk '{print $1}'`
echo done >&2

# 4. Compute the genes in common between the two annot
#######################################################
echo "I am computing the genes in common between the two annotations by simple stranded exonic overlap" >&2
gn1com=`$INTER -a $an1base.exons.ok.$$.gff -b $an2base.exons.ok.$$.gff -s -wao | awk '$NF>0{split($10,b,"\""); ok[b[2]]++; if(ok[b[2]]==1){n++}} END{print n}'`
gn2com=`$INTER -a $an2base.exons.ok.$$.gff -b $an1base.exons.ok.$$.gff -s -wao | awk '$NF>0{split($10,b,"\""); ok[b[2]]++; if(ok[b[2]]==1){n++}} END{print n}'`
gn1pcent=`echo $gn1 $gn1com | awk '{print $2/$1*100}'`
gn2pcent=`echo $gn2 $gn2com | awk '{print $2/$1*100}'`
echo done >&2

# 5. Output the statistics file
###############################
echo "I am making the statistics file" >&2
printf "nrex1\tnrex2\tnrexcommon.nb\tnrexcommon.pcent1\tnrexcommon.pcent2\tgn1\tgn2\tgn1common.nb\tgn1common.pcent\tgn2common.nb\tgn2common.pcent\n"
printf $ex1"\t"$ex2"\t"$excom"\t"$ex1pcent"\t"$ex2pcent"\t"$gn1"\t"$gn2"\t"$gn1com"\t"$gn1pcent"\t"$gn2com"\t"$gn2pcent"\n"
echo done >&2

# 6. Remove intermediate files
##############################
echo "I am removing intermediate files" >&2
rm $an1base.$$.txt $an2base.$$.txt
rm $an1base.exons.ok.$$.gff $an2base.exons.ok.$$.gff
echo done >&2
