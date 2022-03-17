#!/bin/bash
set -Eexo pipefail

# eval.annot.quality.sh
# this script takes as input a gene annotation file in gtf or gff2 format with at least exon rows and with
# at least gene_id and transcript_id in the 9th field, as well as a genome file in fasta format and produces
# statistics to assess the quality of this gene annotation. This statistics include>
# - Number of elements of each kind and several ratios (from make_summary_stat_from_annot.sh)
# - Number and % of monoexonic transcripts (from make_summary_stat_from_annot.sh)
# - Number and % of transcripts longer than several thresholds (from make_summary_stat_from_annot.sh)
# - Number and % of cdnas longer than several thresholds (from make_summary_stat_from_annot.sh)
# - Number and % of transcripts with at least 3 exons with internal exons longer than 500bp and the same for
#   all and distinct internal exons (from compute.trandinternalex.longer500bp.sh)
# - Intron canonicity and quality wrt duplicated sequences (from annot_to_introns_with_ds_and_can.sh)
# - Presence of polya sites at the end of transcripts (from find.polya.in.genome.sh followed by stats.tr.with.polya.at.end.sh)
# !!! no need to make an ok gff file since all the scripts called here do it !!!

# example
#########
# srun --mem=8G --pty bash
# cd ~/fragencode/tools/multi/Renv
# source /work/project/fragencode/tools/modules.genologin.load.2021-10-13.sh
# module load bioinfo/HOMER-v4.10.4
# pgm=~/fragencode/tools/multi/Scripts/eval.annot.quality.sh
# annot=~/fragencode/data/species/sus_scrofa/Sscrofa11.1.102/sus_scrofa.gtf
# genome=~/fragencode/data/species/sus_scrofa/Sscrofa11.1.102/sus_scrofa.fa
# outdir=~/fragencode/workspace/sdjebali/geneswitch/analysis/annotation.quality/complete
# time $pgm $annot $genome 2> $outdir/eval.annot.quality.err > $outdir/eval.annot.quality.out
# real	38m44.574s   *** if we compute the bed file of polya positions, otherwise much less

# files produced in ~/fragencode/data/species/sus_scrofa/Sscrofa11.1.102 with the different statistics
# -rw-rw-r--+ 1 sdjebali FRAGENCODE  829 Mar 11 09:12 sus_scrofa_make_summary_stat_from_annot.out
# -rw-rw-r--+ 1 sdjebali FRAGENCODE  231 Mar 11 09:18 sus_scrofa_compute.trandinternalex.longer500bp.out
# -rw-rw-r--+ 1 sdjebali FRAGENCODE 1.4K Mar 11 09:21 sus_scrofa_annot_to_introns_with_ds_and_can.out
# -rw-rw-r--+ 1 sdjebali FRAGENCODE 236M Mar 11 09:46 sus_scrofa.polyAsites.bed
# -rw-rw-r--+ 1 sdjebali FRAGENCODE  119 Mar 11 09:47 sus_scrofa.tts-50.polyA.stats.tsv


# Check inputs are provided
###########################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: eval.annot.quality.sh annot.gtf genome.fa >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a gtf or gff version 2 file containing the gene annotation to assess and including at least exon rows" >&2
    echo "  and with gene_id and transcript_id in the 9th field (anywhere)" >&2
    echo "- a fasta file for the underlying genome" >&2
    echo "" >&2
    echo "produces as output in the same directory as the annotation file:" >&2
    echo "- many statistics with regard to this annotation file (number and ratio of elements, monoexonic transcripts" >&2
    echo "  very long transcripts and cdnas, transcripts with very long internal exons, intron quality, presence" >&2
    echo "  of polya sites at the end of transcript)" >&2
    echo "" >&2
    echo "Note: since this script needs quite some time (not parallelized though), it is better to use a cluster to run it" >&2
    echo "Note: Needs R, exonerate, bedtools, homer in your path" >&2
    exit 1
fi

# Variable assignment
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
annot=$1
genome=$2
annotbasetmp=${annot%.gtf}
annotbase=${annotbasetmp%.gff}
genbasetmp=${genome%.fa}
genbase=${genbasetmp%.fasta}

# Programs
##########
MAKESUMMARY=$rootDir/make_summary_stat_from_annot.sh
INTERNEX=$rootDir/compute.trandinternalex.longer500bp.sh
INTRONS=$rootDir/annot_to_introns_with_ds_and_can.sh
POLYA=$rootDir/find.polya.in.genome.sh
STATSPA=$rootDir/stats.tr.with.polya.at.end.sh

# Start of the script
#####################

# 1. Call make_summary_stat_from_annot.sh
#########################################
echo "I am calling make_summary_stat_from_annot.sh to get different number, percent and ratios" >&2
echo "as well as number and percent of monoexonic transcripts and number and percent of transcripts" >&2
echo "and cdnas above certain lengths" >&2
$MAKESUMMARY $annot > $annotbase\_make_summary_stat_from_annot.out
echo "done" >&2
# nbex	nbdistinctex	nbtr	nbgn	nbintrons	nbdistinctintrons
# 646004	303420	63041	31908	582963	251017
# distribution of the number of exons per transcript
# # Read 63041 items
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #    1.00    3.00    7.00   10.25   14.00  151.00 
# # argmax 39568
# distribution of the number of transcripts per gene
# # Read 31908 items
# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #   1.000   1.000   1.000   1.976   2.000  10.000 
# # argmax 128
# number and percentage of monoexonic transcripts
# 63041	5452	8.64834   *** for tmerge we had 15% and for stringmerge 10% but maybe no filtering before tmerge, to check
# number and percentage of transcripts longer than 50kb, 100kb, 500kb, 1Mb, 2Mb
# 63041	19261	30.5531	10286	16.3164	611	0.969211	84	0.133247	5	0.00793135
# number and percentage of transcripts with a cdna longer than 2kb, 5kb, 10kb, 50kb, 100kb
# 63041	41419	65.7017	13999	22.2062	428	0.678923	0	0	0	0
# real	1m53.971s   

# 2. Call compute.trandinternalex.longer500bp.sh
################################################
echo "I am calling compute.trandinternalex.longer500bp.sh to get number and percent of transcripts" >&2
echo "with internal exons longer than 500bp and the same for all and distinct internal exons" >&2
$INTERNEX $annot > $annotbase\_compute.trandinternalex.longer500bp.out
echo "done" >&2
# nbtr3ex	withinternalexlonger500bp.nb	withinternalexlonger500bp.pcent	nbinternalex	longer500bp.nb1	longer500bp.pcent1	nbdistinternalex	longer500bp.nb2	longer500bp.pcent2
# 48599	13752	28.2969	525374	17282	3.28947	202258	10318	5.10141
# real	2m34.127s

# 3. Call annot_to_introns_with_ds_and_can.sh
#############################################
echo "I am calling annot_to_introns_with_ds_and_can.sh to know about intron canonicity and quality wrt duplicated sequences" >&2
$INTRONS $annot $genome > $annotbase\_annot_to_introns_with_ds_and_can.out
# long report produced
echo "done" >&2

# 4. Call find.polya.in.genome.sh
##################################
echo "I am calling find.polya.in.genome.sh to get the positions of the polya sites in the genome in case the bed file does not exist" >&2
if [ ! -f "$genbase.polyAsites.bed" ]
then
    $POLYA $genome
fi
# results in a file called $genbase.polyAsites.bed
echo "done" >&2

# 5. Call stats.tr.with.polya.at.end.sh
#######################################
echo "I am calling stats.tr.with.polya.at.end.sh to get statistics about the number and percent" >&2
echo "of transcripts with polya sites at their 3' end" >&2
$STATSPA $genbase.polyAsites.bed $annot
# results in a file called $base.tts-$flank.polyA.stats.tsv
echo "done" >&2
