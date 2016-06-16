#!/bin/bash
# /work2/project/fragencode/tools/peaks_to_distrib.sh

# usage
# peaks_to_distrib.sh peaks.bed annot.tsv

# example
# cd /work2/project/fragencode/workspace/sdjebali/atacseq/pig/liver/ATAC14_ACCACT_L002.q10.rmdup.rmmt.macs_peaks.narrowPeak
# peaks=~kmunyard/fragencode/results/atacseq/sus_scrofa/liver/peaks/ATAC14_ACCACT_L002.q10.rmdup.rmmt.macs_peaks.narrowPeak
# annot=/work2/project/fragencode/workspace/sdjebali/atacseq/pig/annotelt_file.tsv
# /work2/project/fragencode/tools/peaks_to_distrib.sh $peaks $annot

if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo "Usage: peaks_to_distrib.sh peaks.bed annot.tsv > dist_stats.tsv" >&2
    echo "where:" >&2
    echo "- peaks.bed is the list of peaks in bed format" >&2
    echo "- annot.tsv is a 2 column tsv file with no header with the name of the annotated element and the file containing these elements" >&2
    echo "  this should include at least the following categories: prom1Kb, prom5Kb, utr, cds and intron" >&2
    echo "  for example for pig: /work2/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.80/annotelt_file.tsv" >&2
    echo "- dist_stats.tsv is a tsv file containing the distribution of peaks into the annotation" >&2
    echo "BE CAREFUL: do not run several times in the same directory since it uses fixed names for outputs" >&2
    echo "" >&2
    exit 1
fi

# To improve
# do not use fixed names for output

peaks=$1
meta=$2

# Intersect the peaks with each kind of element from the annotation
###################################################################
# assuming the following categories
###################################
# cds     /work2/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.80/sus_scrofa.cds.positions.bed
# exon    /work2/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.80/sus_scrofa.exon.positions.bed
# gene    /work2/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.80/sus_scrofa.gene.positions.bed
# intron  /work2/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.80/sus_scrofa.intron.positions.bed
# prom1Kb /work2/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.80/sus_scrofa.prom1Kb.positions.bed
# prom5Kb /work2/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.80/sus_scrofa.prom5Kb.positions.bed
# utr     /work2/project/fragencode/data/species/sus_scrofa/Sscrofa10.2.80/sus_scrofa.utr.positions.bed
cat $meta | while read elt f
do
intersectBed -a $peaks -b $f -u > peaks_over_$elt.bed
done

# Make the distribution of peaks with 2 different promoter definitions (1Kb and 5KB) and using this order in case of conflict
#############################################################################################################################
# - prom
# - utr
# - cds
# - intron
# - ig
# which means the following rules:
#################################
# - prom = over prom
# - utr = over utr but not in prom
# - cds = over cds but not in prom
# - intron = over intron but not in utr or cds or prom
# - ig = rest
tot=`wc -l $peaks | awk '{print $1}'`
for d in 1Kb 5Kb
do
prom=`wc -l peaks_over_prom$d.bed | awk '{print $1}'`
utr=`awk -v fileRef=peaks_over_prom$d.bed 'BEGIN{while (getline < fileRef >0){ko[$0]=1}} ko[$0]!=1' peaks_over_utr.bed | wc -l | awk '{print $1}'`
cds=`awk -v fileRef=peaks_over_prom$d.bed 'BEGIN{while (getline < fileRef >0){ko[$0]=1}} ko[$0]!=1' peaks_over_cds.bed | wc -l | awk '{print $1}'`
intron=`awk -v fileRef1=peaks_over_prom$d.bed -v fileRef2=peaks_over_utr.bed -v fileRef3=peaks_over_cds.bed 'BEGIN{while (getline < fileRef1 >0){ko[$0]=1}while (getline < fileRef2 >0){ko[$0]=1}while (getline < fileRef3 >0){ko[$0]=1}} ko[$0]!=1' peaks_over_intron.bed | wc -l | awk '{print $1}'`
ig=$((tot-prom-utr-cds-intron))
echo $prom $utr $cds $intron $ig | awk -v d=$d -v tot=$tot 'BEGIN{OFS="\t"}{print d, $1, $1/tot*100, $2, $2/tot*100, $3, $3/tot*100, $4, $4/tot*100, $5, $5/tot*100}'
done | awk 'BEGIN{OFS="\t"; print "dist", "prom_nb", "prom_pcent", "utr_nb", "utr_pcent", "cds_nb", "cds_pcent", "intron_nb", "intron_pcent", "ig_nb", "ig_pcent"}{print}' 
# dist    prom_nb prom_pcent      utr_nb  utr_pcent       cds_nb  cds_pcent       intron_nb       intron_pcent    ig_nb   ig_pcent
# 1Kb     696     36.5546 171     8.98109 123     6.46008 106     5.56723 808     42.437
# 5Kb     788     41.3866 161     8.45588 113     5.93487 100     5.2521  742     38.9706
