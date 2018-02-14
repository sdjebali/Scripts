#!/bin/bash
# refine_comptr_output.sh
# on sept 22nd 2015 make it possible to pass input file with comment
# on dec 15th 2015 make it possible to have annot and transcript file with gene_id and transcript_id anywhere and to have only exons in predictions
# on Feb 3rd 2017 make it possible to use an annotation file which name ends in gff 
# on jan 12th 2018 add a column with intermediate class which is one of the following 4
#   1) annot = predicted transcript that have exactly the same exon structure as a reference tr but that do not
#    extend the annotated transcript on either side (I can use my Exact class but add the constraint that there
#    is no extension)
#   2) extension = predicted transcript that have exactly the same exon structure as a reference tr but that are
#    not in 1) but (I can use my Exact class but that are not in class 1)
#   3) novel_of_annot = new isoform or variant of reference genes = predicted transcripts that are not in 1) or 2) but that
#    have at least one common intron with a reference annotated tr (I have to compute it from scratch by
#    asking for one common intron same strand as the reference and then ask that the transcript is neither
#    in class 1 nor in class 2                                                 
#   4) novel = new transcript = the ones not in 1) or 2) or 3)

# takes as input two mandatory arguments:
########################################
# - a set of transcripts to compare to the annotation as a gff 2 file containing at least exon rows
# - an annotation gtf file that has at least exon and gene rows

# and provides as output:
#########################
# - the comptr tsv output (classes are Monoexonic, Overlap, Inclusion, Extension, Intergenic_or_antisense, Exact) 
#   with additional and broader classes (annot, antisense, extension, intergenic), with the number of exons and
#   with an intermediate class that is annotated, extension, new transcript of annotated gene and new transcript

# Example of usage
##################
# cd ~/ENCODE_AWG/Analyses/Human_Promo_cells/Cufflinks_models/All/Classif/TryScript
# mytr=~/ENCODE_AWG/Analyses/Human_Promo_cells/Cufflinks_models/All/37exp_cuffafterfluxidr0.1_exons_concat.gff_merged_transcripts_stranded_over_cage_ok.gff
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/Long/gencode.v19.annotation.long.gtf
# time refine_comptr_output.sh $mytr $annot 2> refine_comptr_output.err
# real	1m35.986s

# Another example of usage
##########################
# cd /work/project/fragencode/results/rnaseq/bos_taurus/assembled
# pred=/work/project/fragencode/workspace/sdjebali/fragencode/rnaseq/analysis/transcript_models/bos_taurus/new/tpm0.1_2sample/bos_taurus_cuff_tpm0.1_2sample.gtf
# annot=/work/project/fragencode/data/species/bos_taurus/UMD3.1.90/bos_taurus_ok.gff
# time /work/project/fragencode/tools/multi/Scripts/Bash/refine_comptr_output.sh $pred $annot > refine_comptr_output.out 2> refine_comptr_output.err

# Important note: 
#################
# it uses several programs that need to be present on the system
# - comptr
# - overlap
# - make_tss...

# Be careful: 
#############
# cannot be run in parallel in the same directory

# Improvement:
##############
# also consider a distance for the transcripts on the same strand but close
# and consider those as extension and not ig as it is now


# Check the arguments
#####################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: refine_comptr_output.sh mytr.gff annot.gtf >&2
    echo "where:" >&2
    echo "- mytr.gff is a set of transcripts in gff 2 format (exon rows at least) one wants to compare to the annotation" >&2
    echo "- annot.gtf is an annotation file in gtf format (exon and gene rows at least) one wants to compare to the annotation" >&2
    echo "!!! Cannot be run in parallel in the same directory since it produces files with constant name !!!" >&2
    echo "" >&2
    exit 1
fi

path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
mytr=$1
annot=$2
mytrbasetmp=`basename ${mytr%.gff}`
mytrbase=`basename ${mytrbasetmp%.gtf}`
annbasetmp=`basename ${annot%.gtf}`
annbase=${annbasetmp%.gff}

# Programs
##########
MAKEOK=$rootDir/../Awk/make_gff_ok.awk
BOUNDARIES=$rootDir/../Awk/compute_boundaries.awk
MAKETSS=$rootDir/make_TSS_file_from_annotation_simple.sh
MAKESUM=$rootDir/make_summary_stat_from_annot.sh
COMPTR=$rootDir/../bin/comptr
OVERLAP=$rootDir/../bin/overlap
GFF2GFF=$rootDir/../Awk/gff2gff.awk
MAKEINTRONS=$rootDir/../Awk/make_introns.awk 
ADDCLASS=$rootDir/../Awk/add_intermclass_andgnlist.awk
INTER=intersectBed 

# 1. Make gff files of annotated exons, genes and TSS respectively
###################################################################
echo I am making gff files of annotated exons, genes and TSS respectively >&2
awk '$3=="exon"' $annot | awk -f $MAKEOK | awk -f $GFF2GFF | sort -k12,12 -k4,4n -k5,5n > annot_exons.gff
awk '$3=="gene"' $annot | awk -f $MAKEOK > annot_genes.gff
awk -v toadd=transcript -v fldno=12 -f $BOUNDARIES annot_exons.gff | awk -v fileRef=annot_exons.gff 'BEGIN{while (getline < fileRef >0){gnid[$12]=$10}} {print $1, $2, $3, $4, $5, $6, $7, $8, "gene_id", gnid[$10], $9, $10}' | awk -f $GFF2GFF > annot_tr.gff
$MAKETSS $annot
echo done >&2

# 2. Make a complete file from the prediction
#############################################
echo I am making a complete file for the predictions >&2
$MAKESUM $mytr > prediction_sumstat.txt
echo done >&2

# 3. Compute the number of exons of each transcript to be compared to the annotation
#####################################################################################
echo I am computing the number of exons of each transcript to be compared to the annotation >&2
awk '$3=="exon"' $mytr | awk -f $MAKEOK > $mytrbase\_exons.gff
awk '{nbex[$12]++}END{for(t in nbex){print t, nbex[t];}}' $mytrbase\_exons.gff > trid_nbex.txt
echo done >&2

# 4. I am running comptr on the input transcripts and the annotation
####################################################################
echo I am running comptr >&2
$COMPTR $mytrbase\_exons.gff annot_exons.gff -s 5000 -o $mytrbase\_vs_annot.tsv
echo The classes from comptr >&2
awk '{print $2}' $mytrbase\_vs_annot.tsv | sort | uniq -c | sort -k1,1nr >&2
echo done >&2

# 5. Refine the comptr class using the following rules/ideas:
#############################################################
# For the Overlap and the Monoexonic class, I need to know whether their extent overlaps the annotated genes on the same strand
# and if not whether they overlap them on the opposite strand, otherwise they will be considered intergenic if they overlap the 
# annotated genes on the same strand I also need to know whether they are totally included in them or extend them. Also in the
# Intergenic_or_antisense class I need to know how many are ig and how many are antisense so in fact I will collapse the Monoexonic + 
# Overlap + Intergenic_or_antisense into one class and I will overlap their whole extent with annotated gene extents in both a stranded 
# and an unstranded way, and then for the ones strandedly overlapping I will also check whether this is extension or inclusion
# 5.a. Make a transcript file with the monoexonic, overlap and intergenic_antisense class to retrive the ones that are extension or inclusion
#############################################################################################################################################
#      of the annotation on the same strand (but not respecting the structure), the ones that are antisense
###########################################################################################################
echo I am refining the class for the transcripts considered as monoexonic, overlap or intergenic_antisense by comptr >&2
awk '$2=="Monoexonic"||$2=="Overlap"||$2=="Intergenic_or_antisense"{split($1,a,"\""); print a[2]}' $mytrbase\_vs_annot.tsv > monoex_overlap_igas_trids.txt 
# 4
# 140 (1 fields)
echo The number of transcripts that are monoexonic, overlap or intergenic_antisense according to comptr >&2
wc -l monoex_overlap_igas_trids.txt | awk '{print $1}' >&2

awk -v fileRef=monoex_overlap_igas_trids.txt 'BEGIN{while (getline < fileRef >0){ok["\""$1"\";"]=1}} (($3=="transcript")&&(ok[$12]==1))' $mytrbase\_complete.gff > $mytrbase\_tr_monoex_overlap_igas.gff
# chrX	Cufflinks	transcript	41782188	41782387	0.0	-	.	gene_id "SID38216-SID38217"-v.1121.1";";,"SID38222-SID38223"-v.725.1";"; transcript_id 208
# 140 (12 fields)
$OVERLAP $mytrbase\_tr_monoex_overlap_igas.gff annot_genes.gff -st 1 -f gn -o $mytrbase\_tr_monoex_overlap_igas_strover_genes.gff
# chr1	Cufflinks	transcript	22263695	22263843	0.0	-	.	 gene_id "SID38212-SID38213"-v.9.1";";,"SID38204-SID38205"-v.22.1";"; transcript_id 72 ov_gn: 1
# 140 (14 fields)
$OVERLAP $mytrbase\_tr_monoex_overlap_igas_strover_genes.gff annot_genes.gff -f gn -o $mytrbase\_tr_monoex_overlap_igas_strover_genes_over_genes.gff
# chr1	Cufflinks	transcript	22263695	22263843	0.0	-	.	 gene_id "SID38212-SID38213"-v.9.1";";,"SID38204-SID38205"-v.22.1";"; transcript_id 72 ov_gn: 1 ov_gn: 1
# 140 (16 fields)

awk '$(NF-2)==1' $mytrbase\_tr_monoex_overlap_igas_strover_genes_over_genes.gff > $mytrbase\_tr_monoex_overlap_igas_strover_genes_ok_over_genes.gff
# chr1	Cufflinks	transcript	22263695	22263843	0.0	-	.	 gene_id "SID38212-SID38213"-v.9.1";";,"SID38204-SID38205"-v.22.1";"; transcript_id 72 ov_gn: 1 ov_gn: 1
# 122 (16 fields) 
echo The number of those that are overlapping annotated genes on the same strand >&2
wc -l $mytrbase\_tr_monoex_overlap_igas_strover_genes_ok_over_genes.gff | awk '{print $1}' >&2

$OVERLAP $mytrbase\_tr_monoex_overlap_igas_strover_genes_ok_over_genes.gff annot_genes.gff -st 1 -i 1 -f gn -o $mytrbase\_tr_monoex_overlap_igas_strover_genes_ok_over_genes_totstrincl.gff
# chr1	Cufflinks	transcript	22263695	22263843	0.0	-	.	 gene_id "SID38212-SID38213"-v.9.1";";,"SID38204-SID38205"-v.22.1";"; transcript_id 72 ov_gn: 1 ov_gn: 1 i1_gn: 0
# 122 (18 fields)
echo "Their classification into inclusion (1) and extension (0)" >&2
awk '{print $NF}' $mytrbase\_tr_monoex_overlap_igas_strover_genes_ok_over_genes_totstrincl.gff | sort | uniq -c >&2
awk '$NF==1{print $12}' $mytrbase\_tr_monoex_overlap_igas_strover_genes_ok_over_genes_totstrincl.gff > included_but_not_respecting_annot_tr_structure.txt
# 7
# 14 (1 fields)
awk '$NF==0{print $12}' $mytrbase\_tr_monoex_overlap_igas_strover_genes_ok_over_genes_totstrincl.gff > extension_but_not_respecting_annot_tr_structure.txt
# 72
# 108 (1 fields)

awk '(($NF==1)&&($(NF-2)==0)){print $12}' $mytrbase\_tr_monoex_overlap_igas_strover_genes_over_genes.gff > as_overlapping.txt
# 467  
# 2 (1 fields) *** here probably few because we want an overlap by at least 1 bp and do not allow to be at some distance
echo "The number of antisense transcripts (at least 1 bp overlap) from those that are monoexonic, overlap or intergenic_antisense according to comptr" >&2
wc -l as_overlapping.txt | awk '{print $1}' >&2

# 5.b. From the transcripts that are monoexonic, overlap or intergenic_antisense according to comptr but that are not overlapping the annotation in any way 
##########################################################################################################################################################
#      look whether their TSS is close enough (500bp) to an annotated TSS on the other strand in order to classify them as divergent from the annotation
########################################################################################################################################################
# chr1	Cufflinks	transcript	22263695	22263843	0.0	-	.	 gene_id "SID38212-SID38213"-v.9.1";";,"SID38204-SID38205"-v.22.1";"; transcript_id 72 ov_gn: 1 ov_gn: 1
# 140 (16 fields)
awk -v ext=500 'BEGIN{OFS="\t"}(($NF==0)&&($(NF-2)==0)){if($7=="+"){pos=$4}else{pos=$5} if(pos-ext>0){print $1, ".", ".", pos-ext+1, pos+ext, ".", $7, ".", "transcript_id", $12}}' $mytrbase\_tr_monoex_overlap_igas_strover_genes_over_genes.gff | awk -f $GFF2GFF > $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext.gff
# old
# awk -v ext=500 'BEGIN{OFS="\t"}(($NF==0)&&($(NF-2)==0)){if($7=="+"){pos=$4}else{pos=$5} if(pos-ext>0){print $1, pos-ext, pos+ext, $12, ".", $7}}' $mytrbase\_tr_monoex_overlap_igas_strover_genes_over_genes.gff > $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext.bed
# chr1	147993729	147994729	4	.	+
# 16 (6 fields)
$OVERLAP $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext.gff $annbase\_capped_sites_nr.gff -st -1 -f cap -o $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext_overanntss.gff
# old
# $INTER -a $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext.bed -b $annbase\_capped_sites_nr.gff -S -wao -nonamecheck > $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext_overanntss.bed
# chr1	1677232	1678232	1240	.	-	.	.	.	-1	-1	-1	.	.	.	0
# 118 (16 fields)
# 110 (25 fields)
awk '$NF!=0{split($10,a,"\""); print a[2]}' $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext_overanntss.gff | sort | uniq > antisense_500bpwaytssannot.txt
# old
# awk '$NF!=0{print $4}' $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext_overanntss.bed | sort | uniq > antisense_500bpwaytssannot.txt
# 338
# 6 (1 fields)
echo "The number of divergent transcripts from those (500bp from annotated tss)" >&2
wc -l antisense_500bpwaytssannot.txt | awk '{print $1}' >&2
echo done >&2

# 5.c. Gather all the above information into a single tsv file which is an extension of the comptr file
#######################################################################################################
# !!! if the broad class is Unstranded make the refined class unstranded !!!
# !!! would be good to encapsulate this code into a proper awk script !!!
echo I am gathering all the collected information into a single tsv file which is an extension of the comptr file >&2
awk -v fileRef1=included_but_not_respecting_annot_tr_structure.txt -v fileRef2=extension_but_not_respecting_annot_tr_structure.txt -v fileRef3=as_overlapping.txt -v fileRef4=antisense_500bpwaytssannot.txt -v fileRef5=trid_nbex.txt 'BEGIN{OFS="\t"; print "trid", "comptrclass", "annottrlist", "refinedtrclass", "nbex"; while (getline < fileRef1 >0){incl[$1]=1} while (getline < fileRef2 >0){ext[$1]=1} while (getline < fileRef3 >0){asover[$1]=1} while (getline < fileRef4 >0){asaway[$1]=1} while (getline < fileRef5 >0){nbex[$1]=$2}} {split($1,a,"\""); if($2=="Unstranded"){class="unstranded"}else{if(($2=="Exact")||($2=="Inclusion")||(incl[$1]==1)){class="annot"}else{if(($2=="Extension")||(ext[$1]==1)){class="extension"}else{if((asover[$1]==1)||(asaway[$1]==1)){class="antisense"}else{class="intergenic"}}}} gsub(/\;/,"",$3); gsub(/\"/,"",$3); print "t"a[2], $2, $3, class, "n"nbex[$1]}' $mytrbase\_vs_annot.tsv > $mytrbase\_comp_refinedclass_nbex.tsv
# trid	comptrclass	gencv19trlist	trclass	nbex
# t214	Extension	ENST00000429744.2,	extension	n2
# 230 (5 fields)
echo Number of initial transcripts from the different refined classes >&2
awk 'NR>=2{print $4}'  $mytrbase\_comp_refinedclass_nbex.tsv | sort | uniq -c | sort -k1,1nr >&2
echo done >&2

# 6. Add the intermediate class and gene list wrt annot, the class being annot, extension, novel_of_annot or novel: 
###################################################################################################################
# 1) annot = predicted transcript that have exactly the same exon structure as a reference tr but that do not
#    extend the annotated transcript on neither side (I can use my Exact class but add the constraint that there
#    is no extension
# 2) extension = predicted transcript that are not in 1) but that have exactly the same exon structure as    
#    a reference tr (I can use my Exact class but that are not in class 1)
# 3) novel_of_annot = new isoform or variant of reference genes = predicted transcripts that are not in 1) or 2) but that
#    have at least one common intron with a reference annotated tr (I have to compute it from scratch by
#    asking for one common intron same strand as the reference and then ask that the transcript is neither
#    in class 1 nor in class 2                                                 
# 4) novel = new transcript = the ones not in 1) or 2) or 3)
echo I am adding the intermediate class and gene list wrt annot, the class being annot, extension, novel_of_annot or novel >&2
# 6.a. make the introns from the reference gene annotation with gene id
#######################################################################
awk -v fldgn=10 -v fldtr=12 -f $MAKEINTRONS annot_exons.gff | awk '{split($10,a,"\""); print $1":"$4":"$5":"$7, a[2]}' | sort | uniq | awk '{gnlist[$1]=(gnlist[$1])($2)(",")}END{for(i in gnlist){print i, gnlist[i]}}' > annot_introns_coord_gnlist.tsv
# 2:132525125:132526321:+ ENSBTAG00000001082,
# 3:10030765:10030858:- ENSBTAG00000004776,
# 184067 (2 fields) 
# 6.b. make the introns from the predicted transcripts
#######################################################
sort -k12,12 -k4,4n -k5,5n $mytrbase\_exons.gff | awk -v fldgn=10 -v fldtr=12 -f $MAKEINTRONS > $mytrbase\_introns.gff
# 1	Cufflinks	intron	242647	254558	.	+	.	gene_id "XLOC_000002"; transcript_id "TCONS_00000002";  
# 1	Cufflinks	intron	242647	275622	.	+	.	gene_id "XLOC_000002"; transcript_id "TCONS_00000003";  
# 660657 (12 fields)
# 6.c. for each predicted transcript make the list of reference genes sharing an intron with it
###############################################################################################
awk -v fileRef=annot_introns_coord_gnlist.tsv 'BEGIN{OFS="\t"; while (getline < fileRef >0){gnlist[$1]=$2}} {if(gnlist[$1":"$4":"$5":"$7]!=""){split($12,a,"\""); split(gnlist[$1":"$4":"$5":"$7],b,","); k=1; while(b[k]!=""){ok[a[2],b[k]]++; if(ok[a[2],b[k]]==1){nbgn[a[2]]++; gn[a[2],nbgn[a[2]]]=b[k];} k++}}} END{for(t in nbgn){for(i=1; i<=nbgn[t]; i++){s[t]=s[t](gn[t,i])(",")} print t, s[t]}}' $mytrbase\_introns.gff > $mytrbase\_transcriptid_gnlistwithcommonintron.tsv
# TCONS_00083097	ENSBTAG00000019555,
# TCONS_00083098	ENSBTAG00000019555,
# 54411 (2 fields)
# 6.d. add the intermediate class to each transcript based on this info and others as well as the list of associated genes
##########################################################################################################################
awk -v fileRef1=$mytrbase\_complete.gff -v fileRef2=annot_tr.gff -v fileRef3=$mytrbase\_transcriptid_gnlistwithcommonintron.tsv -f $ADDCLASS $mytrbase\_comp_refinedclass_nbex.tsv > $mytrbase\_comp_refinedclass_nbex_intermclass.tsv
# trid	comptrclass	annottrlist	refinedtrclass	nbex	interm_class	interm_gnlist
# tTCONS_00000002	Intergenic_or_antisense	.	intergenic	n2	novel	.
# 84972 (7 fields)

echo Number of initial transcripts from the different intermediate classes >&2
awk 'NR>=2{print $6}'  $mytrbase\_comp_refinedclass_nbex_intermclass.tsv | sort | uniq -c | sort -k1,1nr >&2
echo done >&2

# 7. Delete intermediate files
##############################
echo I am deleting intermediate files >&2
rm annot_exons.gff annot_genes.gff annot_tr.gff 
rm $mytrbase\_exons.gff $mytrbase\_introns.gff
rm $annbase\_capped_sites_nr.gff
rm $annbase\_capped_sites.gff
rm trid_nbex.txt
rm $mytrbase\_vs_annot.tsv
rm monoex_overlap_igas_trids.txt 
rm $mytrbase\_tr_monoex_overlap_igas.gff
rm $mytrbase\_tr_monoex_overlap_igas_strover_genes.gff 
rm $mytrbase\_tr_monoex_overlap_igas_strover_genes_over_genes.gff
rm $mytrbase\_tr_monoex_overlap_igas_strover_genes_ok_over_genes.gff
rm $mytrbase\_tr_monoex_overlap_igas_strover_genes_ok_over_genes_totstrincl.gff
rm included_but_not_respecting_annot_tr_structure.txt
rm extension_but_not_respecting_annot_tr_structure.txt
rm as_overlapping.txt
rm $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext.gff
rm $mytrbase\_tr_monoex_overlap_igas_over_genes_ko_tss_500bpext_overanntss.gff
rm antisense_500bpwaytssannot.txt
rm annot_introns_coord_gnlist.tsv
rm $mytrbase\_introns.gff $mytrbase\_transcriptid_gnlistwithcommonintron.tsv
echo done >&2

