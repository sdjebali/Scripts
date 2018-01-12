# add_intermclass_andgnlist.awk
# this script takes as input
# - a complete file of predicted transcripts in gff format (includes exons and transcripts)
# - an annotated transcript gff file in gff format
# - a 2 column tsv file with for each predicted transcript with an intron common with an annotated transcript, its id and the list of these transcript genes
# - a 5 column tsv file that was the previous output of refine_comptr_output.sh which has
#   predicted transcript id with t in front of it, its comptr class, its associated annotated transcript list, its refined class and its nb ex

# and outputs the last file but with 2 additional fields:
##########################################################
# - its intermediate class which is one of those 4:
#   1) annotated = predicted transcript that have exactly the same exon structure as a reference tr but with
#      +-10bp in 5' and/or 3' (I can use my Exact class but add the constraint that the distance between the
#      5' end of the predicted transcript most 5' exon and the 5' end of the annotated transcript most 5' exon
#      is not more than 10bp and that the distance between the 3' end of the predicted transcript most 3' exon
#      and the 3' end of the annotated transcript most 3' exon is not more than 10bp                                                                         
#   2) extension = predicted transcript that are not in 1) but that have exactly the same exon structure as    
#      a reference tr but that extend more than 10bp in 5' and/or in 3' (I can use my Exact class but that are
#      not in class 1)
#   3) novel_of_annot = new isoform or variant of reference genes = predicted transcripts that are not in 1) or 2) but that
#      have at least one common intron with a reference annotated tr (I have to compute it from scratch by
#      asking for one common intron same strand as the reference and then ask that the transcript is neither
#      in class 1 nor in class 2                                                 
#   4) novel = new transcript = the ones not in 1) or 2) or 3)
# - the list of corresponding annotated genes

# example
#########
# cd /work/project/fragencode/results/rnaseq/bos_taurus/assembled
# time awk -v fileRef1=bos_taurus_cuff_tpm0.1_2sample_complete_complete.gff -v fileRef2=annot_tr.gff -v fileRef3=bos_taurus_cuff_tpm0.1_2sample_complete_transcriptid_gnlistwithcommonintron.tsv -f $ADDCLASS bos_taurus_cuff_tpm0.1_2sample_complete_comp_refinedclass_nbex.tsv > bos_taurus_cuff_tpm0.1_2sample_complete_comp_refinedclass_nbex_intermclass.tsv

# input files
#############
# bos_taurus_cuff_tpm0.1_2sample_complete_complete.gff
# 1	Cufflinks	exon	242243	242646	.	+	.	gene_id "XLOC_000002"; transcript_id "TCONS_00000002"; 
# 1	Cufflinks	exon	254559	261079	.	+	.	gene_id "XLOC_000002"; transcript_id "TCONS_00000002"; 
# 1525552 (12 fields)

# annot_tr.gff
# 1	ensembl	transcript	154923159	154967566	.	+	.	gene_id "ENSBTAG00000008240"; transcript_id "ENSBTAT00000010847";
# 22	ensembl	transcript	50638075	50640896	.	-	.	gene_id "ENSBTAG00000000478"; transcript_id "ENSBTAT00000000605";
# 26740 (12 fields)

# bos_taurus_cuff_tpm0.1_2sample_complete_transcriptid_gnlistwithcommonintron.tsv
# TCONS_00083097	ENSBTAG00000019555
# TCONS_00083098	ENSBTAG00000019555
# 54411 (2 fields)

# bos_taurus_cuff_tpm0.1_2sample_complete_comp_refinedclass_nbex.tsv
# trid	comptrclass	annottrlist	refinedtrclass	nbex
# tTCONS_00000002	Intergenic_or_antisense	.	intergenic	n2
# 84972 (5 fields)

# output file
#############
# trid	comptrclass	annottrlist	refinedtrclass	nbex	interm_class	interm_gnlist
# tTCONS_00000002	Intergenic_or_antisense	.	intergenic	n2	novel	.
# 84972 (7 fields)
# awk 'NR>=2{print $6}' bos_taurus_cuff_tpm0.1_2sample_complete_comp_refinedclass_nbex_intermclass.tsv | sort | uniq -c | sort -k1,1nr
# 40813 novel_of_annot
# 30560 novel
# 11085 annot
# 2513 extension

BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	if($3=="transcript")    # for each predicted transcript, store its gene id as well as its beg and end (beg<end)
	{
	    split($10,a,"\"");
	    split($12,b,"\"");
	    gn["t"b[2]]=a[2];
	    gbeg["t"b[2]]=$4;
	    gend["t"b[2]]=$5;
	}
    }
    while (getline < fileRef2 >0)    # for each annotated transcript, store its gene id as well as its beg and end (beg<end)
    {
	split($10,a,"\"");
	split($12,b,"\"");
	gn2[b[2]]=a[2];
	gbeg2[b[2]]=$4;
	gend2[b[2]]=$5;
    }
    while (getline < fileRef3 >0)   # for each predicted transcript, store the nr list of annotated genes with an intron common with it
    {
	commonintron_gnlist["t"$1]=$2;
    }
}

NR==1{
    # print the header with 2 additional labels
    print $0, "interm_class", "interm_gnlist";   
}

NR>=2{
    # for each predicted transcript, find its intermediate class together with the corresponding gene list
    class="";
    gnlist="";
    if($2=="Exact")   # if the predicted transcript is from the comptr exact class
    {
	split($3,a,",");  
	k=1;
	while(a[k]!="")   # for each annotated transcript corresponding to this predicted transcript
	{
	    if((abs(gbeg[$1]-gbeg2[a[k]])<=10)&&(abs(gend[$1]-gend2[a[k]])<=10))  # look if the two annotated transcript ends are less than 10bp away from the two predicted transcript ends respectively
		                                                                  # in this case it is of the class annot but only associated to this transcript
	    {
		class="annot";   # remember it is of class "annot" but just associated to annotated transcript number k
		ok[$1,k]=1;      # and remember the number of the annotated transcript associated to it
	    }
	    k++;
	}
	if(class=="annot")
	{
	    for(i=1; i<=k; i++)
	    {
		if(ok[$1,i]==1)
		{
		    gnlist=(gnlist)(gn2[a[k]])(",");
		}
	    }
	}
	else
	{
	    class="extension";
	    for(i=1; i<=k; i++)
	    {
		gnlist=(gnlist)(gn2[a[k]])(",");
	    }
	}
    }
    else   # if the predicted transcript is not from the comptr exact class 
    {
	if(commonintron_gnlist[$1]!="") # if it has an intron in common with an annotated gene
	{
	    class="novel_of_annot";
	    gnlist=commonintron_gnlist[$1];
	}
	else
	{
	    class="novel";
	    gnlist=".";
	}
    }
    print $0, class, gnlist;
}


function abs(x){
    return (x>=0 ? x : -1*x);
}
