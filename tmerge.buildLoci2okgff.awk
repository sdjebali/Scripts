# tmerge.buildLoci2okgff.awk
# !!! bug: ref_gene_id of a gene should be the union of all the ref_gene_id of the tr of the gene, and not just of the 1st tr of the gene !!!
# this script takes as input:
#############################
# - an exon gff file from tmerge and buildLoci that has been formatted so that the 9th field
#   starts with gene_id and transcript_id, and where contains is a field that has the set of all transcripts
#   of different experiments and the ref annot that support the transcript of the exon (with an intron chain either
#   included or equal to the one of the current transcript)
# - a ref gene annotation in gtf format
# this script will output a more complete gff file aftering doing the following actions:
########################################################################################
# - adding transcript and gene rows to the file
# - add a ref_gene_id field with . or a gene or a comma separated gene list according to whether the transcript is supported by at least one ref tr
# - report the verbose info of the exon rows to the transcript rows and remove it from the exon rows
# - have gene_id, transcript_id and ref_gene_id in all rows
# - have a transcript id corresponding to a ref tr id whenever the tmerge transcript EXACTLY
#   corresponds to a tr of the ref (taking the longest in case there are many)

# exemple of usage:
###################
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/pipelines/rnaseq/debugging/apr11th2022
# indir=/work/project/fragencode/workspace/ckurylo/geneswitch/TAGADA/test_results/fragencode_stringtie_vs_tmerge/work/05/b1016a06d86b37d258ff41b5066a11
# pgm1=~/fragencode/tools/multi/Scripts/make_gff_ok.awk
# time cat $indir/results/tmerged.gene_id.gtf | awk -f $pgm1 > tmp.gff
# real	3m11.568s *** file is 2.4G
# pgm2=~/fragencode/tools/multi/Scripts/tmerge.buildLoci2okgff.awk
# pgm3=~/fragencode/tools/multi/Scripts/gff2gff.awk
# time awk -v fileRef1=$indir/sus_scrofa.gtf -v fileRef2=tmp.gff -f $pgm2 tmp.gff | awk -f $pgm3 | sort -k1,1 -k4,4n -k5,5n > novel.gtf
# real	3m34.295s
# rm tmp.gff

# fileRef1=$indir/sus_scrofa.gtf
# #!genome-build Sscrofa11.1
# #!genome-version Sscrofa11.1
# 1	ensembl	gene	1	3782	.	+	.	gene_id "ENSSSCG00000048769"; gene_version "1"; gene_source "ensembl"; gene_biotype "lncRNA";
# 1	ensembl	transcript	1	3780	.	+	.	gene_id "ENSSSCG00000048769"; gene_version "1"; transcript_id "ENSSSCT00000066540"; transcript_version "1"; gene_source "ensembl"; gene_biotype "lncRNA"; transcript_source "ensembl"; transcript_biotype "lncRNA";
# 1	ensembl	exon	1	2465	.	+	.	gene_id "ENSSSCG00000048769"; gene_version "1"; transcript_id "ENSSSCT00000066540"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "lncRNA"; transcript_source "ensembl"; transcript_biotype "lncRNA"; exon_id "ENSSSCE00000399651"; exon_version "1";
# 5 (2 fields)
# 15400 (16 fields)
# 16508 (18 fields)
# 42269 (24 fields)
# 17922 (26 fields)
# 164412 (28 fields)
# 225950 (30 fields)
# 3221 (32 fields)
# 1038976 (34 fields)
# 6320 (36 fields)

# after $pgm1 we get tmp.gff that goes as fileRef2 ***and*** input file
# 1	tmerge	exon	1	2465	0	+	.	gene_id "LOC_000000050668"; transcript_id "TM_000000000001"; contains "duodenum_pig4.filtered:STRG.1.# 1	tmerge	exon	1	2465	0	+	.	gene_id "LOC_000000050668"; transcript_id "TM_000000000001"; contains "duodenum_pig4.filtered:STRG.1.1,liver_pig2.filtered:STRG.1.1,cd8_pig1.filtered:STRG.1.1,cd8_pig2.filtered:STRG.1.1,cerebellum_pig2.filtered:STRG.1.1,testis_pig2.filtered:STRG.1.1,ileum_pig4.filtered:STRG.1.1,cerebellum_pig1.filtered:STRG.1.1,muscle_pig4.filtered:STRG.1.1,muscle_pig3.filtered:STRG.1.1,ileum_pig3.filtered:STRG.1.1,kidney_pig4.filtered:STRG.1.1,cd8_pig4.filtered:STRG.1.1,testis_pig1.filtered:STRG.1.1,duodenum_pig1.filtered:STRG.1.1,cd4_pig3.filtered:STRG.1.1,ref:ENSSSCT00000066540,cerebellum_pig3.filtered:STRG.1.1,duodenum_pig3.filtered:STRG.1.1,muscle_pig1.filtered:STRG.1.1,liver_pig4.filtered:STRG.1.1,liver_pig1.filtered:STRG.1.1,kidney_pig3.filtered:STRG.1.1,muscle_pig2.filtered:STRG.1.1,cerebellum_pig4.filtered:STRG.1.1,cd8_pig3.filtered:STRG.1.1,lung_pig2.filtered:STRG.1.1,ileum_pig2.filtered:STRG.1.1,lung_pig1.filtered:STRG.1.1,duodenum_pig2.filtered:STRG.1.1,cd4_pig4.filtered:STRG.1.1,liver_pig3.filtered:STRG.1.1,ileum_pig1.filtered:STRG.1.1,cd4_pig2.filtered:STRG.1.1"; contains_count "34"; 3p_dists_to_3p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; 5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; flrpm "57.190992"; longest "duodenum_pig4.filtered:STRG.1.1,liver_pig2.filtered:STRG.1.1,cd8_pig1.filtered:STRG.1.1,cd8_pig2.filtered:STRG.1.1,cerebellum_pig2.filtered:STRG.1.1,testis_pig2.filtered:STRG.1.1,ileum_pig4.filtered:STRG.1.1,cerebellum_pig1.filtered:STRG.1.1,muscle_pig4.filtered:STRG.1.1,muscle_pig3.filtered:STRG.1.1,ileum_pig3.filtered:STRG.1.1,kidney_pig4.filtered:STRG.1.1,cd8_pig4.filtered:STRG.1.1,testis_pig1.filtered:STRG.1.1,duodenum_pig1.filtered:STRG.1.1,cd4_pig3.filtered:STRG.1.1,ref:ENSSSCT00000066540,cerebellum_pig3.filtered:STRG.1.1,duodenum_pig3.filtered:STRG.1.1,muscle_pig1.filtered:STRG.1.1,liver_pig4.filtered:STRG.1.1,liver_pig1.filtered:STRG.1.1,kidney_pig3.filtered:STRG.1.1,muscle_pig2.filtered:STRG.1.1,cerebellum_pig4.filtered:STRG.1.1,cd8_pig3.filtered:STRG.1.1,lung_pig2.filtered:STRG.1.1,ileum_pig2.filtered:STRG.1.1,lung_pig1.filtered:STRG.1.1,duodenum_pig2.filtered:STRG.1.1,cd4_pig4.filtered:STRG.1.1,liver_pig3.filtered:STRG.1.1,ileum_pig1.filtered:STRG.1.1,cd4_pig2.filtered:STRG.1.1"; longest_FL_supporters "liver_pig2.filtered:STRG.1.1,testis_pig2.filtered:STRG.1.1,cerebellum_pig2.filtered:STRG.1.1,ileum_pig4.filtered:STRG.1.1,cd4_pig2.filtered:STRG.1.1,ileum_pig3.filtered:STRG.1.1,muscle_pig3.filtered:STRG.1.1,kidney_pig4.filtered:STRG.1.1,cd8_pig4.filtered:STRG.1.1,cerebellum_pig3.filtered:STRG.1.1,kidney_pig3.filtered:STRG.1.1,muscle_pig2.filtered:STRG.1.1,cd8_pig3.filtered:STRG.1.1,ileum_pig2.filtered:STRG.1.1,duodenum_pig2.filtered:STRG.1.1,cd4_pig4.filtered:STRG.1.1,liver_pig3.filtered:STRG.1.1,ileum_pig1.filtered:STRG.1.1,duodenum_pig4.filtered:STRG.1.1,cd8_pig1.filtered:STRG.1.1,cd8_pig2.filtered:STRG.1.1,cerebellum_pig1.filtered:STRG.1.1,muscle_pig4.filtered:STRG.1.1,duodenum_pig1.filtered:STRG.1.1,cd4_pig3.filtered:STRG.1.1,testis_pig1.filtered:STRG.1.1,ref:ENSSSCT00000066540,duodenum_pig3.filtered:STRG.1.1,muscle_pig1.filtered:STRG.1.1,liver_pig4.filtered:STRG.1.1,liver_pig1.filtered:STRG.1.1,cerebellum_pig4.filtered:STRG.1.1,lung_pig2.filtered:STRG.1.1,lung_pig1.filtered:STRG.1.1"; longest_FL_supporters_count "34"; mature_RNA_length "3128"; meta_3p_dists_to_5p "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1"; meta_5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; rpm "57.190992"; spliced "1"; 
# 2169536 (38 fields)  *** file is 2.4G and take 3 minutes to get


# output file is like this
# 1	tagada	exon	1	2465	0	+	.	gene_id "LOC_000000050668"; transcript_id "ENSSSCT00000066540"; ref_gene_id "ENSSSCG00000048769"; tmerge_tr_id "TM_000000000001";
# 1	tagada	transcript	1	3780	0	+	.	gene_id "LOC_000000050668"; transcript_id "ENSSSCT00000066540"; contains "duodenum_pig4.filtered:STRG.1.1,liver_pig2.filtered:STRG.1.1,cd8_pig1.filtered:STRG.1.1,cd8_pig2.filtered:STRG.1.1,cerebellum_pig2.filtered:STRG.1.1,testis_pig2.filtered:STRG.1.1,ileum_pig4.filtered:STRG.1.1,cerebellum_pig1.filtered:STRG.1.1,muscle_pig4.filtered:STRG.1.1,muscle_pig3.filtered:STRG.1.1,ileum_pig3.filtered:STRG.1.1,kidney_pig4.filtered:STRG.1.1,cd8_pig4.filtered:STRG.1.1,testis_pig1.filtered:STRG.1.1,duodenum_pig1.filtered:STRG.1.1,cd4_pig3.filtered:STRG.1.1,ref:ENSSSCT00000066540,cerebellum_pig3.filtered:STRG.1.1,duodenum_pig3.filtered:STRG.1.1,muscle_pig1.filtered:STRG.1.1,liver_pig4.filtered:STRG.1.1,liver_pig1.filtered:STRG.1.1,kidney_pig3.filtered:STRG.1.1,muscle_pig2.filtered:STRG.1.1,cerebellum_pig4.filtered:STRG.1.1,cd8_pig3.filtered:STRG.1.1,lung_pig2.filtered:STRG.1.1,ileum_pig2.filtered:STRG.1.1,lung_pig1.filtered:STRG.1.1,duodenum_pig2.filtered:STRG.1.1,cd4_pig4.filtered:STRG.1.1,liver_pig3.filtered:STRG.1.1,ileum_pig1.filtered:STRG.1.1,cd4_pig2.filtered:STRG.1.1"; contains_count "34"; 3p_dists_to_3p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; 5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; flrpm "57.190992"; longest "duodenum_pig4.filtered:STRG.1.1,liver_pig2.filtered:STRG.1.1,cd8_pig1.filtered:STRG.1.1,cd8_pig2.filtered:STRG.1.1,cerebellum_pig2.filtered:STRG.1.1,testis_pig2.filtered:STRG.1.1,ileum_pig4.filtered:STRG.1.1,cerebellum_pig1.filtered:STRG.1.1,muscle_pig4.filtered:STRG.1.1,muscle_pig3.filtered:STRG.1.1,ileum_pig3.filtered:STRG.1.1,kidney_pig4.filtered:STRG.1.1,cd8_pig4.filtered:STRG.1.1,testis_pig1.filtered:STRG.1.1,duodenum_pig1.filtered:STRG.1.1,cd4_pig3.filtered:STRG.1.1,ref:ENSSSCT00000066540,cerebellum_pig3.filtered:STRG.1.1,duodenum_pig3.filtered:STRG.1.1,muscle_pig1.filtered:STRG.1.1,liver_pig4.filtered:STRG.1.1,liver_pig1.filtered:STRG.1.1,kidney_pig3.filtered:STRG.1.1,muscle_pig2.filtered:STRG.1.1,cerebellum_pig4.filtered:STRG.1.1,cd8_pig3.filtered:STRG.1.1,lung_pig2.filtered:STRG.1.1,ileum_pig2.filtered:STRG.1.1,lung_pig1.filtered:STRG.1.1,duodenum_pig2.filtered:STRG.1.1,cd4_pig4.filtered:STRG.1.1,liver_pig3.filtered:STRG.1.1,ileum_pig1.filtered:STRG.1.1,cd4_pig2.filtered:STRG.1.1"; longest_FL_supporters "liver_pig2.filtered:STRG.1.1,testis_pig2.filtered:STRG.1.1,cerebellum_pig2.filtered:STRG.1.1,ileum_pig4.filtered:STRG.1.1,cd4_pig2.filtered:STRG.1.1,ileum_pig3.filtered:STRG.1.1,muscle_pig3.filtered:STRG.1.1,kidney_pig4.filtered:STRG.1.1,cd8_pig4.filtered:STRG.1.1,cerebellum_pig3.filtered:STRG.1.1,kidney_pig3.filtered:STRG.1.1,muscle_pig2.filtered:STRG.1.1,cd8_pig3.filtered:STRG.1.1,ileum_pig2.filtered:STRG.1.1,duodenum_pig2.filtered:STRG.1.1,cd4_pig4.filtered:STRG.1.1,liver_pig3.filtered:STRG.1.1,ileum_pig1.filtered:STRG.1.1,duodenum_pig4.filtered:STRG.1.1,cd8_pig1.filtered:STRG.1.1,cd8_pig2.filtered:STRG.1.1,cerebellum_pig1.filtered:STRG.1.1,muscle_pig4.filtered:STRG.1.1,duodenum_pig1.filtered:STRG.1.1,cd4_pig3.filtered:STRG.1.1,testis_pig1.filtered:STRG.1.1,ref:ENSSSCT00000066540,duodenum_pig3.filtered:STRG.1.1,muscle_pig1.filtered:STRG.1.1,liver_pig4.filtered:STRG.1.1,liver_pig1.filtered:STRG.1.1,cerebellum_pig4.filtered:STRG.1.1,lung_pig2.filtered:STRG.1.1,lung_pig1.filtered:STRG.1.1"; longest_FL_supporters_count "34"; mature_RNA_length "3128"; meta_3p_dists_to_5p "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1"; meta_5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; rpm "57.190992"; spliced "1"; ref_gene_id "ENSSSCG00000048769"; tmerge_tr_id "TM_000000000001";
# 57645 (12 fields)    *** gene rows always have 12 fields, for gene_id and ref_gene_id
# 1553055 (14 fields)
# 616481 (16 fields)   *** exon rows can have 14 or 16 fields, for gene_id, transcript_id, ref_gene_id and same with tmerge_tr_id when tr id becomes from ref
# 197370 (40 fields)
# 60329 (42 fields)    *** transcript rows can have 40 or 42 fields, that are the 38 initial ones plus ref_gene_id and also with tmerge_tr_id when tr id becomes from ref, however gene_id is always a locus id and not a ref gene id, this is given in ref_gene_id field     


BEGIN{
    # here we read the ref gene annotation gtf file to get the correspondance between transcript id and gene id
    # as well as transcript length, begining and end
    while (getline < fileRef1 >0)
    {
	if($3=="transcript")
	{
	    split($0,a,"\t");
	    split(a[9],b," ");
	    k=1;
	    while(b[k]!="")
	    {
		if(b[k]=="transcript_id")
		{
		    split(b[k+1],c,"\"");
		    trid=c[2];
		    chr[trid]=$1;
		    gbeg[trid]=$4;
		    gend[trid]=$5;
		    str[trid]=$7;
		    lg[trid]=$5-$4+1;
		}
		else
		{
		    if(b[k]=="gene_id")
		    {
			split(b[k+1],c,"\"");
			gnid=c[2];
		    }
		}
		k+=2;
	    }
	    refgn[trid]=gnid;
	}
    }

    # here we read the tmerge file that only has exon rows and with gene_id and transcript_id as two first subfields of the 9th field
     while (getline < fileRef2 >0)
    {
	# when we see the transcript for the 1st time only, from its first exon row we put the gbeg and gend and then reading the next exon rows we update the tr beg and end
	split($0,a,"\t");
	split(a[9],b," ");
	k=1;
	while(b[k]!="")
	{
	    if(b[k]=="transcript_id")
	    {
		split(b[k+1],c,"\"");
		chr[c[2]]=$1;
		str[c[2]]=$7;
		if((gbeg[c[2]]=="")||($4<gbeg[c[2]]))
		{
		    gbeg[c[2]]=$4;
		}
		if((gend[c[2]]=="")||($5>gend[c[2]]))
		{
		    gend[c[2]]=$5;
		}
	    }
	    k+=2;
	}
    }
}

# here we read the tmerge file again but this time to print exon rows (tr and gn rows will be printed in the END part)
$3=="exon"{
    # it has to be noted that all the exon rows of a given transcript contain the same information in the 9th field
    # so looking for a ref tr in the 9th field and making the corresponding ref_gene_id can be done for the 1st exon of the transcript only

    # here we read each exon row from tmerge and buildLoci and keep the verbose info for the transcript rows (removing them from the exon rows)
    # we also compute the transcript and gene boundaries to make transcript and gene rows in the END part of the script
    # note that the tmerge file is supposed to only have exon rows but a test is done here just in case
    nbex[$12]++;
    nbex[$10]++;

    # looking for info is done just the 1st time we see the transcript
    # note that here the transcript id is the tmerge tr id ($12)
    if(nbex[$12]==1)
    {
	ct[$12]="";
	# found is here to know whether a reference transcript is present in the contains field
	found[$12]=0;
	reftrlg[$12]=0;
	refgnid[$12]="";
	reftrid[$12]="";
	split($0,a,"\t");
	split(a[9],b," ");
    
	# find out whether there is a ref transcript in the supporting transcripts and in this case change the transcript id to the longest of them
	# only if this transcript has exactly the same beg and end as the tmerge tr (otherwise it will not be exactly identical to it) (long process)
	# here just extract the content of the contains field and assign it to the current tmerge tr id ($12)
	k=1;
	while(found[$12]==0&&b[k]!="")
	{
	    if(b[k]=="contains")
	    {
		found[$12]=1;
		ct[$12]=b[k+1];
	    }
	    k+=2;
	}
    
	# if we find a reference transcript in the contains field, then we take the longest one
	# but we remember the gene ids from all the reference transcripts (only once if we have several tr from the same gene though)
	# the info of ref gene and ref tr is associated to the current tmerge tr (with id in $12)
	if((found[$12]==1)&&(ct[$12]~/ref/))
	{
	    split(ct[$12],c,"\"");
	    split(c[2],d,",");
	    k=1;
	    while(d[k]!="")
	    {
		split(d[k],e,":");
		if(e[1]=="ref")
		{
		    # if the ref gene has not been seen for the current transcript we add it
		    seen[$12,refgn[e[2]]]++;
		    if(seen[$12,refgn[e[2]]]==1)
		    {
			refgnid[$12]=(refgnid[$12])(refgn[e[2]])(",");
		    }
		    if((reftrlg[$12]==0)||(lg[e[2]]>reftrlg[$12]))
		    {
			reftrlg[$12]=lg[e[2]];
			reftrid[$12]=e[2];
		    }
		}
		k++;
	    }
	}
	
	# look whether the current ref tr has the same begining and end as the current tmerge tr (id in $12)
	# in this case it will replace the tmerge tr id in the 12th field and we will remember the tmerge tr id in a tmerge_tr_id field
	# we also remember which of the two has to be printed in the END part
	split($12,a,"\"");
	if((reftrid[$12]!="")&&(gbeg[reftrid[$12]]==gbeg[a[2]])&&(gend[reftrid[$12]]==gend[a[2]]))
	{
	    currtr[$12]="\""reftrid[$12]"\"\;";
	    toprint[currtr[$12]]=1;
	}
	else
	{
	    currtr[$12]=$12;
	    toprint[$12]=1;
	}

	# remember verbose information that is not gene_id and transcript_id to put in future transcript rows
	# for the currtr of the tmerge tr id (either ref tr or tmerge tr, see above)
	split($0,a,"\t");
	split(a[9],b," ");
	k=1;
	while(b[k]!="")
	{
	    if(b[k]!="gene_id"&&b[k]!="transcript_id")
	    {
		info[currtr[$12]]=(info[currtr[$12]])(b[k])(" ")(b[k+1])(" ");
	    }
	    k+=2;
	}
    }

    # print the current exon row with either the tmerge tr id or the ref tr id in $12 (see above)
    # in case the ref tr id is taken, keep the tmerge tr id in a tmerge_transcript_id field
    # first remove the last comma in the refgnid list
    s=(refgnid[$12]=="" ? "." : substr(refgnid[$12],1,(length(refgnid[$12])-1)));
    
    # in case the tmerge tr id is kept
    if(currtr[$12]==$12)
    {
	print $1, "tagada", "exon", $4, $5, $6, $7, $8, $9"\t"$10"\t"$11"\t"$12"\tref_gene_id \""s"\";";
	gene[$12]=$10;
    }
    # otherwise the ref tr id is kept and we remember the tmerge trid as well
    else
    {
	print $1, "tagada", "exon", $4, $5, $6, $7, $8, $9"\t"$10"\t"$11"\t"currtr[$12]"\tref_gene_id \""s"\"; tmerge_tr_id "$12;
	gene[currtr[$12]]=$10;
    }
    
    # compute and update the characteristics of the gene row to be printed in the end part
    chr[$10]=$1;
    str[$10]=$7;
    if((gbeg[$10]=="")||($4<gbeg[$10]))
    {
	gbeg[$10]=$4;
    }
    if((gend[$10]=="")||($5>gend[$10]))
    {
	gend[$10]=$5;
    }
    
    # only once for each gene we update the info of the gene with gene_id and ref_gene_id
    # !!! here this is not correct for ref_gene_id as for the gene it should be the union of all the ref_gene_ids of all the tr of the gene !!!
    if(nbex[$10]==1)
    {
	toprintg[$10]=1;
	info[$10]="gene_id "$10" ref_gene_id \""s"\";";
    }

    # only once for each tr (the one to print, either with a tmerge or a ref tr id) we update the info of the tr
    # with ref_gene_id and tmerge_tr_id if necesary
    if(nbex[$12]==1)
    {
	if(currtr[$12]==$12)
	{
	    info[$12]=info[$12]" ref_gene_id \""s"\";";
	}
	else
	{
	    info[currtr[$12]]=info[currtr[$12]]" ref_gene_id \""s"\"; tmerge_tr_id "$12;
	}
    }   
}

END{
    # tr has double quotes but chr, gbeg, gend and str for tr do not accept double quotes, so need to split
    for(tr in toprint)
    {
	split(tr,a,"\"");
	print chr[a[2]], "tagada", "transcript", gbeg[a[2]], gend[a[2]], 0, str[a[2]], ".", "gene_id "(gene[tr])" transcript_id "(tr)" "info[tr];
    }

    # gn has double quotes and chr, gnbeg, gend and str for genes accept double quotes
    for(gn in toprintg)
    {
	print chr[gn], "tagada", "gene", gbeg[gn], gend[gn], 0, str[gn], ".", info[gn];
    }
}
