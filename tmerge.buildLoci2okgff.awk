# tmerge.buildLoci2okgff.awk
# this script takes as input:
#############################
# - an exon gff file from tmerge and buildLoci that has been formatted so that the 9th field
#   starts with gene_id and transcript_id, and where contains is a field that has the set of all transcripts
#   of different experiments and the ref annot that support the transcript of the exon
# - a ref gene annotation in gtf format
# this script will output a more complete gff file aftering doing the following actions:
########################################################################################
# - adding transcript and gene rows to the file
# - add a ref_gene_id field with . or a gene or a comma separated gene list according to whether the transcript is supported by at least one ref tr
# - report the verbose info of the exon rows to the transcript rows and remove it from the exon rows
# - have gene_id, transcript_id and ref_gene_id in all rows
# - have a transcript id corresponding to a ref tr id whenever the tmerge transcript corresponds to a tr of the ref (taking the longest in case there are many)

# exemple of usage:
###################
# srun --mem=8G --pty bash
# dir=~/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/TAGADA.v1.0.1.GRCg6a.102.21-06-26/annotation/samples
# pgm1=~/fragencode/tools/multi/Scripts/make_gff_ok.awk
# pgm2=~/fragencode/tools/multi/Scripts/tmerge.buildLoci2okgff.awk
# pgm3=~/fragencode/tools/multi/Scripts/gff2gff.awk
# ref=~/fragencode/data/species/gallus_gallus/GRCg6a.102/gallus_gallus.gtf
# cd $dir
# time zcat alltr.supported_intron_chain.and.ref.sorted.tmerge.gene_id.gff.gz | awk -f $pgm1 | awk -v fileRef=$ref -f $pgm2 | awk -f $pgm3  | sort -k1,1 -k4,4n -k5,5n > alltr.supported_intron_chain.and.ref.sorted.tmerge.gene_id.complete.gff
# real	4m40.284s

# after $pgm1 we have this
# 1	tmerge	exon	9314	9364	0	-	.	gene_id "LOC_000000040528"; transcript_id "TM_000000000082"; contains "ref:ENSGALT00000098984"; contains_count "1"; 3p_dists_to_3p "0"; 5p_dists_to_5p "0"; flrpm "2.357812"; longest "ref:ENSGALT00000098984"; longest_FL_supporters "ref:ENSGALT00000098984"; longest_FL_supporters_count "1"; mature_RNA_length "897"; meta_3p_dists_to_5p "1"; meta_5p_dists_to_5p "0"; rpm "2.357812"; spliced "1"; 
# ...
# 1	tmerge	exon	5536	5952	0	-	.	gene_id "LOC_000000040528"; transcript_id "TM_000000000100"; contains "lung_stage3_rep2_female:STRG.1380.4,ileum_stage2_rep2_female:STRG.1729.4,lung_stage2_rep4_male:STRG.1861.3,muscle_stage1_rep4_female:STRG.1409.3,lung_stage1_rep1_male:STRG.1379.4,...,liver_stage3_rep3_female:STRG.1954.4,liver_stage2_rep3_male:STRG.1147.3,cerebellum_stage2_rep3_male:STRG.1505.4,kidney_stage1_rep1_male:STRG.1486.5,cerebellum_stage3_rep4_male:STRG.1211.2,skin_stage1_rep4_female:STRG.2041.3,cerebellum_stage2_rep2_female:STRG.1588.4,ileum_stage2_rep2_female:STRG.1729.4"; longest_FL_supporters_count "57"; mature_RNA_length "594"; meta_3p_dists_to_5p "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1"; meta_5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; rpm "134.395284"; spliced "1"; 

# alltr.supported_intron_chain.and.ref.sorted.tmerge.gene_id.complete.gff
# 1	tagada	exon	5273	5524	0	-	.	gene_id "LOC_000000040528"; transcript_id "ENSGALT00000101177"; ref_gene_id "ENSGALG00000054818";
# 1	tagada	exon	5273	5524	0	-	.	gene_id "LOC_000000040528"; transcript_id "ENSGALT00000101177"; ref_gene_id "ENSGALG00000054818";
# 23628 (12 fields)
# 961840 (14 fields)
# 42973 (28 fields)


BEGIN{
    # here we read the ref gene annotation gtf file to get the correspondance between transcript id and gene id
    while (getline < fileRef >0)
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
}

{
    # here we read each exon row from tmerge and buildLoci and keep the verbose info for the transcript rows
    # we also compute the transcript and gene boundaries to make transcript and gene rows
    ct="";
    found=0;
    reftrlg=0;
    refgnid="";
    reftrid="";
    split($0,a,"\t");
    split(a[9],b," ");
    
    # find out whether there is a ref transcript in the supporting transcripts and in this case change the transcript id to the longest of them
    k=1;
    while(found==0&&b[k]!="")
    {
	if(b[k]=="contains")
	{
	    found=1;
	    ct=b[k+1];
	}
	k+=2;
    }
    if(ct~/ref/)
    {
	split(ct,c,"\"");
	split(c[2],d,",");
	k=1;
	while(d[k]!="")
	{
	    split(d[k],e,":");
	    if(e[1]=="ref")
	    {
		refgnid=(refgnid)(refgn[e[2]])(",");
		if((reftrlg==0)||(lg[e[2]]>reftrlg))
		{
		    reftrlg=lg[e[2]];
		    reftrid=e[2];
		}
	    }
	    k++;
	}
	$12="\""reftrid"\"\;";
    }

    nbex[$10]++;
    nbex[$12]++;
    # remember verbose information for transcript rows
    if(nbex[$12]==1)
    {
	k=1;
	while(b[k]!="")
	{
	    if(b[k]!="gene_id"&&b[k]!="transcript_id")
	    {
		info[$12]=(info[$12])(b[k])(" ")(b[k+1])(" ");
	    }
	    k+=2;
	}
    }

    # print the current exon row
    s=(refgnid=="" ? "." : substr(refgnid,1,(length(refgnid)-1)));
    print $1, "tagada", $3, $4, $5, $6, $7, $8, $9"\t"$10"\t"$11"\t"$12"\tref_gene_id \""s"\";";
    # compute and update the characteristics of transcript and gene rows to be printed in the end part
    chr[$10]=$1;
    chr[$12]=$1;
    str[$10]=$7;
    str[$12]=$7;
    gene[$12]=$10;
    if((gbeg[$10]=="")||($4<gbeg[$10]))
    {
	gbeg[$10]=$4;
    }
    if((gend[$10]=="")||($5>gend[$10]))
    {
	gend[$10]=$5;
    }
    if((tbeg[$12]=="")||($4<tbeg[$12]))
    {
	tbeg[$12]=$4;
    }
    if((tend[$12]=="")||($5>tend[$12]))
    {
	tend[$12]=$5;
    }
    if(nbex[$10]==1)
    {
	info[$10]="gene_id "$10" ref_gene_id \""s"\";";
    }
    if(nbex[$12]==1)
    {
	info[$12]=info[$12]" ref_gene_id \""s"\";";
    }
    
}

END{
    for(tr in tbeg)
    {
	print chr[tr], "tagada", "transcript", tbeg[tr], tend[tr], 0, str[tr], ".", "gene_id "(gene[tr])" transcript_id "(tr)" "info[tr];
    }
    for(gn in gbeg)
    {
	print chr[gn], "tagada", "gene", gbeg[gn], gend[gn], 0, str[gn], ".", info[gn];
    }
}
