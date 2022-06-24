# add_2classes_feelnc_tr_gn_bt_to_tmerge.awk
# The aim is quite similar to add_3classes_feelnc_tr_gn_bt_to_smerge.awk  but only takes the main tagada annotation file as input
# this files is usually produced by tmerge and it already contains feelnc transcript biotype which is either mRNA or lncRNA

# The rules are simple: if the gene has at least one transcript that is considered protein_coding by the ref
# or mRNA by feelnc it will have the feelnc biotype mRNA, otherwise if it has at least one transcript
# labelled as lncRNA by the reference of feelnc then it will be labelled lncRNA, otherwise it will be labelled other


# It takes as input a single file:
##################################
# - the tagada novel gtf file given as input to feelnc as main input file (exon, transcript and gene rows, gene_id
#   as first key in the 9th field, transcript_id as second key for exons and transcripts)
# It outputs a single gtf file which is the same as the tagada novel annot file given as input except that it has
# an additional feelnc_gene_biotype key at the end of the 9th field

# example of usage
# srun --mem=8G --pty bash
# cd /home/cgenet/genrocwork/Carine/MONOPOLY/Bovin/TAGADA_test11/RESULTATS/annotation
# pgm=~/fragencode/tools/multi/Scripts/add_2classes_feelnc_tr_gn_bt_to_tmerge.awk
# outdir=~/fragencode/workspace/sdjebali/monopoly
# time awk -f $pgm novel.gtf > $outdir/novel.with.feelnc.biotypes.gtf


# input file novel.gtf 
# 1	tagada	exon	339070	339312	0	-	.	gene_id "LOC_000000051391"; transcript_id "ENSBTAT00000007786"; ref_gene_id "ENSBTAG00000006648"; tmerge_tr_id "TM_000000000015"; transcript_biotype "protein_coding";
# 1	tagada	exon	339070	339312	0	-	.	gene_id "LOC_000000051391"; transcript_id "ENSBTAT00000008737"; ref_gene_id "ENSBTAG00000006648"; tmerge_tr_id "TM_000000000016"; transcript_biotype "protein_coding";
# 29298 (12 fields)
# 393776 (14 fields)
# 18696 (16 fields)
# 391600 (18 fields)
# 7300 (20 fields)
# 34824 (40 fields)
# 6376 (42 fields)
# 38735 (44 fields)
# 2352 (46 fields)  *** gene_id first in all rows, transcript_id second in exon and transcript rows
#                       there were only rows of 4, 10, 12, 14, 16 or 18 fields in the corresponding file of add_3classes_feelnc_tr_gn_bt_to_smerge.awk 


# the novel annot gtf file has exon, transcript and gene rows and has gene_id as 1st key in the 9th field, and transcript_id
# as 2nd key for exon and transcript rows. It also has the information of feelnc transcript biotype in some exon and tr rows
# here we remember all the exon, transcript and gene rows (exons and tr indexed by tr id, gene rows indexed by gene id in hashtables)
# but we will output the complete file with feelnc gene biotype in the end part of the script
# note that here we also remember for each transcript its reference transcript biotype if it exists as well as its feelnc tr biotype
# if it exists and this will serve to compute the feelnc gene biotype that will be reported for exon, tr and gene rows afterwards
{
    if($3=="exon")
    {
	nbex[$12]++;
	exrow[$12,nbex[$12]]=$0;
    }
    else
    {
	# when it reads a transcript row it tries to look for transcript_biotype (from the ref annot) and for feelnc_biotype (from feelnc) information
	if($3=="transcript")
	{
	    trrow[$12]=$0;
	    trlist[$10]=(trlist[$10])($12)(",");
	    split($0,a,"\t");
	    split(a[9],b," ");
	    k=1;
	    while(b[k]!="")
	    {
		if(b[k]=="transcript_biotype")
		{
		    if(b[k+1]~/protein_coding/)
		    {
			mrna1[$12]=1;
		    }
		    else
		    {
			if(b[k+1]~/lncRNA/)
			{
			    lnc1[$12]=1;
			}
		    }
		}
		else
		{
		    if(b[k]=="feelnc_biotype")
		    {
			if(b[k+1]~/mRNA/)
			{
			    mrna2[$12]=1;
			}
			else
			{
			    if(b[k+1]~/lncRNA/)
			    {
				lnc2[$12]=1;
			    }
			}
		    }
		}
		k+=2;
	    }
	    trbt[$12]=(((mrna1[$12]==1)||(mrna2[$12]==1)) ? "mRNA" : (((lnc1[$12]==1)||(lnc2[$12]==1)) ? "lncRNA" : "other"));
	}
	else
	{
	    if($3=="gene")
	    {
		gnrow[$10]=$0;
	    }
	}
    }
}

END{
    for(t in trrow)
    {
	split(trrow[t],a,"\t");
	split(a[9],b," ");
	# count the nb of times we see the gene in order to only print the gene row once later on
	seen[b[2]]++;
	
	# look for the tr list of the gene of the current transcript t
	# and make the tr bt list out of it in order to compute the gn bt afterwards
	# this feelnc gene bt will go to the tr and its exons but also to the gene of the tr (see gnbt[b[2]]=gnbt[t])
	split(trlist[b[2]],c,",");
	k=1;
	while(c[k]!="")
	{
	    trbtlist[b[2]]=(trbtlist[b[2]])(trbt[c[k]])(",");
	    k++;
	}
	gnbt[t]=((trbtlist[b[2]]~/mRNA/) ? "mRNA" : ((trbtlist[b[2]]~/lncRNA/) ? "lncRNA" : "other"));
	gnbt[b[2]]=gnbt[t];

        # contrary to add_3classes_feelnc_tr_gn_bt_to_smerge.awk only add the feelnc gene bt of the tr to all rows including gene rows (not present before)
	if(seen[b[2]]==1)
	{
	    print gnrow[b[2]]" feelnc_gene_biotype \""(gnbt[b[2]])"\"\;";
	}
	print trrow[t]" feelnc_gene_biotype \""(gnbt[t])"\"\;";
	for(i=1; i<=nbex[t]; i++)
	{
	    print exrow[t,i]" feelnc_gene_biotype \""(gnbt[t])"\"\;";
	}
    }

}

