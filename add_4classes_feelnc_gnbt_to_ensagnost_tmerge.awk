# add_4classes_feelnc_gnbt_to_ensagnost_tmerge.awk
# The aim here is simply to add a feelnc_gene_biotype to a gtf file produced by tagada (with tmerge) on a ref file that is not ensembl
# (for instance isoseq) and on illumina reads. The input file therefore already has a feelnc_biotype for transcripts from the 4 possible
# classes mRNA, lncRNA, TUCp and noORF, and the output file contains the same info but with a feelnc_gene_biotype only made from the
# transcript feelnc_biotype, and that is made very simply this way:
# - if one of the tr of the gene is mrna then the gene is mrna
# - otherwise if one of the tr of the gene is lncrna then the gene is lncrna
# - otherwise if one of the tr of the gene is noORF then the gene is noORF
# - otherwise the gene is TUCp


# It takes as input a single file:
##################################
# - the tagada novel gtf file given as input to feelnc as main input file (exon, transcript and gene rows, gene_id
#   as first key in the 9th field, transcript_id as second key for exons and transcripts)
# It outputs:
# - a single gtf file which is the same as the tagada novel annot file given as input except that it has
#   an additional feelnc_gene_biotype key at the end of the 9th field

# Example of usage
##################
# srun --x11 --mem=8G --pty bash
# dir=~/fragencode/workspace/geneswitch/results/rnaseq
# pgm=~/fragencode/tools/multi/Scripts/add_4classes_feelnc_gnbt_to_ensagnost_tmerge.awk
# sp=gallus_gallus
# cd $dir/$sp/nf-feelnc.tagiso.23-03-03/results/annotation
# time awk -f $pgm novel.clean.gtf > novel.clean.withgnbt.gtf 
# real	0m37.063s   *** checked on two genes with 1 tr and on two genes with 2 tr that it was working fine
# also as we see below we have the same number of rows but just 1 (key,value) pair added for each one

# input file novel.clean.gtf 
# 1	tagada	exon	5315	5512	0	-	.	gene_id "LOC_000000158014"; transcript_id "TM_000000000012"; ref_gene_id "."; ensembl_gene_id "."; feelnc_biotype "mRNA";
# 1	tagada	exon	5315	5524	0	-	.	gene_id "LOC_000000158014"; transcript_id "TM_000000000026"; ref_gene_id "."; ensembl_gene_id "."; feelnc_biotype "mRNA";
# 14764 (14 fields)
# 19168 (16 fields)
# 49409 (18 fields)
# 840371 (20 fields)
# 23230 (22 fields)
# 958338 (24 fields)
# 18584 (44 fields)
# 71047 (46 fields)
# 9706 (48 fields)
# 122339 (50 fields)*** gene_id first in all rows, transcript_id second in exon and transcript rows


# output file novel.clean.withgnbt.gtf 
# 14	tagada	gene	6815112	6843877	0	+	.	gene_id "LOC_000000004151"; ref_gene_id "G8864"; ensembl_gene_id "ENSGALG00000005767"; best_ref "ENSGALG00000005767"; feelnc_gene_biotype "mRNA";
# 14	tagada	transcript	6830505	6843874	0	+	.	gene_id "LOC_000000004151"; transcript_id "G8864.17"; contains "ref:G8864.17"; contains_count "1"; 3p_dists_to_3p "0"; 5p_dists_to_5p "0"; flrpm "2.398841"; longest "ref:G8864.17"; longest_FL_supporters "ref:G8864.17"; longest_FL_supporters_count "1"; mature_RNA_length "2413"; meta_3p_dists_to_5p "1"; meta_5p_dists_to_5p "0"; rpm "2.398841"; spliced "1"; ref_gene_id "G8864"; tmerge_tr_id "TM_000000587832"; transcript_biotype ""; ensembl_gene_id "ENSGALG00000005767"; best_ref "ENSGALG00000005767"; feelnc_biotype "mRNA"; feelnc_gene_biotype "mRNA";
# 14764 (16 fields)
# 19168 (18 fields)
# 49409 (20 fields)
# 840371 (22 fields)
# 23230 (24 fields)
# 958338 (26 fields)
# 18584 (46 fields)
# 71047 (48 fields)
# 9706 (50 fields)
# 122339 (52 fields)  


# the novel annot gtf file has exon, transcript and gene rows and has gene_id as 1st key in the 9th field, and transcript_id
# as 2nd key for exon and transcript rows. It also has the information of feelnc transcript biotype in some exon and tr rows
# but not in all. here we remember all the exon, transcript and gene rows (exons and tr indexed by tr id, gene rows indexed
# by gn id, in hashtables), but we will output the complete file with feelnc gene biotype in the end part of the script
# note that here we also remember for each transcript its feelnc tr biotype since this will serve to compute the feelnc gene
# biotype that will be reported for exon, transcript and gene rows at the end
{
    if($3=="exon")
    {
	nbex[$12]++;
	exrow[$12,nbex[$12]]=$0;
    }
    else
    {
	# when it reads a transcript row it looks for feelnc_biotype (from feelnc) information
	# note that tr and gn indices include both double quotes and the semi colon
	if($3=="transcript")
	{
	    found=0;
	    trrow[$12]=$0;
	    trlist[$10]=(trlist[$10])($12)(",");
	    split($0,a,"\t");
	    split(a[9],b," ");
	    k=1;
	    while((found==0)&&(b[k]!=""))
	    {
		if(b[k]=="feelnc_biotype")
		{
		    split(b[k+1],c,"\"");
		    trbt[$12]=c[2];
		    found=1;
		}
		k+=2;
	    }
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
	# this feelnc gene bt will go to the tr and its exons but also to the gene of the tr (see gnbt[b[2]]=gnbt[t] below)
	split(trlist[b[2]],c,",");
	k=1;
	while(c[k]!="")
	{
	    trbtlist[b[2]]=(trbtlist[b[2]])(trbt[c[k]])(",");
	    k++;
	}
	if(trbtlist[b[2]]~/mRNA/)
	{
	    gnbt[t]="mRNA";
	}
	else
	{
	    if(trbtlist[b[2]]~/lncRNA/)
	    {
		gnbt[t]="lncRNA";
	    }
	    else
	    {
		if(trbtlist[b[2]]~/noORF/)
		{
		    gnbt[t]="noORF";
		}
		else
		{
		    if(trbtlist[b[2]]~/TUCp/)
		    {
			gnbt[t]="TUCp";
		    }
		    else
		    {
			# but in theory we should not have other
			gnbt[t]="other";
		    }
		}
	    }
	}
	gnbt[b[2]]=gnbt[t];

        # contrary to add_3classes_feelnc_tr_gn_bt_to_smerge.awk but like add_2classes_feelnc_tr_gn_bt_to_tmerge.awk
	# only add the feelnc gene bt of the tr to all rows including gene rows
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

