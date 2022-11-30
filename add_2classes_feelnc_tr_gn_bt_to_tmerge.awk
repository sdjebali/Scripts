# add_2classes_feelnc_tr_gn_bt_to_tmerge.awk
# The aim is quite similar to add_3classes_feelnc_tr_gn_bt_to_smerge.awk  but only takes the main tagada annotation file as input
# this file is usually produced by tmerge and already contains feelnc transcript biotype which is either mRNA or lncRNA

# The rules are therefore very simple: if the gene has at least one transcript that is considered protein_coding by the ref
# or mRNA by feelnc it will have the feelnc biotype mRNA, otherwise if it has at least one transcript labelled as lncRNA by
# the reference or feelnc then it will be labelled lncRNA, otherwise it will be labelled other

# It takes as input a single file:
##################################
# - the tagada novel gtf file given as input to feelnc as main input file (exon, transcript and gene rows, gene_id
#   as first key in the 9th field, transcript_id as second key for exons and transcripts)
# It outputs:
# - a single gtf file which is the same as the tagada novel annot file given as input except that it has
#   an additional feelnc_gene_biotype key at the end of the 9th field

# Example of usage
##################
# srun --mem=8G --pty bash
# cd /work/project/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/TAGADA.v2.1.0.GRCg6a.102.22-07-30/annotation
# pgm=~/fragencode/tools/multi/Scripts/add_2classes_feelnc_tr_gn_bt_to_tmerge.awk
# time awk -f $pgm novel.gtf > novel.with.gnbt.gtf
# real	0m16.596s

# input file novel.gtf 
# 2	tagada	gene	60230459	60323312	0	-	.	gene_id "LOC_000000004120"; ref_gene_id "."; feelnc_gene_biotype "mRNA";
# 2	tagada	transcript	60234249	60323312	0	-	.	gene_id "LOC_000000004120"; transcript_id "TM_000000407027"; contains "ileum_stage3.filtered:STRG.8622.1,muscle_stage2.filtered:STRG.7929.1,lung_stage1.filtered:STRG.8599.3,muscle_stage1.filtered:STRG.8172.4"; contains_count "4"; 3p_dists_to_3p "0,0,2,0"; 5p_dists_to_5p "0,5,31,43"; flrpm "4.35556"; longest "ileum_stage3.filtered:STRG.8622.1"; longest_FL_supporters "ileum_stage3.filtered:STRG.8622.1,lung_stage1.filtered:STRG.8599.3,muscle_stage1.filtered:STRG.8172.4,muscle_stage2.filtered:STRG.7929.1"; longest_FL_supporters_count "4"; mature_RNA_length "8548"; meta_3p_dists_to_5p "1,1,0.999766027140852,1"; meta_5p_dists_to_5p "0,0.000584932147870847,0.00362657931679925,0.00503041647168928"; rpm "4.35556"; spliced "1"; ref_gene_id "."; feelnc_gene_biotype "mRNA";
# 34712 (14 fields)
# 812243 (16 fields)
# 129578 (18 fields)
# 280908 (20 fields)
# 33399 (22 fields)
# 64717 (42 fields)
# 39140 (44 fields)
# 25713 (46 fields)
# 8702 (48 fields)  *** gene_id first in all rows, transcript_id second in exon and transcript rows
#                       there were only rows of 4, 10, 12, 14, 16 or 18 fields in the corresponding file of add_3classes_feelnc_tr_gn_bt_to_smerge.awk 


# output file novel.with.gnbt.gtf
# 2	tagada	gene	60230459	60323312	0	-	.	gene_id "LOC_000000004120"; ref_gene_id "."; feelnc_gene_biotype "mRNA";
# 2	tagada	transcript	60234249	60323312	0	-	.	gene_id "LOC_000000004120"; transcript_id "TM_000000407027"; contains "ileum_stage3.filtered:STRG.8622.1,muscle_stage2.filtered:STRG.7929.1,lung_stage1.filtered:STRG.8599.3,muscle_stage1.filtered:STRG.8172.4"; contains_count "4"; 3p_dists_to_3p "0,0,2,0"; 5p_dists_to_5p "0,5,31,43"; flrpm "4.35556"; longest "ileum_stage3.filtered:STRG.8622.1"; longest_FL_supporters "ileum_stage3.filtered:STRG.8622.1,lung_stage1.filtered:STRG.8599.3,muscle_stage1.filtered:STRG.8172.4,muscle_stage2.filtered:STRG.7929.1"; longest_FL_supporters_count "4"; mature_RNA_length "8548"; meta_3p_dists_to_5p "1,1,0.999766027140852,1"; meta_5p_dists_to_5p "0,0.000584932147870847,0.00362657931679925,0.00503041647168928"; rpm "4.35556"; spliced "1"; ref_gene_id "."; feelnc_gene_biotype "mRNA";
# 34712 (14 fields)
# 812243 (16 fields)
# 129578 (18 fields)
# 280908 (20 fields)
# 33399 (22 fields)
# 64717 (42 fields)
# 39140 (44 fields)
# 25713 (46 fields)
# 8702 (48 fields) 



# the novel annot gtf file has exon, transcript and gene rows and has gene_id as 1st key in the 9th field, and transcript_id
# as 2nd key for exon and transcript rows. It also has the information of feelnc transcript biotype in some exon and tr rows
# but not in all. here we remember all the exon, transcript and gene rows (exons and tr indexed by tr id, gene rows indexed
# by gn id, in hashtables), but we will output the complete file with feelnc gene biotype in the end part of the script
# note that here we also remember for each transcript its reference transcript biotype if it exists as well as its feelnc tr biotype
# if it exists, and this will serve to compute the feelnc gene biotype that will be reported for exon, transcript and gene rows at the end
{
    if($3=="exon")
    {
	nbex[$12]++;
	exrow[$12,nbex[$12]]=$0;
    }
    else
    {
	# when it reads a transcript row it tries to look for transcript_biotype (from the ref annot) and for feelnc_biotype (from feelnc) information
	# note that tr and gn indices include both double quotes and the semi colon
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
		k++;
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
	# this feelnc gene bt will go to the tr and its exons but also to the gene of the tr (see gnbt[b[2]]=gnbt[t] below)
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

