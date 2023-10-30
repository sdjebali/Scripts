# make_tr_gn_file_with_bt_and_nov.awk
# !!! supposes that gnid and transcript id come first in the 9th field (no tr id for gene rows)

# example
# srun --x11 --mem=8G --pty bash
# cd ~/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/TAGADA.v2.1.0.GRCg6a.102.22-07-30/annotation
# time awk -f $pgm novel.with.gnbt.gnnov.gtf > novel.tr.gn.id.bt.nov.tsv
# real	0m19.251s

# input file
# 2	tagada	transcript	60234249	60323312	0	-	.	gene_id "LOC_000000004120"; transcript_id "TM_000000407027"; contains "ileum_stage3.filtered:STRG.8622.1,muscle_stage2.filtered:STRG.7929.1,lung_stage1.filtered:STRG.8599.3,muscle_stage1.filtered:STRG.8172.4"; contains_count "4"; 3p_dists_to_3p "0,0,2,0"; 5p_dists_to_5p "0,5,31,43"; flrpm "4.35556"; longest "ileum_stage3.filtered:STRG.8622.1"; longest_FL_supporters "ileum_stage3.filtered:STRG.8622.1,lung_stage1.filtered:STRG.8599.3,muscle_stage1.filtered:STRG.8172.4,muscle_stage2.filtered:STRG.7929.1"; longest_FL_supporters_count "4"; mature_RNA_length "8548"; meta_3p_dists_to_5p "1,1,0.999766027140852,1"; meta_5p_dists_to_5p "0,0.000584932147870847,0.00362657931679925,0.00503041647168928"; rpm "4.35556"; spliced "1"; ref_gene_id "."; feelnc_gene_biotype "mRNA"; novelty "unknown";
# 2	tagada	exon	60234249	60241926	0	-	.	gene_id "LOC_000000004120"; transcript_id "TM_000000407027"; ref_gene_id "."; feelnc_gene_biotype "mRNA"; novelty "unknown";
# 34712 (16 fields)
# 812243 (18 fields)
# 129578 (20 fields)
# 280908 (22 fields)
# 33399 (24 fields)
# 64717 (44 fields)
# 39140 (46 fields)
# 25713 (48 fields)
# 8702 (50 fields)

# output file
# TM_000000692181	LOC_000000004035	other	mRNA	unknown	unknown
# TM_000000957964	LOC_000000022007	lncRNA	lncRNA	unknown	mix
# 138272 (6 fields)

BEGIN{
    OFS="\t";
}

{
    split($0,a,"\t");
    gnid="";
    trid="";
    split(a[9],b," ");
    k=1;
    if(a[3]=="transcript")
    {
	while(b[k]!="")
	{
	    split(b[k+1],c,"\"");
	    if(b[k]=="gene_id")
	    {
		gnid=c[2];
	    }
	    else
	    {
		if(b[k]=="transcript_id")
		{
		    trid=c[2];
		    gn[trid]=gnid;
		}
		else
		{
		    if(b[k]=="feelnc_gene_biotype")
		    {
			gnbt[gnid]=c[2];
		    }
		    else
		    {
			if(b[k]=="transcript_biotype")
			{
			    enstrbt[trid]=c[2];
			}
			else
			{
			    if(b[k]=="feelnc_biotype")
			    {
				feeltrbt[trid]=c[2];
			    }
			    else
			    {
				if(b[k]=="novelty")
				{
				    novtr[trid]=c[2];
				}
			    }
			}
		    }
		}
	    }
	    k++;
	}
    }
    else
    {
	if(a[3]=="gene")
	{
	    while(b[k]!="")
	    {
		split(b[k+1],c,"\"");
		if(b[k]=="gene_id")
		{
		    gnid=c[2];
		}
		else
		{
		    if(b[k]=="novelty")
		    {
			novgn[gnid]=c[2];
		    }
		}
		k++;
	    }
	}
    }
}

END{
    for(t in gn)
    {
	# here in case a transcript is said to be coding by feelnc but not by ensembl it will be coding
	# and in case a transcript is said to be coding by ensembl but not by feelnc it will be coding
	trbt=((enstrbt[t]=="protein_coding"||feeltrbt[t]=="mRNA") ? "mRNA" : ((enstrbt[t]=="lncRNA"||feeltrbt[t]=="lncRNA") ? "lncRNA" : "other"));
	print t, gn[t], trbt, gnbt[gn[t]], novtr[t], novgn[gn[t]];
    }
}
