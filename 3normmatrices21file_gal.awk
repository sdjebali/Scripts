# 3normmatrices21file_gal.awk
# the idea is to make a single tsv file with header for ggplot2 from a set of experiments present in 3 matrices
# one for tmm norm, one for loess norm and one for the raw expr values (as input)
# !!! be careful need to only print the expt,gene rows where the gene expr value exists in the 3 files !!!
# !!! since the normalization from NV's script could filter out some genes !!!
# !!! this is a generalisation of 3normmatrices21file.awk that was made for diogenes data specifically !!!
# !!! but also here the header and the body of the matrices have the same number of columns, unlike in diogenes for some reason !!!
# !!! and here we have additional columns in the main input at the beg so need a from param to say where to start the raw expr !!!
# !!! but in fact not since we only take the columns where the samples we want are !!!
# !!! and the gene id is not in the 1st column of the main input so need another param to say where the gene id is, called colgnid !!!
# !!! also since the sample id can be very long, here the 1st fileRef file will have two columns for the samples we want to include in the plot !!!
# !!! the sample id present in the matrices and a small sample id for the plot !!!

# srun --x11 --mem=32G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/analysis/srnaseq/normalisation/gallus_gallus
# indir=~/fragencode/workspace/geneswitch/results/srnaseq/gallus_gallus/shortstack_quantification-4.0.4.GRCg6a.25-02-19
# pgm=~/fragencode/tools/multi/Scripts/3normmatrices21file_gal.awk
# time awk -v colgnid=2 -v fileRef1=labExpId.long.short.tsv -v fileRef2=$indir/normalization/sRNA_features_normc_TMM.tsv -v fileRef3=$indir/normalization/sRNA_features_normc_loess.tsv -f $pgm $indir/counts_PAD_20_Dicer_20_30_mincov_0.1.tsv > lid.gnid.normmeth.exprval.tsv
# real	1m14.656s

# fileRef1=labExpId.long.short.tsv 
# chicken_stage1_embryoday8_cerebellum_rep1_1_trimmed	cerebellum_1_1
# 84 (2 fields)

# fileRef2=$indir/normalization/sRNA_features_normc_TMM.tsv
# ""	"chicken_stage1_embryoday8_cerebellum_rep1_1_trimmed"	"chicken_stage1_embryoday8_cerebellum_rep2_1_trimmed"	...	"chicken_stage3_newbornday1_skin_rep3_1_trimmed"	"chicken_stage3_newbornday1_skin_rep4_1_trimmed"
# "ShortStack_merge:Unknown_sRNA_locus_237061"	159.432674935645	131.850541976771	...	22.2986361493367	33.2293960965707
# 119376 (85 fields) 

# fileRef3=$indir/normalization/sRNA_features_normc_loess.tsv 
# ""	"chicken_stage1_embryoday8_cerebellum_rep1_1_trimmed"	"chicken_stage1_embryoday8_cerebellum_rep2_1_trimmed"	...	"chicken_stage3_newbornday1_skin_rep3_1_trimmed"	"chicken_stage3_newbornday1_skin_rep4_1_trimmed"
# "ShortStack_merge:Unknown_sRNA_locus_237061"	409.661213130473	82.1172022187131 ...	14.3107191267776	19.6942895001063
# 119376 (85 fields)

# main input = $indir/counts_PAD_20_Dicer_20_30_mincov_0.1.tsv
# Coords	Name	MIRNA	chicken_stage1_embryoday8_cerebellum_rep1_1_trimmed	chicken_stage1_embryoday8_cerebellum_rep2_1_trimmed	...	chicken_stage3_newbornday1_skin_rep3_1_trimmed	chicken_stage3_newbornday1_skin_rep4_1_trimmed
# MT:1-978	ShortStack_merge:Unknown_sRNA_locus_237061	N	14209	11413	...	1571	2483
# 119376 (87 fields)

# output = lid.gnid.normmeth.exprval.tsv
# lid	gnid	normmeth	exprval
# cerebellum_1_2	ShortStack_merge:Unknown_sRNA_locus_237061	raw	11413
# 29724376 (4 fields)

BEGIN{
    if(colgnid=="")
    {
	colgnid=1;
    }
    OFS="\t";
    print "lid", "gnid", "normmeth", "exprval";

# fileRef1=labExpId.long.short.tsv 
# chicken_stage1_embryoday8_cerebellum_rep1_1_trimmed	cerebellum_1_1
# 84 (2 fields)
    while (getline < fileRef1 >0)
    {
	nbexp++;
	expt[nbexp]=$1;
	corr[$1]=$2;
	okexp[$1]=1;
    }
    
# fileRef2=$indir/normalization/sRNA_features_normc_TMM.tsv
# ""	"chicken_stage1_embryoday8_cerebellum_rep1_1_trimmed"	"chicken_stage1_embryoday8_cerebellum_rep2_1_trimmed"	...	"chicken_stage3_newbornday1_skin_rep3_1_trimmed"	"chicken_stage3_newbornday1_skin_rep4_1_trimmed"
# "ShortStack_merge:Unknown_sRNA_locus_237061"	159.432674935645	131.850541976771	...	22.2986361493367	33.2293960965707
# 119376 (85 fields) 
   while (getline < fileRef2 >0)
    {
	n++;
	gsub(/\"/,"",$0);
	if(n==1)
	{
	    for(i=1; i<=NF; i++)
	    {
		if(okexp[$i]==1)
		{
		    j++;
		    idx[j]=i+1;
		    meta[j]=$i;
		}
	    }
	}
	else
	{
	    for(k=1; k<=j; k++)
	    {
		okgn[$1]=1;
		tmm[$1,meta[k]]=$(idx[k]);
	    }
	}
    }
    n=0;
    j=0;
    
# fileRef3=$indir/normalization/sRNA_features_normc_loess.tsv 
# ""	"chicken_stage1_embryoday8_cerebellum_rep1_1_trimmed"	"chicken_stage1_embryoday8_cerebellum_rep2_1_trimmed"	...	"chicken_stage3_newbornday1_skin_rep3_1_trimmed"	"chicken_stage3_newbornday1_skin_rep4_1_trimmed"
# "ShortStack_merge:Unknown_sRNA_locus_237061"	409.661213130473	82.1172022187131 ...	14.3107191267776	19.6942895001063
# 119376 (85 fields)
   while (getline < fileRef3 >0)
    {
	n++;
	gsub(/\"/,"",$0);
	if(n==1)
	{
	    for(i=1; i<=NF; i++)
	    {
		if(okexp[$i]==1)
		{
		    j++;
		    idx[j]=i+1;
		    meta[j]=$i;
		}
	    }
	}
	else
	{
	    for(k=1; k<=j; k++)
	    {
		okgn[$1]=1;
		loess[$1,meta[k]]=$(idx[k]);
	    }
	}
    }
}

# main input = $indir/counts_PAD_20_Dicer_20_30_mincov_0.1.tsv
# Coords	Name	MIRNA	chicken_stage1_embryoday8_cerebellum_rep1_1_trimmed	chicken_stage1_embryoday8_cerebellum_rep2_1_trimmed	...	chicken_stage3_newbornday1_skin_rep3_1_trimmed	chicken_stage3_newbornday1_skin_rep4_1_trimmed
# MT:1-978	ShortStack_merge:Unknown_sRNA_locus_237061	N	14209	11413	...	1571	2483
# 119376 (87 fields)
NR==1{
    j=0;
    for(i=2; i<=NF; i++)
    {
	if(okexp[$i]==1)
	{
	    j++;
	    idx[j]=i;
	    meta[j]=$i;
	}
    }
}

NR>=2{
    nbgn++;
    gene[nbgn]=$colgnid;
    for(k=1; k<=j; k++)
    {
	raw[$colgnid,meta[k]]=$(idx[k]);
    }
}

END{
    for(j=1; j<=nbexp; j++)
    {
	e=expt[j];
	for(i=1; i<=nbgn; i++)
	{
	    g=gene[i];
	    if(raw[g,e]!=""&&tmm[g,e]!=""&&loess[g,e]!="")
	    {
		print corr[e], g, "raw", raw[g,e];
		print corr[e], g, "tmm", tmm[g,e];
		print corr[e], g, "loess", loess[g,e];
	    }
	}
    }
}
