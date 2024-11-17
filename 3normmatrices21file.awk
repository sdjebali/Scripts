# 3normmatrices21file.awk
# the idea is to make a single tsv file with header for ggplot2 from a set of experiments present in 3 matrices
# one for tmm norm, one for loess norm and one for the raw expr values (as input)
# !!! be careful need to only print the expt,gene rows where the gene expr value exists in the 3 files !!!
# !!! since the normalization from NV's script could filter out some genes !!!

# srun --x11 --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/mixomics/Viguerie.Moro.Obesity/rnaseq.qpcr
# pgm=~/fragencode/tools/multi/Scripts/3normmatrices21file.awk
# time awk -v fileRef1=50_random_rnaseq_expt.txt -v fileRef2=normfeatures_normc_TMM.tsv -v fileRef3=normfeatures_normc_loess.tsv -f $pgm QB_GENE_rawCounts.okheader.tsv > 50_random_rnaseq_expt.id.cidno.gnid.normmeth.exprval.tsv


# fileRef1 = 50_random_rnaseq_expt.txt
# 31122_1
# 10662_1
# 50 (1 fields)

# fileRef2=normfeatures_normc_TMM.tsv 
# "X10011_1"	"X10012_3" ... "X91092_1"	"X91101_1"
# "ENSG00000000003"	104.046202062656	80.0779698983842 ... 100.268628528504	62.958319363479
# 1 (1004 fields) = header
# 54043 (1005 fields) = body

# fileRef3=normfeatures_normc_loess.tsv
# "X10011_1"	"X10012_3" ...
# "ENSG00000000003"	93.9980484082609	78.760602172892	 ...
# 1 (1004 fields) = header
# 54043 (1005 fields) = body

# main input = QB_GENE_rawCounts.okheader.tsv 
# gnid	10011_1	10012_3 ...
# ENSG00000000003	1784	1305	...
# 54044 (1005 fields) = header + body

# output = 50_random_rnaseq_expt.id.cidno.gnid.normmeth.exprval.tsv
# exptid	cidno	gnid	normmeth	exprval
# 31122_1	1	ENSG00000000003	raw	687
# 8106451 (5 fields) *** real	0m20.178s

BEGIN{
    OFS="\t";
    print "exptid", "cidno", "gnid", "normmeth", "exprval";
    while (getline < fileRef1 >0)
    {
	nbexp++;
	expt[nbexp]=$1;
	okexp[$1]=1;
    }
# fileRef2=normfeatures_normc_TMM.tsv 
# "X10011_1"	"X10012_3" ... "X91092_1"	"X91101_1"
# "ENSG00000000003"	104.046202062656	80.0779698983842 ... 100.268628528504	62.958319363479
# 1 (1004 fields) = header
# 54043 (1005 fields) = body
    while (getline < fileRef2 >0)
    {
	n++;
	gsub(/\"/,"",$0);
	if(n==1)
	{
	    for(i=1; i<=NF; i++)
	    {
		gsub(/X/,"",$0);
		if(okexp[$i]==1)
		{
		    j++;
		    idx[j]=(i+1);
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
# fileRef3=normfeatures_normc_loess.tsv
# "X10011_1"	"X10012_3" ...
# "ENSG00000000003"	93.9980484082609	78.760602172892	 ...
# 1 (1004 fields) = header
# 54043 (1005 fields) = body
    while (getline < fileRef3 >0)
    {
	n++;
	gsub(/\"/,"",$0);
	if(n==1)
	{
	    gsub(/X/,"",$0);
	    for(i=1; i<=NF; i++)
	    {
		if(okexp[$i]==1)
		{
		    j++;
		    idx[j]=(i+1);
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

# main input = QB_GENE_rawCounts.okheader.tsv 
# gnid	10011_1	10012_3 ...
# ENSG00000000003	1784	1305	...
# 54044 (1005 fields) = header + body
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
    gene[nbgn]=$1;
    for(k=1; k<=j; k++)
    {
	raw[$1,meta[k]]=$(idx[k]);
    }
}

END{
    for(j=1; j<=nbexp; j++)
    {
	e=expt[j];
	split(e,a,"_");
	for(i=1; i<=nbgn; i++)
	{
	    g=gene[i];
	    if(raw[g,e]!=""&&tmm[g,e]!=""&&loess[g,e]!="")
	    {
		print e, a[2], g, "raw", raw[g,e];
		print e, a[2], g, "tmm", tmm[g,e];
		print e, a[2], g, "loess", loess[g,e];
	    }
	}
    }
}
