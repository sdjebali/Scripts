# gather_3_norm_expr_on_50_expt.awk
# The aim here is to produce for the diogenes project that has quantseq expr in the adipose tissue
# a file to give to ggplot2 to represent boxplots of gene expression for 50 experiments
# in 3 flavors: with raw counts, with tmm norm and with loess norm
# So here given:
################
# - a 1 column file without header providing the ids of 50 experiments (-v fileRef1)
# - 3 tsv files with headers with genes in rows (possibly with double quotes) and experiments in columns (possibly with double quotes)
#   and corresponding to gene expression matrices with the same nb of genes in the 2 first files and possible less than in the 3rd one:
#   * -v fileRef2 = tmm normalised expression
#   * -v fileRef3 = loess normalised expression
#   * main input file = raw non normalised expression
# produces a tsv file with header that has for each of the 50 experiments and each of the genes of the 3 matrices
# - experiment id
# - cid no
# - gene id
# - kind of norm from 3 (raw counts, tmm, loess)
# - expr value

# Example
# srun --x11 --mem=64G--time=10:00:00 --pty bash
# cd ~/work/mixomics/Viguerie.Moro.Obesity/quantseq
# pgm=~/fragencode/tools/multi/Scripts/gather_3_norm_expr_on_50_expt.awk
# time awk -v fileRef1=50_random_quantseq_expt.txt -v fileRef2=normfeatures_normc_TMM.tsv -v fileRef3=normfeatures_normc_loess.tsv -f $pgm cleaned_count.patientsok.cid23.geneswithexpr.headerforDE.tsv > 50_random_quantseq_expt.id.cidno.gnid.normmeth.exprval.tsv
# real	0m8.333s

# fileRef1=50_random_quantseq_expt.txt
# cid2_rep133
# 50 (1 fields)

# fileRef2=normfeatures_normc_TMM.tsv
# ""	"cid2_rep1"	"cid2_rep2"	...	"cid2_rep307"	"cid3_rep223"
# "ENSG00000000003"	37.058146377209	46.1916045737977	...	35.9552424813924	47.9796880958686
# 25769 (531 fields)

# fileRef3=normfeatures_normc_loess.tsv
# ""	"cid2_rep1"	"cid2_rep2"	...	"cid2_rep307"	"cid3_rep223"
# "ENSG00000000003"	37.1051697117774	45.6429271524703	...	34.7596937053737	46.6403401366107
# 25769 (531 fields)

# main input file = cleaned_count.patientsok.cid23.geneswithexpr.headerforDE.tsv
# gnid	cid2_rep1	cid2_rep2	...	cid2_rep307	cid3_rep223
# ENSG00000000003	195	217	...	159	180
# 32042 (531 fields)

# output file = 50_random_quantseq_expt.id.cidno.gnid.normmeth.exprval.tsv
# exptid	cidno	gnid	normmeth	exprval
# cid3_rep3	3	ENSG00000211642	1.raw	0
# 3865201 (5 fields)   *** 166M file


# in the begin reads the 3 fileRef file and store information
BEGIN{
    OFS="\t";
    # flag the 50 expts of interest
    # cid2_rep133
    # 50 (1 fields)
    while (getline < fileRef1 >0)
    {
	ok[$1]=1;
	ok["\""$1"\""]=1;
    }
    
    # for each of the experiments of interest, remember their index and the name of their experiment for the tmm matrix
    # when reading its header, and then the tmm expressions when reading its body and just for the expt of interest
    # ""	"cid2_rep1"	"cid2_rep2"	...	"cid2_rep307"	"cid3_rep223"
    # "ENSG00000000003"	37.058146377209	46.1916045737977	...	35.9552424813924	47.9796880958686
    # 25769 (531 fields)
    while (getline < fileRef2 >0)
    {
	n++;
	if(n==1)
	{
	    for(i=2; i<=NF; i++)
	    {
		if(ok[$i]==1)
		{
		    j++;
		    idx1[j]=i;
		    expt1[j]=$i;
		}
	    }
	}
	else
	{
	    for(m=1; m<=j; m++)
	    {
		ok1[$1]=1;
		tmm[$1,expt1[m]]=$(idx1[m]);
	    }
	}
    }

    n=0;
    # for each of the experiments of interest, remember their index and the name of their experiment for the loess matrix
    # when reading its header, and then the loess expressions when reading its body and just for the expt of interest
    while (getline < fileRef3 >0)
    {
	n++;
	if(n==1)
	{
	    for(i=2; i<=NF; i++)
	    {
		if(ok[$i]==1)
		{
		    k++;
		    idx2[k]=i;
		    expt2[k]=$i;
		}
	    }
	}
	else
	{
	    for(m=1; m<=k; m++)
	    {
		ok2[$1]=1;
		loess[$1,expt2[m]]=$(idx2[m]);
	    }
	}
    }
}

# when reading the main input file of raw count, do the same as for the tmm and loess matrix
# gnid	cid2_rep1	cid2_rep2	...	cid2_rep307	cid3_rep223
# ENSG00000000003	195	217	...	159	180
# 32042 (531 fields)
NR==1{
    for(i=2; i<=NF; i++)
    {
	if(ok[$i]==1)
	{
	    l++;
	    idx3[l]=i;
	    expt3[l]=$i;
	}
    }
}

# since the 2 norm expr matrices have a gene id with double quotes, add some to the gene id in the raw matrix
# and do the same for the expt
NR>=2{
    for(m=1; m<=l; m++)
    {
	gid="\""$1"\"";
	raw[gid,"\""expt3[m]"\""]=$(idx3[m]);
    }
}

# After having read the 4 files, write what we want = for each gene and each of the 50 exp write the raw, tmm and loess expr of the gene in this expt
END{
    print "exptid", "cidno", "gnid", "normmeth", "exprval";
    for(g in ok1)
    {
	if(ok2[g]==1)
	{
	    split(g,d,"\"");
	    gnid=d[2];
	    for(m=1; m<=j; m++)
	    {
		e=expt1[m];
		split(e,a,"\"");
		split(a[2],b,"_");
		split(b[1],c,"cid");
		cidno=c[2];
		print a[2], cidno, gnid, "1.raw", raw[g,e];
		print a[2], cidno, gnid, "2.tmm", tmm[g,e];
		print a[2], cidno, gnid, "3.loess", loess[g,e];
	    }
	}
    }
}
