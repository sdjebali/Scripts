# add_mean_sd_cv_to_gnexprmatrix.awk
# this scripts takes as input a matrix of gene expression values with a header (with n or n-1 elements compared to the body)
# that has genes as rows and samples as columns, and with gene id in the 1st column, and adds 3 columns to this file
# that are for a given gene:
# - mean expression across samples
# - expression standard deviation across samples
# - expression coefficient of variation across samples
# in case we want it on log10+1e-3 we have to add -v applog=1 to the command, otherwise it does not do it
# TODO next:
# - make it possible to end at column c2 for expression values (in case we have something more than expr values in the last columns)
# Note that I checked the below formula gave exactly the same SD as doing sqrt(sum((xi-mean)^2)/(N-1)) 
   

# example
# srun --x11 --mem=8G --pty bash
# cd ~/fragencode/workspace/geneswitch/results/multi/rnaseq.srnaseq/gallus_gallus/PCIT
# pgm=~/fragencode/tools/multi/Scripts/add_mean_sd_cv_to_gnexprmatrix.awk
# time awk -f $pgm tagiso.and.shortstack2exp.normexpr.tsv > tmp
# awk 'BEGIN{OFS="\t"} {print (NR==1 ? "gnid" : $1), $(NF-2), $(NF-1), $NF}' tmp > tagiso.and.shortstack2exp.mean.sd.cv.onnormexpr.tsv
# rm tmp

# input tagiso.and.shortstack2exp.normexpr.tsv
# cerebellum_stage1_rep1	cerebellum_stage1_rep2	...	skin_stage3_rep3	skin_stage3_rep4
# LOC_000000000000	52.5538746178964	57.7705037477386	...	14.3229316708071	14.0720588638608
# 1 (84 fields)
# 153304 (85 fields)

# output tmp
# cerebellum_stage1_rep1	cerebellum_stage1_rep2	 ...	skin_stage3_rep3	skin_stage3_rep4	mean	sd	cv
# LOC_000000000000	52.5538746178964	57.7705037477386	...	14.3229316708071	14.0720588638608	16.1693	12.8104	0.792265
# 1 (87 fields)
# 153304 (88 fields)

BEGIN{
    if(c1=="")
    {
	c1=2;
    }
}

NR==1{
    OFS="\t";
    print $0, "mean", "sd", "cv";
}

NR>=2{
    s=0;
    s2=0;
    for(i=c1; i<=NF; i++)
    {
	if(applog!=1)
	{
	    v=$i;
	}
	else
	{
	    v=log10($i);
	}
	s+=v;
	s2+=v*v;
    }
    tot=(NF-c1+1);
    me=s/tot;
    sd=sqrt((tot*s2-s*s)/((tot*(tot-1))));
    print $0, me, sd, sd/me;
}

function log10(x)
{
    return (log(x+0.001)/log(10))
}

