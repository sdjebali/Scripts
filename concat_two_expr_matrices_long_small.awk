# concat_two_expr_matrices_long_small.awk

# srun --x11 --mem=8G --pty bash
# cd /work/project/fragencode/workspace/geneswitch/results/multi/rnaseq.srnaseq/gallus_gallus
# tagiso=~/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/TAGADA.v2.1.1.GRCg6a.isoseq.23-01-13/quantification/normalization/normfeatures_normc_loess.tsv
# mirna=~/fragencode/workspace/geneswitch/results/srnaseq/gallus_gallus/mmquant-1.0.6.GRCg6a.novannot.23-09-07/normalization/normfeatures_normc_TMM.tsv
# pgm=~/fragencode/tools/multi/Scripts/concat_two_expr_matrices_long_small.awk
# time awk -v fileRef=$mirna -f $pgm $tagiso > tagiso.and.mirdeep2.normexpr.tsv

# fileRef=$mirna 
# "chicken_stage1_embryoday8_cerebellum_rep1_1"	"chicken_stage1_embryoday8_cerebellum_rep2_1"	...	"chicken_stage3_newbornday1_skin_rep3_1"	"chicken_stage3_newbornday1_skin_rep4_1"
# "mirdeep2.1"	0	0	...	0	0
# 1 (84 fields)
# 6725 (85 fields)

# main input = $tagiso
# "cerebellum_stage1_rep1"	"cerebellum_stage1_rep2"	...	"skin_stage3_rep3"	"skin_stage3_rep4"
# "LOC_000000000000"	52.5538746178964	57.7705037477386	...	14.3229316708071	14.0720588638608
# 1 (84 fields)
# 33929 (85 fields)

# output
# cerebellum_stage1_rep1	cerebellum_stage1_rep2	...	skin_stage3_rep3	skin_stage3_rep4
# LOC_000000000000	52.5538746178964	57.7705037477386	...	14.3229316708071	14.0720588638608
# 1 (84 fields)
# 40654 (85 fields)

BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	gsub(/\"/,"",$0);
	n++;
	if(n==1)
	{
	    for(i=1; i<=NF; i++)
	    {
		split($i,a,"_");
		expt[i]=a[4]"_"a[2]"_"a[5];
	    }
	}
	else
	{
	    for(i=2; i<=NF; i++)
	    {
		ok[$1]=1;
		expr[$1,expt[i-1]]=$i;
	    }
	}
    }
}

{
    gsub(/\"/,"",$0);
    m++;
    if(m==1)
    {
	tot=NF;
	for(i=1; i<=NF; i++)
	{
	    experiment[i]=$i;
	}
    }
    print;
}


END{
    for(small in ok)
    {
	s=(small)("\t");
	for(i=1; i<tot; i++)
	{
	    s=(s)(expr[small,experiment[i]])("\t");
	}
	print (s)(expr[small,experiment[i]]);
    }
}
