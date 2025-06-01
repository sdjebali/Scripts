# compute_euclidean_dist.awk
# As its name indicates, this script takes as input a tsv file with header that has indication of at least two vectors in two columns
# (by default 1 and 2) and computes the euclidean distance between these vectors (same kind of input as AB's R script Rscripts/ggpoints.R
# that traces the two vectors and computes pearson and spearman correlations). In case we want the columns with the two vectors to be
# different from 1 and 2 we can specify them with -v col1=x -v col2=y
# and in case we want it on the log10+0.001 of the values we just need to do -v applog=1

# example 
# srun --x11 --mem=8G --pty bash
# cd ~/fragencode/workspace/geneswitch/results/rnaseq/multi/projection.to.human/geneswitch/multi/examples
# pgm=~/fragencode/tools/multi/Scripts/compute_euclidean_dist.awk
# time awk -v col1=2 -v col2=3 -f $pgm lid.LOC_000000000017.LOC_000000002021LOESSon1to1divgnlg.avgacrossbiorep.tissue.stage.tsv
# real	0m0.016s


# input lid.LOC_000000000017.LOC_000000002021LOESSon1to1divgnlg.avgacrossbiorep.tissue.stage.tsv
# lid	gn1expr	gn2expr	tissue	stage
# cerebellum_stage1	1799.63	1652.87	cerebellum	stage1
# 22 (5 fields)

# output
# 6766.43

BEGIN{
    OFS="\t";
    if(col1=="")
    {
	col1=1;
    }
    if(col2=="")
    {
	col2=2;
    }
}

NR>=2{
    if(applog!=1)
    {
	s+=($col1-$col2)^2;
    }
    else
    {
	s+=(log10($col1)-log10($col2))^2;
    }
}

END{
    print sqrt(s);
}

function log10(x)
{
    return (log(x+0.001)/log(10))
}
