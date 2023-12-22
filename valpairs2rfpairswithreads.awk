# from a hicpro valid pair file makes a file of numerically ordered restrictions fragments with numeric ids provided in fileRef

# example
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/pipelines/phic/chicago/pig.liver_1_11_S2
# valpairs=/home/sfoissac/fragencode/workspace/catchi/results/capturehic1/hic_results/data/Capture_HiC_liver_1_11_S2_all.allValidPairs
# pgm=~/fragencode/tools/multi/Scripts/valpairs2rfpairswithreads.awk
# time awk -v fileRef=pig.liver_1_11_S2.rf.id.no.tsv -f $pgm $valpairs | sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' > Capture_HiC_liver_1_11_S2_all.numordered.rfpairs.with.nbvalidpairs.tsv
# real	3m27.308s

# fileRef
# HIC_1_1	1
# HIC_1_2	2
# 15188239 (2 fields)
# $valpairs
# NS500446:722:HGFT7AFX2:1:11101:1019:7264    X    109329233    -    X    109334928    +    199    HIC_X_660414    HIC_X_660448    42    23
# NS500446:722:HGFT7AFX2:1:11101:1019:7499    16    22642992    +    16    22643184    -    131    HIC_16_132117    HIC_16_132119    23    42
# 34 338 589 (12 fields)


# output before the sort
# 14581415	14581449

# output Capture_HiC_liver_1_11_S2_all.numordered.rfpairs.with.nbvalidpairs.tsv
# 1	135	1
# 1	2840	2
# real	3m27.308s


BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	id[$1]=$2;
    }
}

{
    i=id[$9];
    j=id[$10];
    if(i<=j)
    {
	print i, j;
    }
    else
    {
	print j, i;
    }
}
