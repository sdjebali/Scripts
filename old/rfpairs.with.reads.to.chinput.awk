# rfpairs.with.reads.to.chinput.awk
# from a 3 column tsv file without header of numerically ordered restriction fragment (rf) pairs (num ids) with number of supporting reads
# and two auxiliary files in bed kind of format with bait rf coordinates and rf coordinates, makes a chinput file for chicago

# takes as input
# - as fileRef1 a bait file in bed kind of format (in fact chicago bait file)
# - as fileRef2 a rf file in bed kind of format (in fact chicago rmap file) (super set of bait file)
# - a 3 column tsv file without header of numerically ordered restriction fragment (rf) pairs with number of supporting reads
# provides as output a chinput file for chicago
# Note1: in case none of the two rf is a bait does not output the row
# Note2: can take quite some ram (10G or so) because of rf file that has to be read first and coord information stored !!!

# example
# cd ~/fragencode/workspace/sdjebali/geneswitch/pipelines/phic/chicago
# bait=~/fragencode/workspace/catchi/data/annotation/restriction_fragments_withprobes.bed
# rf=~/fragencode/workspace/catchi/data/annotation/restriction_fragments.bed
# pgm=~/fragencode/tools/multi/Scripts/rfpairs.with.reads.to.chinput.awk
# echo "
#      cd ~/fragencode/workspace/sdjebali/geneswitch/pipelines/phic/chicago
#      awk -v fileRef1=$bait -v fileRef2=$rf -f $pgm Capture_HiC_liver_1_11_S2_all.alphaordered.rfpairs.with.nbvalidpairs.tsv > Capture_HiC_liver_1_11_S2_allvalidpairs.chinput
# " | awk 'BEGIN{print "\#\!\/bin\/sh"} {print}' > Capture_HiC_liver_1_11_S2_all.to.chinput.sh
# sbatch --mem=12G --cpus-per-task=1 -J tochinput --mail-user=sarah.djebali@inserm.fr --mail-type=END,FAIL --workdir=$PWD --export=ALL -p workq Capture_HiC_liver_1_11_S2_all.to.chinput.sh
# SLURM Job_id=22612443 Name=tochinput Ended, Run time 00:02:15, COMPLETED, ExitCode 0   *** and used 10G of ram

# inputs
# 1	21124	21492	HIC_1_126	0	+
# 16981 (6 fields)
# 1	0	330	HIC_1_1	0	+
# 15 188 239 (6 fields)
# HIC_10_10000	HIC_5_628875	2
# 28 080 654 (3 fields)

# output
# baitID	otherEndID	N	otherEndLen	distSign
# HIC_10_9862	HIC_10_10001	3	205	21841
# 3323046 (5 fields)


BEGIN{
    OFS="\t";
    print "baitID", "otherEndID", "N", "otherEndLen", "distSign";
    while (getline < fileRef1 >0)
    {
	bait[$4]=1;
    }
    while (getline < fileRef2 >0)
    {
	chr[$4]=$1;
	gbeg[$4]=$2;
	gend[$4]=$3;
    }
}

{
    if(bait[$1]==1)
    {
	first=$1;
	second=$2;
	len=gend[$2]-gbeg[$2];
	if(chr[$1]==chr[$2])
	{
	    m1=(gbeg[$1]+gend[$1])/2;
	    m2=(gbeg[$2]+gend[$2])/2;
	    dist=m2-m1;
	}
	else
	{
	    dist="NA";
	}
    }
    else
    {
	if(bait[$2]==1)
	{
	    first=$2;
	    second=$1;
	    len=gend[$1]-gbeg[$1];
	    if(chr[$1]==chr[$2])
	    {
		m1=(gbeg[$2]+gend[$2])/2;
		m2=(gbeg[$1]+gend[$1])/2;
		dist=m2-m1;
	    }
	    else
	    {
		dist="NA";
	    }
	}
    }
    if(bait[$1]==1||bait[$2]==1)
    {
	print first, second, $3, len, dist;
    }
}
