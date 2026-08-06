# pchic.or.hichip.to.sc.dist.opttype.partition.class.awk

# Given as input
################
# - 4 variables passed with -v corresponding for a given partition of pchic or hichip relations of:
###################################################################################################
#   * the partition name (with variable partname)
#   * the column in the input bedpe file of pchic relations where the classes of the partition are (with variable partcol)
#   * the comma separated list of ordered class (with a number before that is not in the bedpe file) of the partition (with variable partclasslist)
#   * the column where the score is in the input bed file
#   * the column where the type is in the input bed file in case we want it
# - a bedpe file of pchic or hichip relations with indication of:
#################################################################
# - score in the given column
# - p-p or p-o type in the given column and in case we want it 
# - class of the relation according to the given partition in column provided
# Provides as output a tsv file with header to be used for ggplot2 next that has
################################################################################
# - score
# - distance
# - optional type
# - partition name
# - class in this partition

# Example of usage
##################
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018
# pgm=~/fragencode/tools/multi/Scripts/pchic.or.hichip.to.sc.dist.opttype.partition.class.awk
# partname=twofrag.part
# partcol=37
# partclasslist=01_ppode,02_ppop,03_ppoe,04_pdeode,05_peoe,06_pdeop,07_peop,08_pona,09_pnao,10_pnaona
# time awk -v partname=$partname -v partcol=$partcol -v partclasslist=$partclasslist -v scorecol=7 -v typecol=11 -f $pgm pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.fragextto5kb.distinct.sorted.overpls.els.dels.first.second.part.16catwrtplsels.3partitions.bedpe > tmp

# input file = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.fragextto5kb.distinct.sorted.overpls.els.dels.first.second.part.16catwrtplsels.3partitions.bedpe
# chr1	832757	837756	chr1	892456	897456	6.6	NOC2L*NM	1	0	p-o	rev.order	0	1	1	1	1	01	0	1	1	0	1	1	1	0	0	0	1	0	0	0	0	pp	ode	ppode
# 378259 (37 fields)

# output file = tmp
# score	dist	type	partition	class
# 6.6	59699	p-o	twofrag.part	01_ppode
# 378260 (5 fields)  *** 1+nb relations since each relation is in a single class with a partition

BEGIN{
    OFS="\t";
    # make the correspondence between class in the bedpe file and class provided with a number before (for ggplot2)
    split(partclasslist,a,",");
    k=1;
    while(a[k]!="")
    {
	split(a[k],b,"_");
	corr[b[2]]=a[k];
	k++;
    }
    # in case the column of the type of relation is provided we output it, otherwise not
    if(typecol!="")
    {
	print "score", "dist", "type", "partition", "class";
    }
    else
    {
	print "score", "dist", "partition", "class";
    }

}

{
    if($1==$4)
    {
	d=dist($2,$3,$5,$6);
    }
    else
    {
	d="NA";
    }
    # in case the column of the type of relation is provided we output it, otherwise not
    if(typecol!="")
    {
	print $scorecol, d, $typecol, partname, corr[$partcol];
    }
    else
    {
	print $scorecol, d, partname, corr[$partcol];
    }
}

# computes distance between two fragments of (x,y) and (z,t) coordinates on the same chromosome
# as middle to middle
function dist(x,y,z,t){
    m1=(x+y)/2;
    m2=(z+t)/2;
    return ((m2-m1) > 0 ? int(m2-m1) : int(m1-m2));
}
