# extend_pchic_to_5kb.awk

# as its name indicate takes as input a bedpe file of pchic (or other 3d) relations
# and extends each fragment to 5kb when not already above that

# example of usage
##################
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd /work/project/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018
# pgm=~/fragencode/tools/multi/Scripts/extend_pchic_to_5kb.awk
# time awk -f $pgm pchic.ipsccm.montefiori2018.hg19.distinct.sorted.overgencv49tss.part1.part2.reltype.bedpe | sort -V -k1,1 -k2,2n -k3,3n -k5,5n -k6,6n > pchic.ipsccm.montefiori2018.hg19.distinct.sorted.overgencv49tss.part1.part2.reltype.fragextto5kb.bedpe

# input file = pchic.ipsccm.montefiori2018.hg19.distinct.sorted.overgencv49tss.part1.part2.reltype.bedpe
# chr1	834807	896871	chr1:834807_835706:893041_896871:59699.5:NOC2L*NM:p-o:rev.order:0:1	6.6	.
# 378259 (6 fields)

# output file = pchic.ipsccm.montefiori2018.hg19.distinct.sorted.overgencv49tss.part1.part2.reltype.fragextto5kb.bedpe
# chr1	832757	837756	chr1	892456	897456	6.6	NOC2L*NM	rev.order	0	1	p-o
# 378259 (12 fields) *** real	0m11.438s

BEGIN{
    OFS="\t";
}

{
    if(($3-$2)<5000)
    {
	d1=(5000-($3-$2));
	toadd1=int(d1/2);
	$2=$2-toadd1;
	$3=$3+toadd1;
    }
    if(($6-$5)<5000)
    {
	d2=(5000-($6-$5));
	toadd2=int(d2/2);
	$5=$5-toadd2;
	$6=$6+toadd2;
    }
    print;
}
