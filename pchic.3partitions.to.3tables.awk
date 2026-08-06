# pchic.3partitions.to.3tables.awk
# given a bedpe file of pchic relations and the number of a column in this file representing a partition
# (based on promfrag or on otherend frag or on both frag), outputs a tsv file of 4 fields without header
# that has for each class of the partition
# - the total nb of relations
# - the class
# - the number of relations in this class
# - the corresponding %

# Example of usage
##################
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018
# pgm=~/fragencode/tools/multi/Scripts/pchic.3partitions.to.3tables.awk
# col=37
# awk -v col=$col -f $pgm pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.3partitions.bedpe | sort -k2,2 

# input file = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.3partitions.bedpe 
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	01	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	pp	ona	pona
# 401098 (37 fields)

# output 
# 401098  01_ppode   121850  30.3791
# 401098  02_ppop    55733   13.8951
# 401098  03_ppoe    29162   7.27054
# 401098  04_pdeode  1928    0.480681
# 401098  05_peoe    10097   2.51734
# 401098  06_pdeop   40      0.00997263
# 401098  07_peop    190     0.04737
# 401098  08_pona    167370  41.728
# 401098  09_pnao    7032    1.75319
# 401098  10_pnaona  7696    1.91873
# 10 (4 fields)

BEGIN{
    cl["pp"]="1_pp";
    cl["pde"]="2_pde";
    cl["pe"]="3_pe";
    cl["pna"]="4_pna";
    cl["ode"]="1_ode";
    cl["oe"]="2_oe";
    cl["op"]="3_op";
    cl["ona"]="4_ona";
    cl["ppode"]="01_ppode";
    cl["ppop"]="02_ppop";
    cl["ppoe"]="03_ppoe";
    cl["pdeode"]="04_pdeode";
    cl["peoe"]="05_peoe";
    cl["pdeop"]="06_pdeop";
    cl["peop"]="07_peop";
    cl["pona"]="08_pona";
    cl["pnao"]="09_pnao";
    cl["pnaona"]="10_pnaona";
}

{
    tot++;
    tot1[cl[$col]]++;
}

END{
    OFS="\t";
    for(c in tot1)
    {
	print tot, c, tot1[c], tot1[c]/tot*100;
    }
}
