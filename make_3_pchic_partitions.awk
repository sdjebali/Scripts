# make_3_pchic_partitions.awk

# Given a pchic bedpe file with 32 fields of which the last 16 correspond to categories of overlap of the two fragments
# with ccre pls, els and dels, provides the same file but with 3 additional columns for 3 partitions of the pchic relations:
# - one based on the prom frag overlap with pls, els and dels
# - one based on the other end frag overlap with pls, els and dels
# - one based on the two frag overlap with pls, els and dels

# In the input file we have the last 16 columns, from col 19 to 34 that represent the following information
# - col 19 (a) = pp (prom frag over pls) / col 21 (b) = pe (prom frag over els) / col 23 (c) = pde (prom frag over dels)
# - col 20 (d) = op (otherend frag over pls) / col 22 (e) = oe (otherend frag over els) / col 24 (f) = ode (otherend frag over dels)
# - col 25 (g) = ppoe (prom over pls and other over els) / col 26 (h) = ppode (...) / col 27 (i) = peop /
#   col 28 (j) = pdeop / col 29 (k) = ppop / col 30 (l) = peoe /
#   col 31 (m) = pdeode / col 32 (n) = pona / col 33 (o) = pnao / col 34 (p) = pnaona

# The 3 paritions will be made using these 3 priority orders
############################################################
# - Pchic relation partitioning based on promoter fragments:
#   pp (a) > pde (c) > pe (b) > pna   (4 classes, pde before pe otherwise no pde class)
# - Pchic relation partitioning based on other end fragments:
#   ode (f) > oe (e) > op (d) > ona  (4 classes)
# - Pchic relation partitioning based on prom and other end fragments:
#   ppode (h) > ppop (k) > ppoe (g) > pdeode (m) > peoe (l) > pdeop (j) > peop (i) > pona (n) > pnao (o) > pnaona (p)   (10 classes)

# Example of usage:
###################
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018
# pgm=~/fragencode/tools/multi/Scripts/make_3_pchic_partitions.awk
# time awk -f $pgm pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.bedpe > pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.3partitions.bedpe

# input file = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	01	0	1	0	0	0	0	0	0	0	0	0	0	0
# 401098 (32 fields)

# output file = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.3partitions.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	0	1	0	1	0	0	0	00	0	0	0	0	0	1	0	0	pp	ona	pona
# 401098 (37 fields)


BEGIN{
    OFS="\t";
}

{
    print $0, ppart($19,$21,$23), opart($20,$22,$24), popart($25,$26,$27,$28,$29,$30,$31,$32,$33,$34);
}

function ppart(a,b,c)
{
    return (a==1 ? "pp" : (c==1 ? "pde" : (b==1 ? "pe" : "pna")));
}

function opart(d,e,f)
{
    return (f==1 ? "ode" : (e==1 ? "oe" : (d==1 ? "op" : "ona")));
}

function popart(g,h,i,j,k,l,m,n,o,p)
{
    return (h==1 ? "ppode" : (k==1 ? "ppop" : (g==1 ? "ppoe" : (m==1 ? "pdeode" : (l==1 ? "peoe" : (j==1 ? "pdeop" : (i==1 ? "peop" : (n==1 ? "pona" : (o==1 ? "pnao" : (p==1 ? "pnaona" : "NA"))))))))));
}
