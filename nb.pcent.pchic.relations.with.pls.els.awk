# nb.pcent.pchic.relations.with.pls.els.awk
# given a pchic bedpe file where the two fragments are given in gx order
# but with indication of init.order or rev.order in column 12 and p-p or p-o in column 11
# and with 6 last fields (columns 13 to 18) including booleans saying whether or not the
# 1st part (in gx order) and the 2nd part (in gx order) overlap ccre pls, els and dels
# respectively 
# provides
# 1. in a file called pchic.9catwrtplsels.nb.pcent.tsv, the number and pcent of all relations
#    falling in each of the 9 following categories:
#    - number and % of the relations with a pls on the prom side = ppls
#    - number and % of the relations with a pls on the other end side = opls
#    - number and % of the relations with an els on the prom side = pels
#    - number and % of the relations with an els on the other end side = oels
#    - number and % of the relations with a pls on the prom side and an els on the other end side = pplsoels
#    - number and % of the relations with a pls on the other end side and an els on the prom side = oplspels
#    - number and % of the relations with a pls on both sides = pplsopls
#    - number and % of the relations with an els on both sides = pelsoels
#    - number and % of the relations with neither a pls nor an els on any side = pnaona
# 2. as output the same file as the input but with 9 additional columns that are booleans about these 9 categories
# !!! here no use of the dels information for the moment !!!

# example
# cd ~/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018
# pgm=~/fragencode/tools/multi/Scripts/nb.pcent.pchic.relations.with.pls.els.awk
# time awk -f $pgm pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.bedpe > pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.9catwrtplsels.bedpe
# real	0m0.667s

# input = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	0
# chr1	834807	835706	chr1	894208	897150	6.58	KLHL17*NM	1	0	p-o	rev.order	0	0	0	1	1	0
# 401098 (18 fields)


# pchic.9catwrtplsels.nb.pcent.tsv
# ppls	401098	364021	90.7561
# opls	401098	55436	13.8211
# 9 (4 fields)

# main output = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.9catwrtplsels.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	0	1	0	1	0	0	0	00	0
# chr1	834807	835706	chr1	894208	897150	6.58	KLHL17*NM	1	0	p-o	rev.order	0	0	0	1	1	0	1	0	1	0	0	0	00	0
# 401098 (27 fields) 

BEGIN{
    OFS="\t";
}

{
    boolppls=0;
    boolopls=0;
    boolpels=0;
    booloels=0;
    boolpplsoels=0;
    booloplspels=0;
    boolpplsopls=0;
    boolpelsoels=0;
    boolpnaona=0;
    n++;
    # check first whether no overlap with pls and els on each side, in which case pnaona is incremented
    s=$13+$14+$16+$17;
    if(s==0)
    {
	boolpnaona=1;
	pnaona++;
    }
    else
    {
	# if in init order means the 1st frag is the prom frag and the 2nd frag is the other end frag
	if($12=="init.order")
	{
	    if($13==1)   # prom frag over pls
	    {
		ppls++;
		boolppls=1;
		if($17==1)  # otherend frag over els
		{
		    pplsoels++;
		    boolpplsoels=1;
		}
		if($16==1)  # otherend frag over pls
		{
		    pplsopls++;
		    boolpplsopls=1;
		}
	    }
	    if($16==1)  # otherend frag over pls
	    {
		opls++;
		boolopls=1;
		if($14==1)  # prom frag over els
		{
		    oplspels++;
		    booloplspels=1;
		}
	    }
	    if($14==1)   # prom frag over els
	    {
		pels++;
		boolpels=1;
		if($17==1)  # otherend frag over els
		{
		    pelsoels++;
		    boolpelsoels=1;
		}
	    }
	    if($17==1)   # otherend frag over els
	    {
		oels++;
		booloels=1;
	    }
	}
	# if in rev order means the 1st frag is the other end frag and the 2nd frag is the prom frag
	else
	{
	    if($12=="rev.order")
	    {
		if($16==1)   # prom frag over pls
		{
		    ppls++;
		    boolppls=1;
		    if($14==1)  # other end frag over els
		    {
			pplsoels++;
			boolpplsoels=1;
		    }
		    if($13==1)  # other end frag over pls
		    {
			pplsopls++;
			boolpplsopls=1;
		    }
		}
		if($13==1)  # other end frag over pls
		{
		    opls++;
		    boolopls=1;
		    if($17==1)  # prom frag over els
		    {
			oplspels++;
			booloplspels=1;
		    }
		}
		if($17==1)   # prom frag over els
		{
		    pels++;
		    boolpels=1;
		    if($14==1)  # other end frag over els
		    {
			pelsoels++;
			boolpelsoels=1;
		    }
		}
		if($14==1)  # other end frag over els
		{
		    oels++;
		    booloels=1;
		}
	    }
	}
    }
    print $0, boolppls, boolopls, boolpels, booloels, boolpplsoels, booloplspels, boolpplsopls, boolpelsoels, boolpnaona;
}

# - number and % of the relations with a pls on the prom side = ppls
# - number and % of the relations with a pls on the other end side = opls
# - number and % of the relations with an els on the prom side = pels
# - number and % of the relations with an els on the other end side = oels
# - number and % of the relations with a pls on the prom side and an els on the other end side = pplsoels
# - number and % of the relations with a pls on the other end side and an els on the prom side = oplspels
# - number and % of the relations with a pls on both sides = pplsopls
# - number and % of the relations with an els on both sides = pelsoels
# - number and % of the relations with neither a pls nor an els on any side = pnaona
END{
    print "ppls", n, ppls, ppls/n*100 > "pchic.9catwrtplsels.nb.pcent.tsv";
    print "opls", n, opls, opls/n*100 > "pchic.9catwrtplsels.nb.pcent.tsv";
    print "pels", n, pels, pels/n*100 > "pchic.9catwrtplsels.nb.pcent.tsv";
    print "oels", n, oels, oels/n*100 > "pchic.9catwrtplsels.nb.pcent.tsv";
    print "pplsoels", n, pplsoels, pplsoels/n*100 > "pchic.9catwrtplsels.nb.pcent.tsv";
    print "oplspels", n, oplspels, oplspels/n*100 > "pchic.9catwrtplsels.nb.pcent.tsv";
    print "pplsopls", n, pplsopls, pplsopls/n*100 > "pchic.9catwrtplsels.nb.pcent.tsv";
    print "pelsoels", n, pelsoels, pelsoels/n*100 > "pchic.9catwrtplsels.nb.pcent.tsv";
    print "pnaona", n, pnaona, pnaona/n*100 > "pchic.9catwrtplsels.nb.pcent.tsv";
} 
