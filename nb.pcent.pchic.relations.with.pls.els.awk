# nb.pcent.pchic.relations.with.pls.els.awk
# given a pchic bedpe file where the two fragments are given in gx order
# but with indication of init.order or rev.order in column 12 and p-p or p-o in column 11
# and with 6 last fields (columns 13 to 18) including booleans saying whether or not the
# 1st part (in gx order) and the 2nd part (in gx order) overlap ccre pls, els and dels
# (distal els) respectively, will provide:
# 1. in a file called pchic.14catwrtplsels.nb.pcent.tsv, the number and pcent of all relations
#    falling in each of the 14 following categories:
#    - number and % of the relations with a pls on the prom side = ppls ($19)
#    - number and % of the relations with a pls on the other end side = opls ($20)
#    - number and % of the relations with an els on the prom side = pels ($21)
#    - number and % of the relations with an els on the other end side = oels ($22)
#    - number and % of the relations with a dels on the other end side = pdels ($23)
#    - number and % of the relations with a dels on the other end side = odels ($24)
#    - number and % of the relations with a pls on the prom side and an els on the other end side = pplsoels ($25)
#    - number and % of the relations with a pls on the prom side and a dels on the other end side = pplsodels ($26)
#    - number and % of the relations with an els on the prom side and a pls on the other end side = pelsopls ($27)
#    - number and % of the relations with a dels on the prom side and a pls on the other end side = pdelsopls ($28)
#    - number and % of the relations with a pls on both sides = pplsopls ($29)
#    - number and % of the relations with an els on both sides = pelsoels ($30)
#    - number and % of the relations with a dels on both sides = pdelsodels ($31)
#    - number and % of the relations with neither a pls nor an els on any side = pnaona ($32)
# 2. as output the same file as the input but with 14 additional columns that are booleans about these 14 categories
# !!! Note that we have 1/20 and 1/31 of the prom and other end frag respectively that overlap both a pels and a dels !!!
# !!! and since this is not very high we will not specifically ask that when over dels it does not overlap pels !!!

# example
# cd ~/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018
# pgm=~/fragencode/tools/multi/Scripts/nb.pcent.pchic.relations.with.pls.els.awk
# time awk -f $pgm pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.bedpe > pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.14catwrtplsels.bedpe
# real	0m0.843s

# input = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	0
# chr1	834807	835706	chr1	894208	897150	6.58	KLHL17*NM	1	0	p-o	rev.order	0	0	0	1	1	0
# 401098 (18 fields)

# pchic.14catwrtplsels.nb.pcent.tsv
# ppls        401098  364021  90.7561 
# opls        401098  63220   15.7617 
# pels        401098  372990  92.9922 
# oels        401098  222243  55.4087 
# pdels       401098  24064   5.99953  
# odels       401098  133551  33.2964 
# pplsoels    401098  203362  50.7013 
# pplsodels   401098  121850  30.3791 
# pelsopls    401098  58984   14.7056 
# pdelsopls   401098  3332    0.83072 
# pplsopls    401098  58306   14.5366 
# pelsoels    401098  208149  51.8948 
# pdelsodels  401098  8748    2.18101 
# pnaona      401098  7696    1.91873 
# 14 (4 fields)                       

# main output = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.14catwrtplsels.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	0	1	0	1	0	0	0	00	0	0	0	0	0	0
# chr1	834807	835706	chr1	894208	897150	6.58	KLHL17*NM	1	0	p-o	rev.order	0	0	0	1	1	0	1	0	1	0	0	0	00	0	0	0	0	0	0
# 401098 (32 fields) 

BEGIN{
    OFS="\t";
}

{
    boolppls=0;
    boolopls=0;
    boolpels=0;
    boolpdels=0;
    booloels=0;
    boolodels=0;
    boolpplsoels=0;
    boolpplsodels=0;
    boolpelsopls=0;
    boolpdelsopls=0;
    boolpplsopls=0;
    boolpelsoels=0;
    boolpdelsodels=0;
    boolpnaona=0;
    n++;

    # $13=ppls / $14=pels / $15=pdels / $16=opls / $17=oels / $18=odels 
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
		if($18==1)  # otherend frag over dels
		{
		   pplsodels++;
		   boolpplsodels=1; 
		}
		if($16==1)  # otherend frag over pls
		{
		    pplsopls++;
		    boolpplsopls=1;
		}
	    }
	    # $13=ppls / $14=pels / $15=pdels / $16=opls / $17=oels / $18=odels 
	    if($16==1)  # otherend frag over pls
	    {
		opls++;
		boolopls=1;
		if($14==1)  # prom frag over els
		{
		    pelsopls++;
		    boolpelsopls=1;
		}
		if($15==1)  # prom frag over dels
		{
		    pdelsopls++;
		    boolpdelsopls=1;
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
	    if($15==1)   # prom frag over dels
	    {
		pdels++;
		boolpdels=1;
		if($18==1)  # otherend frag over dels
		{
		    pdelsodels++;
		    boolpdelsodels=1;
		}
	    }
	    if($17==1)   # otherend frag over els
	    {
		oels++;
		booloels=1;
	    }
	    if($18==1)   # otherend frag over dels
	    {
		odels++;
		boolodels=1;
	    }
	}
	# $13=ppls / $14=pels / $15=pdels / $16=opls / $17=oels / $18=odels 
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
		    if($15==1)  # other end frag over dels
		    {
			pplsodels++;
			boolpplsodels=1;
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
			pelsopls++;
			boolpelsopls=1;
		    }
		    if($18==1)  # prom frag over dels
		    {
			pdelsopls++;
			boolpdelsopls=1;
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
		if($18==1)   # prom frag over dels
		{
		    pdels++;
		    boolpdels=1;
		    if($15==1)  # other end frag over dels
		    {
			pdelsodels++;
			boolpdelsodels=1;
		    }
		}
		if($14==1)  # other end frag over els
		{
		    oels++;
		    booloels=1;
		}
		if($15==1)  # other end frag over dels
		{
		    odels++;
		    boolodels=1;
		}
	    }
	}
    }
    print $0, boolppls, boolopls, boolpels, booloels, boolpdels, boolodels, boolpplsoels, boolpplsodels, boolpelsopls, boolpdelsopls, boolpplsopls, boolpelsoels, boolpdelsodels, boolpnaona;
}

# - number and % of the relations with a pls on the prom side = ppls ($19)
# - number and % of the relations with a pls on the other end side = opls ($20)
# - number and % of the relations with an els on the prom side = pels ($21)
# - number and % of the relations with an els on the other end side = oels ($22)
# - number and % of the relations with a dels on the other end side = pdels ($23)
# - number and % of the relations with a dels on the other end side = odels ($24)
# - number and % of the relations with a pls on the prom side and an els on the other end side = pplsoels ($25)
# - number and % of the relations with a pls on the prom side and a dels on the other end side = pplsodels ($26)
# - number and % of the relations with an els on the prom side and a pls on the other end side = pelsopls ($27)
# - number and % of the relations with a dels on the prom side and a pls on the other end side = pdelsopls ($28)
# - number and % of the relations with a pls on both sides = pplsopls ($29)
# - number and % of the relations with an els on both sides = pelsoels ($30)
# - number and % of the relations with a dels on both sides = pdelsodels ($31)
# - number and % of the relations with neither a pls nor an els on any side = pnaona ($32)
END{
    print "ppls", n, ppls, ppls/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "opls", n, opls, opls/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pels", n, pels, pels/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "oels", n, oels, oels/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pdels", n, pdels, pdels/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "odels", n, odels, odels/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pplsoels", n, pplsoels, pplsoels/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pplsodels", n, pplsodels, pplsodels/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pelsopls", n, pelsopls, pelsopls/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pdelsopls", n, pdelsopls, pdelsopls/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pplsopls", n, pplsopls, pplsopls/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pelsoels", n, pelsoels, pelsoels/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pdelsodels", n, pdelsodels, pdelsodels/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
    print "pnaona", n, pnaona, pnaona/n*100 > "pchic.14catwrtplsels.nb.pcent.tsv";
} 
