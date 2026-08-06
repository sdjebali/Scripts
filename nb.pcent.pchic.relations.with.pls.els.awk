# nb.pcent.pchic.relations.with.pls.els.awk
# given a pchic bedpe file where the two fragments are given in gx order
# but with indication of init.order or rev.order in column 12 and p-p or p-o in column 11
# and with 6 last fields (columns 13 to 18) including booleans saying whether or not the
# 1st part (in gx order) and the 2nd part (in gx order) overlap ccre pls, els and dels
# (distal els) respectively, will provide:
# 1. in a file called pchic.16catwrtplsels.nb.pcent.tsv, the number and pcent of all relations
#    falling in each of the 16 following categories:
#    - number and % of the relations with a pls on the prom side = pp ($19)
#    - number and % of the relations with a pls on the other end side = op ($20)
#    - number and % of the relations with an els on the prom side = pe ($21)
#    - number and % of the relations with an els on the other end side = oe ($22)
#    - number and % of the relations with a dels on the other end side = pde ($23)
#    - number and % of the relations with a dels on the other end side = ode ($24)
#    - number and % of the relations with a pls on the prom side and an els on the other end side = ppoe ($25)
#    - number and % of the relations with a pls on the prom side and a dels on the other end side = ppode ($26)
#    - number and % of the relations with an els on the prom side and a pls on the other end side = peop ($27)
#    - number and % of the relations with a dels on the prom side and a pls on the other end side = pdeop ($28)
#    - number and % of the relations with a pls on both sides = ppop ($29)
#    - number and % of the relations with an els on both sides = peoe ($30)
#    - number and % of the relations with a dels on both sides = pdeode ($31)
#    - number and % of the relations with an overlap on the prom side but no overlap on the other end side = pona ($32)
#    - number and % of the relations with an overlap on the other end side but no overlap on the prom side = pnao ($33)
#    - number and % of the relations with neither a pls nor an els on any side = pnaona ($34)
# 2. as output the same file as the input but with 16 additional columns that are booleans about these 16 categories
# !!! Note that we have 1/20 and 1/31 of the prom and other end frag respectively that overlap both a pels and a dels !!!
# !!! and since this is not very high we will not specifically ask that when over dels it does not overlap pels !!!

# example
# cd ~/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018
# pgm=~/fragencode/tools/multi/Scripts/nb.pcent.pchic.relations.with.pls.els.awk
# time awk -f $pgm pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.bedpe > pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.bedpe
# real	0m0.843s

# input = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	0
# chr1	834807	835706	chr1	894208	897150	6.58	KLHL17*NM	1	0	p-o	rev.order	0	0	0	1	1	0
# 401098 (18 fields)

# pchic.16catwrtplsels.nb.pcent.tsv
# pp      401098  364021  90.7561
# op      401098  63220   15.7617
# pe      401098  372990  92.9922
# oe      401098  222243  55.4087
# pde     401098  24064   5.99953
# ode     401098  133551  33.2964
# ppoe    401098  203362  50.7013
# ppode   401098  121850  30.3791
# peop    401098  117968  29.4113
# pdeop   401098  6664    1.66144
# ppop    401098  58306   14.5366
# peoe    401098  208149  51.8948
# pdeode  401098  8748    2.18101
# pona    401098  167370  41.728
# pnao    401098  7032    1.75319
# pnaona  401098  7696    1.91873
# 16 (4 fields)                       

# main output = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	0	1	0	1	0	0	0	00	0	0	0	0	0	1	0	0
# chr1	834807	835706	chr1	894208	897150	6.58	KLHL17*NM	1	0	p-o	rev.order	0	0	0	1	1	0	1	0	1	0	0	0	00	0	0	0	0	0	1	0	0
# 401098 (34 fields)

BEGIN{
    OFS="\t";
}

{
    boolpp=0;
    boolop=0;
    boolpe=0;
    boolpde=0;
    booloe=0;
    boolode=0;
    boolppoe=0;
    boolppode=0;
    boolpeop=0;
    boolpdeop=0;
    boolppop=0;
    boolpeoe=0;
    boolpdeode=0;
    boolpona=0;
    boolpnao=0;
    boolpnaona=0;
    n++;

    # $13=pp / $14=pe / $15=pde / $16=op / $17=oe / $18=ode 
    # check first whether no overlap with pls and els on each side, in which case pnaona is incremented
    # $15 and $18 are not used because dels is included in els and when s is 0, $15+$18 is also 0 (checked)
    s=$13+$14+$16+$17;
    if(s==0)
    {
	pnaona++;
	boolpnaona=1;
    }
    else
    {
	# if in init order means the 1st frag is the prom frag and the 2nd frag is the other end frag
	if($12=="init.order")
	{
	    if($13==1)   # prom frag over pls
	    {
		pp++;
		boolpp=1;
		if($17==1)  # otherend frag over els
		{
		    ppoe++;
		    boolppoe=1;
		}
		if($18==1)  # otherend frag over dels
		{
		   ppode++;
		   boolppode=1; 
		}
		if($16==1)  # otherend frag over pls
		{
		    ppop++;
		    boolppop=1;
		}
	    }
	    # $13=pp / $14=pe / $15=pde / $16=op / $17=oe / $18=ode
	    if($14==1)   # prom frag over els
	    {
		pe++;
		boolpe=1;
		if($16==1)   # otherend frag over pls
		{
		    peop++;
		    boolpeop=1;
		}
		if($17==1)   # otherend frag over els
		{
		    peoe++;
		    boolpeoe=1;
		}
	    }
	    if($15==1)    # prom frag over dels
	    {
		pde++;
		boolpde=1;
		if($16==1)  # otherend frag over pls
		{
		    pdeop++;
		    boolpdeop=1;
		}
		if($18==1)   # otherend frag over dels
		{
		    pdeode++;
		    boolpdeode=1;
		}
	    }
	    if($16==1)   # otherend frag over pls
	    {
		op++;
		boolop=1;
		if($14==1)   # prom frag over els
		{
		    peop++;
		    boolpeop=1;
		}
		if($15==1)  # prom frag over dels
		{
		    pdeop++;
		    boolpdeop=1;
		}
	    }
	    if($17==1)   # otherend frag over els
	    {
		oe++;
		booloe=1;
	    }
	    if($18==1)  # otherend frag over dels
	    {
		ode++;
		boolode=1;
	    }
	}
	# $13=pp / $14=pe / $15=pde / $16=op / $17=oe / $18=ode 
	# if in rev order means the 1st frag is the other end frag and the 2nd frag is the prom frag
	else
	{
	    if($12=="rev.order")
	    {
		if($16==1)   # prom frag over pls
		{
		    pp++;
		    boolpp=1;
		    if($14==1)  # other end frag over els
		    {
			ppoe++;
			boolppoe=1;
		    }
		    if($15==1)  # other end frag over dels
		    {
			ppode++;
			boolppode=1;
		    }
		    if($13==1)  # other end frag over pls
		    {
			ppop++;
			boolppop=1;
		    }
		}
		if($17==1)   # prom frag over els
		{
		    pe++;
		    boolpe=1;
		    if($13==1)   # otherend frag over pls
		    {
			peop++;
			boolpeop=1;
		    }
		    if($14==1)   # otherend frag over els
		    {
			peoe++;
			boolpeoe=1;
		    }
		}
		if($18==1)    # prom frag over dels
		{
		    pde++;
		    boolpde=1;
		    if($13==1)  # otherend frag over pls
		    {
			pdeop++;
			boolpdeop=1;
		    }
		    if($15==1)   # otherend frag over dels
		    {
		    pdeode++;
		    boolpdeode=1;
		    }
		}
		if($13==1)   # otherend frag over pls
		{
		    op++;
		    boolop=1;
		    if($17==1)   # prom frag over els
		    {
			peop++;
			boolpeop=1;
		    }
		    if($18==1)  # prom frag over dels
		    {
			pdeop++;
			boolpdeop=1;
		    }
		}
		if($14==1)   # otherend frag over els
		{
		    oe++;
		    booloe=1;
		}
		if($15==1)  # otherend frag over dels
		{
		    ode++;
		    boolode=1;
		}
	    }
	}
	# at this point we could still have relations where the sum of all 2 frag cat is 0 while not being pnaona
	# those are cases where either the 1st part or the 2nd part overlaps a ccre but not the other, so two
	# categories are added for them, pona in case it is the 2nd part that does not overlap anything and pnao
	# in case it is the 1st part that does not overlap anything
	s2=boolppoe+boolppode+boolpeop+boolpdeop+boolppop+boolpeoe+boolpdeode;
	s3=boolpp+boolpe;
	s4=boolop+booloe;
	if(s2==0)
	{
	    if(s3!=0&&s4==0)
	    {
		pona++;
		boolpona=1;
	    }
	    else
	    {
		if(s3==0&&s4!=0)
		{
		    pnao++;
		    boolpnao=1;
		}
	    }
	}
    }
    print $0, boolpp, boolop, boolpe, booloe, boolpde, boolode, boolppoe, boolppode, boolpeop, boolpdeop, boolppop, boolpeoe, boolpdeode, boolpona, boolpnao, boolpnaona;
}

# - number and % of the relations with a pls on the prom side = pp ($19)
# - number and % of the relations with a pls on the other end side = op ($20)
# - number and % of the relations with an els on the prom side = pe ($21)
# - number and % of the relations with an els on the other end side = oe ($22)
# - number and % of the relations with a dels on the other end side = pde ($23)
# - number and % of the relations with a dels on the other end side = ode ($24)
# - number and % of the relations with a pls on the prom side and an els on the other end side = ppoe ($25)
# - number and % of the relations with a pls on the prom side and a dels on the other end side = ppode ($26)
# - number and % of the relations with an els on the prom side and a pls on the other end side = peop ($27)
# - number and % of the relations with a dels on the prom side and a pls on the other end side = pdeop ($28)
# - number and % of the relations with a pls on both sides = ppop ($29)
# - number and % of the relations with an els on both sides = peoe ($30)
# - number and % of the relations with a dels on both sides = pdeode ($31)
# - number and % of the relations with an overlap on the prom side but no overlap on the other end side = pona ($32)
# - number and % of the relations with an overlap on the other end side but no overlap on the prom side = pnao ($33)
# - number and % of the relations with neither a pls nor an els on any side = pnaona ($34)
END{
    print "pp", n, pp, pp/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "op", n, op, op/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "pe", n, pe, pe/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "oe", n, oe, oe/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "pde", n, pde, pde/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "ode", n, ode, ode/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "ppoe", n, ppoe, ppoe/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "ppode", n, ppode, ppode/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "peop", n, peop, peop/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "pdeop", n, pdeop, pdeop/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "ppop", n, ppop, ppop/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "peoe", n, peoe, peoe/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "pdeode", n, pdeode, pdeode/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "pona", n, pona, pona/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "pnao", n, pnao, pnao/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
    print "pnaona", n, pnaona, pnaona/n*100 > "pchic.16catwrtplsels.nb.pcent.tsv";
} 
