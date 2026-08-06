# nb.pcent.hichip.relations.with.pls.els.awk
# given a hichip bedpe file where the two fragments are given in gx order
# and with 6 last fields (columns 9 to 14) including booleans saying whether or not the
# 1st part (in gx order) and the 2nd part (in gx order) overlap ccre pls, els and dels
# (distal els) respectively, will provide:
# 1. in a file called hichip.10catwrtplsels.nb.pcent.tsv, the number and pcent of all relations
#    falling in each of the 10 following categories:
#    - number and % of the relations with a pls on one side = p ($15)
#    - number and % of the relations with an els on one side = e ($16)
#    - number and % of the relations with a dels on one side = de ($17)
#    - number and % of the relations with a pls on one side and an els on the other side = pe ($18)
#    - number and % of the relations with a pls on one side and a dels on the other side = pde ($19)
#    - number and % of the relations with a pls on both sides = pp ($20)
#    - number and % of the relations with an els on both sides = ee ($21)
#    - number and % of the relations with a dels on both sides = dede ($22)
#    - number and % of the relations with a pls or els on one side but nothing on the other = na ($23)
#    - number and % of the relations with neither a pls nor an els on any side = nana ($24)
# 2. as output the same file as the input but with 10 additional columns that are booleans about these 10 categories
# !!! Note that we have 1/5 of the frag that overlap an els that overlap both a pels and a dels !!!
# !!! this is high but the hichip frag are quite big here (5kb) !!!
# !!! need to see if we need to specifically ask that when over dels it does not overlap pels !!!

# example
# cd ~/bridge/results/hichip/homo_sapiens/hg19/wilson_2020
# pgm=~/fragencode/tools/multi/Scripts/nb.pcent.hichip.relations.with.pls.els.awk
# time awk -f $pgm hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.bedpe > hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.10catwrtplsels.bedpe
# real	0m0.150s

# input = hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.bedpe
# chr1	1250000	1255000	chr1	1310000	1315000	chr1:1250000-1255000?chr1:1310000-1315000	7.139784871	0	1	1	1	1	0
# chr1	1290000	1295000	chr1	1365000	1370000	chr1:1290000-1295000?chr1:1365000-1370000	11.04504867	1	1	1	1	1	1
# 61136 (14 fields)

# main output = hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.10catwrtplsels.bedpe
# chr1	1250000	1255000	chr1	1310000	1315000	chr1:1250000-1255000?chr1:1310000-1315000	7.139784871	0	1	1	1	1	0	1	1	1	1	1	0	10	0	0
# chr1	1290000	1295000	chr1	1365000	1370000	chr1:1290000-1295000?chr1:1365000-1370000	11.04504867	1	1	1	1	1	1	1	1	1	1	1	1	11	0	0
# 61136 (24 fields)

# secondary output = hichip.10catwrtplsels.nb.pcent.tsv
# p     61136  36580  59.8338
# e     61136  60281  98.6015
# de    61136  56263  92.0292
# pe    61136  32977  53.9404
# pde   61136  30100  49.2345
# pp    61136  8174   13.3702
# ee    61136  52971  86.6445
# dede  61136  35794  58.5482
# na    61136  7147   11.6903
# nana  61136  831    1.35926


BEGIN{
    OFS="\t";
}

{
    boolp=0;
    boole=0;
    boolde=0;    # 3 classes for the 1 frag categories (no cat with na here)
    boolpe=0;
    boolpde=0;
    boolpp=0;
    boolee=0;
    booldede=0;
    boolna=0;
    boolnana=0;   # 7 calsses for the 2 frag categories (2 last cat with nas here)
    n++;

    # $9=1stp / $10=1ste / $11=1stde / $12=2ndp / $13=2nde / $14=2ndde 
    # check first whether no overlap with pls and els on each side, in which case naona is incremented
    s=$9+$10+$12+$13;
    if(s==0)
    {
	boolnana=1;
	nana++;
    }
    else
    {
	if($9==1||$12==1)
	{
	    boolp=1;
	    p++;
	    if($9==1&&$12==1)
	    {
		boolpp=1;
		pp++;
	    }
	}
	if($10==1||$13==1)
	{
	    boole=1;
	    e++;
	    if($10==1&&$13==1)
	    {
		boolee=1;
		ee++;
	    }
	}
	if($11==1||$14==1)
	{
	    boolde=1;
	    de++;
	    if($11==1&&$14==1)
	    {
		booldede=1;
		dede++;
	    }
	}
	if(($9==1&&$13==1)||($10==1&&$12==1))
	{
	    boolpe=1;
	    pe++;
	}
	if(($9==1&&$14==1)||($11==1&&$12==1))
	{
	    boolpde=1;
	    pde++;
	}
	# at this point we could still have relations where the sum of all 2 frag categories is 0 while not being nana
	# those are cases where either the 1st part or the 2nd part overlaps a ccre but not the other, so one category
	# is added for them, called na, that we compute here
	s2=boolpe+oolpde+boolpp+boolee+booldede;
	s3=boolp+boole;
	if(s2==0)
	{
	    if(s3!=0)
	    {
		na++;
		boolna=1;
	    }
	}
    }
    print $0, boolp, boole, boolde, boolpe, boolpde, boolpp, boolee, booldede, boolna, boolnana;
}

# - number and % of the relations with a pls on one side = p ($15)
# - number and % of the relations with an els on one side = e ($16)
# - number and % of the relations with a dels on one side = de ($17)
# - number and % of the relations with a pls on one side and an els on the other side = pe ($18)
# - number and % of the relations with a pls on one side and a dels on the other side = pde ($19)
# - number and % of the relations with a pls on both sides = pp ($20)
# - number and % of the relations with an els on both sides = ee ($21)
# - number and % of the relations with a dels on both sides = dede ($22)
# - number and % of the relations with a pls or els on one side but nothing on the other = na ($23)
# - number and % of the relations with neither a pls nor an els on any side = nana ($24)
END{
    print "p", n, p, p/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
    print "e", n, e, e/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
    print "de", n, de, de/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
    print "pe", n, pe, pe/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
    print "pde", n, pde, pde/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
    print "pp", n, pp, pp/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
    print "ee", n, ee, ee/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
    print "dede", n, dede, dede/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
    print "na", n, na, na/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
    print "nana", n, nana, nana/n*100 > "hichip.10catwrtplsels.nb.pcent.tsv";
} 
