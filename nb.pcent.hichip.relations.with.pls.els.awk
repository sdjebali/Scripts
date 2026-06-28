# nb.pcent.hichip.relations.with.pls.els.awk
# given a hichip bedpe file where the two fragments are given in gx order
# and with 6 last fields (columns 9 to 14) including booleans saying whether or not the
# 1st part (in gx order) and the 2nd part (in gx order) overlap ccre pls, els and dels
# (distal els) respectively, will provide:
# 1. in a file called hichip.9catwrtplsels.nb.pcent.tsv, the number and pcent of all relations
#    falling in each of the 9 following categories:
#    - number and % of the relations with a pls on one side = pls ($15)
#    - number and % of the relations with an els on one side = els ($16)
#    - number and % of the relations with a dels on one side = dels ($17)
#    - number and % of the relations with a pls on one side and an els on the other side = plsels ($18)
#    - number and % of the relations with a pls on one side and a dels on the other side = plsdels ($19)
#    - number and % of the relations with a pls on both sides = plspls ($20)
#    - number and % of the relations with an els on both sides = elsels ($21)
#    - number and % of the relations with a dels on both sides = delsdels ($22)
#    - number and % of the relations with neither a pls nor an els on any side = nana ($23)
# 2. as output the same file as the input but with 9 additional columns that are booleans about these 9 categories
# !!! Note that we have 1/5 of the frag that overlap an els that overlap both a pels and a dels !!!
# !!! this is high but the hichip frag are quite big here (5kb) !!!
# !!! need to see if we need to specifically ask that when over dels it does not overlap pels !!!

# example
# cd ~/bridge/results/hichip/homo_sapiens/hg19/wilson_2020
# pgm=~/fragencode/tools/multi/Scripts/nb.pcent.hichip.relations.with.pls.els.awk
# time awk -f $pgm hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.bedpe > hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.9catwrtplsels.bedpe
# real	0m0.150s

# input = hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.bedpe
# chr1	1250000	1255000	chr1	1310000	1315000	chr1:1250000-1255000?chr1:1310000-1315000	7.139784871	0	1	1	1	1	0
# chr1	1290000	1295000	chr1	1365000	1370000	chr1:1290000-1295000?chr1:1365000-1370000	11.04504867	1	1	1	1	1	1
# 61136 (14 fields)

# main output = hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.9catwrtplsels.bedpe
# chr1	1250000	1255000	chr1	1310000	1315000	chr1:1250000-1255000?chr1:1310000-1315000	7.139784871	0	1	1	1	1	01	1	1	1	1	0	1	0	0
# chr1	1290000	1295000	chr1	1365000	1370000	chr1:1290000-1295000?chr1:1365000-1370000	11.04504867	1	1	1	1	1	11	1	1	1	1	1	1	1	0
# 61136 (23 fields)

# secondary output = hichip.9catwrtplsels.nb.pcent.tsv
# pls       61136  36580  59.8338
# els       61136  60281  98.6015
# dels      61136  56263  92.0292
# plsels    61136  32977  53.9404
# plsdels   61136  30100  49.2345
# plspls    61136  8174   13.3702
# elsels    61136  52971  86.6445
# delsdels  61136  35794  58.5482
# nana      61136  831    1.35926


BEGIN{
    OFS="\t";
}

{
    boolpls=0;
    boolels=0;
    booldels=0;
    boolplsels=0;
    boolplsdels=0;
    boolplspls=0;
    boolelsels=0;
    booldelsdels=0;
    boolnana=0;
    n++;

    # $9=1stpls / $10=1stels / $11=1stdels / $12=2ndpls / $13=2ndels / $14=2nddels 
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
	    boolpls=1;
	    pls++;
	    if($9==1&&$12==1)
	    {
		boolplspls=1;
		plspls++;
	    }
	}
	if($10==1||$13==1)
	{
	    boolels=1;
	    els++;
	    if($10==1&&$13==1)
	    {
		boolelsels=1;
		elsels++;
	    }
	}
	if($11==1||$14==1)
	{
	    booldels=1;
	    dels++;
	    if($11==1&&$14==1)
	    {
		booldelsdels=1;
		delsdels++;
	    }
	}
	if(($9==1&&$13==1)||($10==1&&$12==1))
	{
	    boolplsels=1;
	    plsels++;
	}
	if(($9==1&&$14==1)||($11==1&&$12==1))
	{
	    boolplsdels=1;
	    plsdels++;
	}
    }
    print $0, boolpls, boolels, booldels, boolplsels, boolplsdels, boolplspls, boolelsels, booldelsdels, boolnana;
}

# - number and % of the relations with a pls on one side = pls ($15)
# - number and % of the relations with an els on one side = els ($16)
# - number and % of the relations with a dels on one side = dels ($17)
# - number and % of the relations with a pls on one side and an els on the other side = plsels ($18)
# - number and % of the relations with a pls on one side and a dels on the other side = plsdels ($19)
# - number and % of the relations with a pls on both sides = plspls ($20)
# - number and % of the relations with an els on both sides = elsels ($21)
# - number and % of the relations with a dels on both sides = delsdels ($22)
# - number and % of the relations with neither a pls nor an els on any side = nana ($23)
END{
    print "pls", n, pls, pls/n*100 > "hichip.9catwrtplsels.nb.pcent.tsv";
    print "els", n, els, els/n*100 > "hichip.9catwrtplsels.nb.pcent.tsv";
    print "dels", n, dels, dels/n*100 > "hichip.9catwrtplsels.nb.pcent.tsv";
    print "plsels", n, plsels, plsels/n*100 > "hichip.9catwrtplsels.nb.pcent.tsv";
    print "plsdels", n, plsdels, plsdels/n*100 > "hichip.9catwrtplsels.nb.pcent.tsv";
    print "plspls", n, plspls, plspls/n*100 > "hichip.9catwrtplsels.nb.pcent.tsv";
    print "elsels", n, elsels, elsels/n*100 > "hichip.9catwrtplsels.nb.pcent.tsv";
    print "delsdels", n, delsdels, delsdels/n*100 > "hichip.9catwrtplsels.nb.pcent.tsv";
    print "nana", n, nana, nana/n*100 > "hichip.9catwrtplsels.nb.pcent.tsv";
} 
