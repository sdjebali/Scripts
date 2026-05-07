# add_4booleans_to_pchicbedpe_relative_to_hichip_overlap.awk
# In the context of the comparison between pchic and hichip relations found on the same sample (with the same restriction enzyme)
# I need to add to the pchic bedpe file (only cis relations) ordered by gx order, and using two additional files, 4 booleans:
# - whether each pchic fragment had to be extended to reach 5kb (in case it was less than that, 2 booleans)
# - whether the only or the two promoter fragments extended to 5kb were overlapping a hichip fragment (1 bool)
# - whether it overlaps an hichip relation (1 bool)
# the additional files used to add the 2 last booleans are
# - a bed file with the pchic promoter fragments extended to 5kb (when needed) that overlap a hichip fragment
# - a kind of bedpe file of pchic whole segments that actually overlapped a hichip relation when extended to 5kb when needed

# example of application on genotoul
# srun --x11 --mem=8G --pty bash
# cd /work/project/bridge/workspace/sdjebali/enhancer.gene/3D.comparison/ipsc_cardiomyocytes
# pgm=~/fragencode/tools/multi/Scripts/add_4booleans_to_pchicbedpe_relative_to_hichip_overlap.awk
# pchic=/work/project/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018/pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.bedpe
# time awk -v fileRef1=pchicfragext5kb.over.hichipfrag.bed -v fileRef2=pchic.ipsccm.montefiori2018.hg19.overgencv49tss.fragextto5kb.onesegment.over.hichip.over.othertech.tsv -f $pgm $pchic > pchic.all.relations.onefragoverpchicpromext5kb.overpchic.bedpe
# real	0m1.822s  

# fileRef1=pchicfragext5kb.over.hichipfrag.bed
# chr1	1291367	1296366
# 8203 (3 fields)   *** here we have prom pchic fragments extended by 5kb (when needed) that overlap a hichip fragment

# fileRef2=pchic.ipsccm.montefiori2018.hg19.overgencv49tss.fragextto5kb.onesegment.over.hichip.over.othertech.tsv
# chr15	74419473	74600381	chr15:74419473_74424472:74595381_74600381:175908:ISLR2*NM:p-o:init.order:1:0	20.07	.	chr15,74420000,74605000,chr15:74420000_74425000:74600000_74605000:180000,8.303520593,.	180381
# 26871 (8 fields)   *** here we have the pchic whole segments that are overlapped by hichip relations with extension of 5kb

# main input file = $pchic 
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order
# 401098 (12 fields)  *** here we have the initial pchic relations (with no extension to 5kb of the fragments) with all information

# main output file = pchic.all.relations.onefragoverpchicpromext5kb.overpchic.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	1	1	0	0
# chr1	834807	835706	chr1	894208	897150	6.58	KLHL17*NM	1	0	p-o	rev.order	1	1	0	0
# 401098 (16 fields) 

BEGIN{
    OFS="\t";
    # from the 1st file we get the coordinates of the pchic fragments extended to 5kb that overlapped hichip fragments
    while (getline < fileRef1 >0)
    {
	okpromfrag[$1":"$2"_"$3]=1;
    }
    # from the 2nd file we get the coordinates of the pchic relations in gx order that are actually overlapped by hichip when extended to 5kb
    while (getline < fileRef2 >0)
    {
	split($4,a,":");
	over[a[1]":"a[2]":"a[3]]=1;
    }
}

# here we read the actual bedpe file of pchic in gx order and we add the 4 booleans to it
# in particular those 2 more difficult at the end
# - whether the only or the two promoter fragments extended to 5kb were overlapping a hichip fragment (o1, 1 bool)
# - whether it overlaps an hichip relation (o2, 1 bool)
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order
{
    # we first extend the two fragments since this is needed many times
    e1=ext($2,$3);
    e2=ext($5,$6);
    
    # we first try and compute o1, based on whether we have a p-o or a p-p relation and in the 1st case whether we are in the init or rev order
    # in order to know which fragment is the prom one
    o1=0;
    if($11=="p-o")
    {
	if($12=="init.order")
	{
	    o1=(o1||okpromfrag[$1":"e1]);
	}
	else
	{
	    if($12=="rev.order")
	    {
		o1=(o1||okpromfrag[$4":"e2]);
	    }
	}
    }
    else
    {
	if($11=="p-p")
	{
	    o1=(o1||okpromfrag[$1":"e1]||okpromfrag[$4":"e2]);
	}
    }

    # o2 is given by the over hashtable, correct formatting of the key (order should be the same as in bedpe file) and extension of the fragments when needed
    # we expect something like this chr15:74419473_74424472:74595381_74600381 (gx order with 5kb extension)
    # and start from the bedpe file with no extension in the same order (gx order, w/o 5kb extension)
    o2=(over[$1":"e1":"e2]==1 ? 1 : 0);
    
    # and then we print the pchic relation with all info
    print $0, isext($3-$2), isext($6-$5), o1, o2;
}


# this is a function that says whether a segment has to be extended to 5kb
function isext(x){
    return (x>=5000 ? 0 : 1);
}

# this is a function that returns a 5kb extended fragment when it needs it, exactly in the same way as was done for the intersection
function ext(b,e){
    if((e-b)>=5000)
    {
	return b"_"e;
    }
    else
    {
	d=(5000-(e-b));
	toadd=int(d/2);
	return (b-toadd)"_"(e+toadd);
    }
}
