# conn2connwithsnplists.awk
# takes as input
# - an ext variable that indicates that for the fileRef file there has been an extension of ext on each side of each element wrt the coord in the 11 column tsv file
# - a fileRef file that is a tsv file with two columns with no header that has a prom or elt2 coord as 1st column (chr:beg:end, bed coord)
#   and the list of snps (chr:pos) overlapped by it as 2nd column
# - an 11+ column tsv file without header that has
#   * prom coord list
#   * number of prom
#   * elt2 coord list
#   * number of elt2
#   * prom-elt2 score list
#   * prom cumul length
#   * elt2 cumul length
#   * prom and elt2 cumul length
#   * min dist between a prom and an elt2
#   * max dist between a prom and an elt2
#   * elt2 coord list for elt2 equal to prom in the object
# and outputs the same 11+ column tsv file as the one given as input but with 3 additional columns representing:
#  * comma separated list of unique snps (chr:pos) in proms according to fileRef
#  * comma separated list of unique snps (chr:pos) in elt2s according to fileRef
#  * comma separated list of unique snps in proms and elt2s according to fileRef
# Practical note:
#################
# the fileRef has the extension done on the coord while the 11+ field tsv input file does not
# and while it would have been easier and less time consuming for the script to put the coord
# of fileRef right (adding ext on beg and substracting ext on end), it is not possible since we lost info of
# the real beg when substracting more than the beg and has to put 0 as beg in bed file and therefore fileRef,
# therefore we will extend the coord of the 11+ field input tsv file although more complicated and more time consuming since in loops

# example:
##########
# outdir=~/egingwas/parkinson/intersection.with.jung.ren.2019
# vcf=~sdjebali/data/parkinson/2017/HRC_Imputations/SANGER/allchr.basic.INFOin0.5-1.ACin60-5086.notstrambig.uniqchrompos.ok.minor.unfreq.info.vcf.gz
# pgm=~/tools/multi/Scripts/conn2connwithsnplists.awk
# ct=HCmerge
# obj=allpromallelt2
# tsv=$outdir/$ct/$ct.pall.score2.ens.$obj.prom.elt2.list.nb.sclist.prom.elt2.cumullength.sum.dist.min.max.elt2equprom.allpostqcsnp.nb1.nb2.nb3.0kb.10kb.100kb.tsv
# time awk -v ext=0 -v fileRef=$outdir/$ct/$ct.pall.score2.ens.promelt2.prom.elt2.uniq.segcoord.allpostqcsnplist.tsv -f $pgm $tsv > $outdir/$ct/$ct.pall.score2.ens.$obj.prom.elt2.list.nb.sclist.prom.elt2.cumullength.sum.dist.min.max.elt2equprom.allpostqcsnp.nb1.nb2.nb3.0kb.10kb.100kb.list1.list2.list3.0kb.tsv

# inputs:
#########
# $outdir/$ct/$ct.pall.score2.ens.promelt2.prom.elt2.uniq.segcoord.allpostqcsnplist.tsv
#######################################################################################
# 16:64054543:64059837	16:64055436,16:64055592,16:64055734,16:64056384,16:64056536,16:64056640,16:64056924,16:64057076,16:64057244,16:64058296,16:64058883,16:64059022,16:64059054,16:64059158,
# 8:40731183:40737464	8:40732000,8:40732477,8:40733020,8:40733289,8:40733345,8:40733705,8:40734130,8:40734250,8:40734576,8:40734662,8:40734706,8:40734746,8:40735728,8:40736242,8:40736347,8:40736382,8:40736395,8:40736663,
# 78323 (2 fields)
# $tsv is
#########
# 1:1183838:1209657,1:1331882:1357608,1:1379293:1385838,1:1435680:1451318,1:1502828:1517794,1:1558274:1582294,1:1584751:1621263,1:1647936:1674371,1:1674372:1684389,1:1704140:1729970,1:1819361:1831600,	11	1:1491937:1496017,1:1647936:1674371,1:1435680:1451318,1:1757202:1764600,1:1814692:1819360,1:1331882:1357608,1:1496018:1502827,1:1752561:1757201,1:1183838:1209657,1:1379293:1385838,1:1502828:1517794,1:1539937:1549357,1:1558274:1582294,1:1584751:1621263,1:1621264:1634847,1:1896933:1910676,1:1739815:1747632,	17	2.18497234006452,2.18346841284655,2.50607490345802,2.45370189969025,2.36165800236146,3.32816882381377,3.14906502956605,2.5106741974425,2.29507495792811,2.35460879825405,6.00647209644406,3.33034565728985,4.68596905408866,2.08116070475416,2.41291622150509,2.1717770163861,2.18026919137198,3.32463272367633,2.35083687076021,3.32164185204039,2.46013680123846,2.07802431392853,2.41005829485489,92.3459964117904,4.00235507786993,2.43148040906616,12.4760314900563,2.11419909995048,2.49565164085615,3.86588299256036,	223747	247819	295905	1	687276	1:1183838:1209657,1:1331882:1357608,1:1379293:1385838,1:1435680:1451318,1:1502828:1517794,1:1558274:1582294,1:1584751:1621263,1:1647936:1674371,	267	312	360	544	641	733	1588	1764	1764
# 4640 (20 fields)


# output
########
# $outdir/$ct/$ct.pall.score2.ens.$obj.prom.elt2.list.nb.sclist.prom.elt2.cumullength.sum.dist.min.max.elt2equprom.allpostqcsnp.nb1.nb2.nb3.0kb.10kb.100kb.tsv
# 1:1621264:1634847,	1	1:1674372:1684389,	1	92.2880184499934,	13583	10017	23600	39525	39525	NA	2	5	7	15	41	56	140	198	199	1:1624396,1:1627805,	1:1674767,1:1675151,1:1680219,1:1683565,1:1683900,	1:1624396,1:1627805,1:1674767,1:1675151,1:1680219,1:1683565,1:1683900,
# 4640 (23 fields)


BEGIN{
    # read the file of segments with list of snps inside them and to each segment coordinates associate its list of overlapping snps
    OFS="\t";
    if(ext=="")
    {
	ext=0;
    }
    while (getline < fileRef >0)
    {
	
	snplist[$1]=$2;
    }
}

{
    # initialize the comma separated list of unique snps in proms, elt2s and proms and elt2s
    s1="";
    s2="";
    s3="";
    # split the prom list into indiv prom (coord)
    split($1,a,",");
    # go over the indiv prom of the prom list
    k=1;
    while(a[k]!="")
    {
	# extending the prom coord of ext on each side and putting beg to 0 if beg becomes negative doing so
	split(a[k],a0,":");
	newcoord=a0[1]":"((a0[2]-ext)<0 ? 0 : (a0[2]-ext))":"(a0[3]+ext);
	# go over the list of snps overlapped by the current indiv prom extended on each side
	split(snplist[newcoord],a1,",");
	l=1;
	while(a1[l]!="")
	{
	    # increment the lists of unique snps of proms and proms and elt2s only if it is the 1st time the snp is seen for these categories
	    seen1[NR,a1[l]]++;
	    seen3[NR,a1[l]]++;
	    if(seen1[NR,a1[l]]==1)
	    {
		s1=(s1)(a1[l])(",");
	    }
	    if(seen3[NR,a1[l]]==1)
	    {
		s3=(s3)(a1[l])(",");
	    }
	    l++;
	}
	k++;
    }
    # in the same way as for prom, split the elt2 list into indiv elt2 (coord)
    split($3,b,",");
    # go over the indiv elt2 of the elt2 list
    k=1;
    while(b[k]!="")
    {
	# extending the elt2 coord of ext on each side and putting beg to 0 if beg becomes negative doing so
	split(b[k],b0,":");
	newcoord=b0[1]":"((b0[2]-ext)<0 ? 0 : (b0[2]-ext))":"(b0[3]+ext);
	# go over the list of snps overlapped by the current indiv elt2 extended on each side
	split(snplist[newcoord],b1,",");
	l=1;
	while(b1[l]!="")
	{
	    # increment the lists of unique snps of elt2 and proms and elt2s only if it is the 1st time the snp is seen for these categories
	    seen2[NR,b1[l]]++;
	    seen3[NR,b1[l]]++;
	    if(seen2[NR,b1[l]]==1)
	    {
		s2=(s2)(b1[l])(",");
	    }
	    if(seen3[NR,b1[l]]==1)
	    {
		s3=(s3)(b1[l])(",");
	    }
	    l++;
	}
	k++;
    }
    print $0, na(s1), na(s2), na(s3);
}

function na(x){
    return (x!="" ? x : "NA")
}
   
