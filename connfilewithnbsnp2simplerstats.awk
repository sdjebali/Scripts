# connfilewithnbsnp2simplerstats.awk
# from a tsv file of connections with nb of snps in first part, in second part and in total (with no redund)
# in columns fldno1, fldno2 and fldnotot respectively, a cell type ct, a type of object obj, it outputs a tsv row with
# - cell type
# - object type
# - total number of objects
# - number and % of objects with 0 snp
# - number and % of objects with 1 snp
# - number and % of objects with 2 snps but on one side
# - number and % of objects with 3+ snps but on one side
# - number and % of objects with 2 snps with one on each side
# - number and % of objects with 3+ snps with one on each side


# outdir=~/egingwas/parkinson/intersection.with.jung.ren.2019
# fldno=14
# ext=1.0kb
# filt=all
# ct=HCmerge
# obj=promelt2
# pgm=~/tools/multi/Scripts/connfilewithnbsnp2simplerstats.awk
# cd $outdir/$ct
# time awk -v ct=$ct -v obj=$obj -v fldno1=$((fldno-2)) -v fldno2=$((fldno-1)) -v fldnotot=$fldno -f $pgm $ct.pall.score2.ens.$obj.prom.elt2.list.nb.sclist.prom.elt2.cumullength.sum.dist.min.max.elt2equprom.$filt\postqcsnp.nb1.nb2.nb3.0kb.10kb.100kb.tsv

# input is like this
# 1:1183838:1209657	1	1:1491937:1496017	1	2.18497234006452,	25819	4080	29899	282280	282280	NA	56	4	60	121	53	174	665	212	877
# 1:1183838:1209657	1	1:1647936:1674371	1	2.18346841284655,	25819	26435	52254	438279	438279	NA	56	36	92	121	47	168	665	182	847
# 64514 (20 fields)  *** 11 first fields = promlist, nbprom, elt2list, nbelt2, scorelist, promlistcumullength, elt2listcumullength, sumofthetwo, mindist, maxdist, listofelt2equaltoprom
#                    *** 9 next fields = nbsnpinprom, nbsnpinelt2, nbsnpinboth at 0kb extension and then the same for 10kb extension and for 100kb extension


{
    tot++;
    if($fldnotot==0)
    {
	nb0++;
    }
    else
    {
	if($fldnotot==1)
	{
	    nb1++;
	}
	else
	{
	    if($fldnotot==2)
	    {
		if($fldno1>=1&&$fldno2>=1)
		{
		    nb22++;
		}
		else
		{
		    nb21++
		}
	    }
	    else
	    {
		if($fldnotot>=3)
		{
		   if($fldno1>=1&&$fldno2>=1)
		   {
		       nb32++;
		   }
		   else
		   {
		       nb31++;
		   }
		}
	    }
	}
    }
}

END{
    OFS="\t";
    print ct"\t"obj"\t"tot"\t"(nn(nb0))("\t")(nb0/tot*100)("\t")(nn(nb1))("\t")(nb1/tot*100)("\t")(nn(nb21))("\t")(nb21/tot*100)("\t")(nn(nb31))("\t")(nb31/tot*100)("\t")(nn(nb22))("\t")(nb22/tot*100)("\t")(nn(nb32))("\t")(nb32/tot*100);
}

function nn(x){
    return (x=="" ? 0 : x);
}

