# add_union_inter_enhlists.awk

# Given a tsv file with header that has gene pairs separated by a dist in column 5 (otherwise NA)
# and in the last two columns the nr comma separated list of enhancers
# contacted by each gene at a reasonable distance on the same chr
# for the gene pairs separated by less than 4Mb will add the list of nr enhancers in the union
# and in the intersection between those two lists

# example
# cd ~/work/duplicons/common_enhancers
# pgm=~/fragencode/tools/multi/Scripts/add_union_inter_enhlists.awk
# time awk -f $pgm realparalogs.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.tsv >  realparalogs.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.union.inter.tsv

# main input file = realparalogs.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.tsv
# gnid1	gnid2	tsspos1	tsspos2	tssdist	nrenhdistok1	nrenhdistok2
# ENSMUSG00000000056	ENSMUSG00000002280	11:121237253	17:25773776	NA	121266911:121267061,121265411:121265561,121442711:121442861,121439511:121439661,121442011:121442161,121269511:121269661,121263611:121263761,121439911:121440261,121446411:121446561,121343211:121343561,121342811:121342961,121380811:121381061,121343611:121343761,121341411:121341561,121376711:121376861,121377311:121377461,121368111:121368261,121366511:121366661,121346211:121346361,121345811:121345961,121367411:121367561,121351311:121351461,121347311:121347461,121349411:121349561,121348311:121348461,121346811:121346961,121808011:121808161,121808511:121808661,121430611:121430761,121362111:121362261,120982811:120982961,120983011:120983161,120981811:120981961,120982611:120982761,120982411:120982561,	25744780:25745030,25567280:25567430,25467580:25467830,25468180:25468330,25819280:25819430,25885380:25885530,26256480:26256630,26256180:26256330,26261480:26261630,
# 1421 (7 fields)

NR==1{
    OFS="\t";
    print $0, "enhunionlist", "enhinterlist";
}

NR>=2{
    # if the genes are on different chromosomes or more distant than 4Mb we dont even compute the union and inter lists of enhancers
    if($5=="NA"||$5>4000000)
    {
	print $0, "NA", "NA";
    }
    
    else
    {
	# the union and inter lists are initialised to empty
	ulist="";
	ilist="";

	# we go over the enh list of gene1 and for the current gene pair (row here) we store whether each enhancer was seen for gene1 in seen1
	# and we add to the union list only if we are seeing this enhancer for the 1st time for this gene pair
	split($6,a,",");
	k=1;
	while(a[k]!="")
	{
	    seen[NR,a[k]]++;
	    seen1[NR,a[k]]=1;
	    if(seen[NR,a[k]]==1)
	    {
		ulist=(ulist)(a[k])(",");
	    }
	    k++;
	}

	# we go over the enh list of gene2 and for the current gene pair (row here) we store whether each enhancer was seen for gene2 in seen2
	# we add to the union list only if we are seeing this enhancer for the 1st time for this gene pair
	# and then only if we have seen this enh in gene1 and we have not added it to the inter list, do we put it in this inter list 
	split($7,a,",");
	k=1;
	while(a[k]!="")
	{
	    seen[NR,a[k]]++;
	    if(seen[NR,a[k]]==1)
	    {
		ulist=(ulist)(a[k])(",");
	    }
	   if(seen1[NR,a[k]]==1&&seenboth[NR,a[k]]=="")
	    {
		seenboth[NR,a[k]]++;
		ilist=(ilist)(a[k])(",");
	    }
	   k++;
	}
	print $0, nn(ulist), nn(ilist);
    }
}

function nn(x)
{
    return (x!="" ? x : "NA");
}

