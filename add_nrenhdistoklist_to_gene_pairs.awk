# add_nrenhdistoklist_to_gene_pairs.awk

# Given:
########
# - a tsv file of gene tss coordinates with header (fileRef1)
# - a tsv file of contacted enhancers in many cell types for each gene with header (fileRef2)
# - a main tsv file with header with gene pairs
# adds to the main tsv file
###########################
# - the comma separated list of contacted enhancers of gene 1 that are at a distance between 25kb and 2M
#   from the tss of gene1 and at a distance between 25kb and 2Mb from the tss of gene2
# - the same for gene2

# example
# cd ~/work/duplicons/common_enhancers
# genetsscoord=~/work/duplicons/paper/supplementary_datasets_and_tables/supplementary_datasets10.tsv
# contactedenhancers=~/work/duplicons/paper/supplementary_datasets_and_tables/supplementary_datasets6.tsv
# pgm=~/fragencode/tools/multi/Scripts/add_nrenhdistoklist_to_gene_pairs.awk
# time awk -v fileRef1=$genetsscoord -v fileRef2=$contactedenhancers -f $pgm realparalogs.gn1.gn2.tss1.tss2.dist.tsv > realparalogs.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.tsv 
# real	0m16.363s

# fileRef1=$genetsscoord
# ID	Biotype	Chr	TSSPosition	MeanTPM
# ENSMUSG00000064372	Mt_tRNA	MT	15422	NA
# 56306 (5 fields)

# fileRef2=$contactedenhancers
# gene	EpiSC	ESC	ESC_18	ESC_NKO	ESC_wild	ESd_starved	ESd_TPO	FLC	preB_aged	preB_young	TSC	preadip_D0	preadip_D2	preadip_4H
# ENSMUSG00000000001	chr3:108101507:108101657;chr3:...

# main input file = realparalogs.gn1.gn2.tss1.tss2.dist.tsv
# gnid1	gnid2	tsspos1	tsspos2	tssdist
# ENSMUSG00000000056	ENSMUSG00000002280	11:121237253	17:25773776	NA
# 1421 (5 fields)

# main output file = realparalogs.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.tsv 
# gnid1	gnid2	tsspos1	tsspos2	tssdist	nrenhdistok1	nrenhdistok2
# ENSMUSG00000000056	ENSMUSG00000002280	11:121237253	17:25773776	NA	121266911:121267061,121265411:121265561,121442711:121442861,121439511:121439661,121442011:121442161,121269511:121269661,121263611:121263761,121439911:121440261,121446411:121446561,121343211:121343561,121342811:121342961,121380811:121381061,121343611:121343761,121341411:121341561,121376711:121376861,121377311:121377461,121368111:121368261,121366511:121366661,121346211:121346361,121345811:121345961,121367411:121367561,121351311:121351461,121347311:121347461,121349411:121349561,121348311:121348461,121346811:121346961,121808011:121808161,121808511:121808661,121430611:121430761,121362111:121362261,120982811:120982961,120983011:120983161,120981811:120981961,120982611:120982761,120982411:120982561,	25744780:25745030,25567280:25567430,25467580:25467830,25468180:25468330,25819280:25819430,25885380:25885530,26256480:26256630,26256180:26256330,26261480:26261630,
# 1421 (7 fields)

BEGIN{
    OFS="\t";
    # store the tss position of each gene
    while (getline < fileRef1 >0)
    {
	tsschr[$1]=($1!="MT" ? "chr"$3 : "chrM");
	tsspos[$1]=$4;
    }
    
    # store the nr list of enhancers at the ok distance from each gene, using the above tss info
    while (getline < fileRef2 >0)
    {
	c=tsschr[$1];
	p=tsspos[$1];
	n=split($0,a,"\t");
	for(i=2; i<=n; i++)
	{
	    if(a[i]!="")
	    {
		split(a[i],b,";");
		k=1;
		while(b[k]!="")
		{
		    # compute the position of the enhancer as its middle and then the distance from this middle to the gene tss
		    # and only keep the enhancers that are between 25kb and 2Mb from the gene in question
		    split(b[k],d,":");
		    pe=avg(d[2],d[3]);
		    dist=abs(pe-p);
		    if(d[1]==c&&(dist>=25000)&&(dist<=2000000))
		    {
			seen[$1,b[k]]++;
			if(seen[$1,b[k]]==1)
			{
			    nrenhlist[$1]=(nrenhlist[$1])(d[2]":"d[3])(",");
			}
		    }
		    k++;
		}
	    }
	}
    }
}

NR==1{
    print $0, "nrenhdistok1", "nrenhdistok2";
}

NR>=2{
    # here need to futher filter the two enhancer lists so that:
    # - the enhancers of gene1 are also between 25kb and 2Mb of the tss of gene2
    # - the enhancers of gene2 are also between 25kb and 2Mb of the tss of gene1

    # for the enh list of gene2 we ask the enh to also be at the right distance of the tss of gene 1
    s1="";
    split(nrenhlist[$1],a,",");
    split($4,b,":");
    p=b[2];
    k=1;
    while(a[k]!="")
    {
	split(a[k],a1,":");
	pe=avg(a1[1],a1[2]);
	dist=abs(pe-p);
	if(dist>=25000&&dist<=2000000)
	{
	    s1=(s1)(a[k])(",");
	}
	k++;
    }
    
    # for the enh list of gene2 we ask the enh to also be at the right distance of the tss of gene 1
    s2="";
    split(nrenhlist[$2],a,",");
    split($3,b,":");
    p=b[2];
    k=1;
    while(a[k]!="")
    {
	split(a[k],a1,":");
	pe=avg(a1[1],a1[2]);
	dist=abs(pe-p);
	if(dist>=25000&&dist<=2000000)
	{
	    s2=(s2)(a[k])(",");
	}
	k++;
    }
    print $0, nn(s1), nn(s2);
}


function avg(x1,x2)
{
    return (x1+x2)/2;
}

function abs(x)
{
    return (x>=0 ? x : (-1*x));
} 

function nn(x)
{
    return (x!="" ? x : "NA")
}
