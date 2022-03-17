# compute.internalexlonger500bp.awk

# this script takes as input a gff file of internal exons of transcripts sorted by transcript id and then exon coordinates
# and with trid as second key in the 9th field and compute the number and % of transcripts of this file that have an internal
# exon longer than 500bp as well as the number and % of all and of distinct internal exons that are longer than 500bp.
# Note that the transcripts of this file all had the number of exons in this file plus two (the 5' and 3' not present here)

# example
# cd ~/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/TAGADA.v1.0.1.GRCg6a.102.21-06-26/annotation
# pgm=~/fragencode/tools/multi/Scripts/compute.internalexlonger500bp.awk
# time zcat novel.genes_supported_intron_chains.tr3exons.internalexons.gff.gz | awk -f $pgm  
# nbtr3ex	withinternalexlonger500bp.nb	withinternalexlonger500bp.pcent	nbinternalex	longer500bp.nb1	longer500bp.pcent1	nbdistinternalex	longer500bp.nb2	longer500bp.pcent2
# 75060	21153	28.1815	831345	27650	3.32594	187810	10421	5.54869
# real	0m2.274s


{
    # compute the total number of transcripts in this file (ie of transcripts with at least 3 exons)
    seentr[$12]++;
    if(seentr[$12]==1)
    {
	nbtr++;
    }

    # compute the total number of exons in this file (ie of internal exons of the tr of this file)
    nbex++;

    # compute the number of distinct exons in this file (ie of distinct internal exons of the tr of this file)
    seenex[$1":"$4":"$5]++;
    if(seenex[$1":"$4":"$5]==1)
    {
	nbdistex++;
    }

    if(($5-$4+1)>=500)
    {
	nbex1++;
	seenex1[$1":"$4":"$5]++;
	if(seenex1[$1":"$4":"$5]==1)
	{
	    nbdistex1++;
	}
	seentr1[$12]++;
	if(seentr1[$12]==1)
	{
	    nbtr1++;
	}
    }
}

END{
    OFS="\t";
    print "nbtr3ex", "withinternalexlonger500bp.nb", "withinternalexlonger500bp.pcent", "nbinternalex", "longer500bp.nb1", "longer500bp.pcent1", "nbdistinternalex", "longer500bp.nb2", "longer500bp.pcent2";
    print nbtr, nbtr1, nbtr1/nbtr*100, nbex, nbex1, nbex1/nbex*100, nbdistex, nbdistex1, nbdistex1/nbdistex*100;
}
