# add_utrs_to_gtffile.awk
# this is to add 5' and 3' utrs to a gtf file output by transdecoder using an intermediate gff3 file also produced upstreamd by transdecoder
# the particularity of the gtf file is that in the 9th field we first have transcript_id and then gene_id and also the trid has .pX after its real id
# and the gnid has ^... after its real id, so need to remove those, also I can remove all thekey values after trid and gnid since they have spaces

# !!! be careful: even if we want to eliminate .pX in trid, need to consider it to make the correspondance with the utrs !!!
# !!! be careful: some trid can have . for example G4.11 !!!

# example
# srun --x11 --mem=8G --pty bash
# cd /work/project/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/nf-feelnc.tagiso.23-03-03/results/annotation/transdecoder
# pgm1=~/fragencode/tools/multi/Scripts/add_utrs_to_gtffile.awk
# pgm2=~/fragencode/tools/multi/Scripts/gff2gff.awk
# pgm3=~/fragencode/tools/multi/Scripts/make_gff_ok.awk
# awk -v fileRef=tagisotr.fasta.transdecoder.genome.gff3 -f $pgm1 tagisotr.fasta.transdecoder.genome.bestcds.gtf | awk -f $pgm2 | awk -f $pgm3 > tagisotr.fasta.transdecoder.genome.bestcds.withutrs.gtf

# fileRef=tagisotr.fasta.transdecoder.genome.gff3
# 13	transdecoder	five_prime_UTR	65721467	65721530	.	+	.	ID=G9282.1.p1.utr5p1;Parent=G9282.1.p1

# main input = tagisotr.fasta.transdecoder.genome.bestcds.gtf
# 1	transdecoder	CDS	27072	27346	.	-	2	transcript_id "G3.1.p2"; gene_id "LOC_000000090338^1^-"; gene_name "ORF type:5prime_partial (+),score=17.47";
# 1	transdecoder	CDS	27608	27767	.	-	0	transcript_id "G3.1.p2"; gene_id "LOC_000000090338^1^-"; gene_name "ORF type:5prime_partial (+),score=17.47";
# 3187630 (16 fields)

# output
# 1	transdecoder	CDS	27072	27346	.	-	2	gene_id "LOC_000000090338"; transcript_id "G3.1"; 
# 1	transdecoder	three_prime_UTR	26468	27071	.	-	.	gene_id "LOC_000000090338"; transcript_id "G3.1"; 
# 3753767 (12 fields) 



BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	if($3~/UTR/)
	{
	    split($9,a,";");
	    split(a[1],b,"=");
	    split(b[2],c,".p");
	    split(c[2],d,".");
	    trid=c[1]".p"d[1];
	    
	    if($3~/three/)
	    {
		nbtputr[trid]++;
		tputr[trid,nbtputr[trid]]=$1":"$4":"$5":"$7;
	    }
	    else
	    {
		if($3~/five/)
		{
		    nbfputr[trid]++;
		    fputr[trid,nbfputr[trid]]=$1":"$4":"$5":"$7;
		}
	    }
	}
    }
}


{
    # print the exon and cds rows first
    split($10,a,"\"");
    split(a[2],b,".p");
    $10="\""b[1]"\"\;";
    trid=b[1]".p"b[2];
     
    split($12,a,"\"");
    split(a[2],b,"^");
    $12="\""b[1]"\"\;";
    
    for(i=13; i<=NF; i++)
    {
	$i="";
    }
    print;

    # and now prints the 5' and 3' utr rows
    seen[$10]++;
    if(seen[$10]==1)
    {
	if(nbtputr[trid]!="")
	{
	    for(i=1; i<=nbtputr[trid]; i++)
	    {
		split(tputr[trid,i],b,":");
		print b[1], "transdecoder", "three_prime_UTR", b[2], b[3], ".", b[4], ".", "transcript_id", $10, "gene_id", $12;
	    }
	}
	if(nbfputr[trid]!="")
	{
	    for(i=1; i<=nbfputr[trid]; i++)
	    {
		split(fputr[trid,i],b,":");
		print b[1], "transdecoder", "five_prime_UTR", b[2], b[3], ".", b[4], ".", "transcript_id", $10, "gene_id", $12;
	    }
	}
    }
}

