# read.through.intron.to.n-1to2ndex.awk
# takes as input a gff file of potential read-through introns with indication of 1st and 2nd exon it is connecting in column 18 and 20
# as well as a fileRef tsv file of 5 columns for unique exons with indication of list of transcript ids, list of number of exons of these
# transcripts, list of gx positions of the exon in these transcripts and list of biological (5' to 3') positions of the exons in these
# transcripts, and outputs the same intron file but with indication of whether the intron is connecting the n-1th exon of a transcript
# (in biological order, i.e 5' to 3') and the 2nd exon of another transcript
# !!! be careful to deal with the monoex tr with no strand as well !!!

# example
# cd ~/fragencode/workspace/sdjebali/geneswitch/analysis/read.throughs/gallus_gallus
# pgm=~/fragencode/tools/multi/Scripts/read.through.intron.to.n-1to2ndex.awk
# time awk -v fileRef=excoord.tr.nbex.gxpos.biopos.lists.tsv -f $pgm stringtie_introns_onerefgneachside_diff_can_ds_adjrefgn12.gff > stringtie_introns_onerefgneachside_diff_can_ds_adjrefgn12_n-1to2ndex.gff

# fileRef file
# 17:5446192:5446265:+	MSTRG.5999.7,	16,	1,	1,
# 19:6113577:6121620:-	MSTRG.6746.6,	4,	1,	4,
# 346836 (5 fields)

# intron file
# 1	StringTie	intron	39868	258714	.	-	.	gene_id "MSTRG.1"; transcript_id "MSTRG.1.10"; refgn1 "ENSGALG00000042023"; refgn2 "ENSGALG00000053510"; refex1 "1:38214:39867:-"; refex2 "1:258715:258821:-";	can 1 dsmaxlg 0 adjref12 0
# 6091 (26 fields)

# output file
# 1	StringTie	intron	39868	258714	.	-	.	gene_id "MSTRG.1"; transcript_id "MSTRG.1.10"; refgn1 "ENSGALG00000042023"; refgn2 "ENSGALG00000053510"; refex1 "1:38214:39867:-"; refex2 "1:258715:258821:-";	can 1 dsmaxlg 0 adjref12 0 n-1to2nd 0
# 6091 (28 fields)


BEGIN{
    OFS="\t";
    # read the fileRef tsv file and remember the total nb of exons of each transcripts to which an exon belongs
    # and the biological position of this exon in each of those transcripts
    while (getline < fileRef >0)
    {
	nbexlist[$1]=$3;
	bioposlist[$1]=$5;
    }
}

{
    # found5p and found3p are the booleans that indicate whether the most 5' and the most 3' exons are indeed the n-1th and the 2nd exon of a transcript
    found5p=0;
    found3p=0;
    if($7=="+"||$7==".")
    {
	split($18,a5p,"\"");
	split($20,a3p,"\"");
    }
    else
    {
	if($7=="-")
	{
	    split($20,a5p,"\"");
	    split($18,a3p,"\"");
	}
    }
    split(nbexlist[a5p[2]],a15p,",");
    split(bioposlist[a5p[2]],a25p,",");
    split(bioposlist[a3p[2]],a23p,",");
    k=1;
    while(found5p==0&&a25p[k]!="")
    {
	if(a25p[k]==(a15p[k]-1))
	{
	    found5p=1;
	}
	k++;
    }
    k=1;
    while(found3p==0&&a23p[k]!="")
    {
	if(a23p[k]==2)
	{
	    found3p=1;
	}
	k++;
    }
    print $0" n-1to2nd "(found5p==1&&found3p==1 ? 1 : 0);
}
