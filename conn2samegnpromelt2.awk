# conn2samegnpromelt2.awk
# Given a bedpe file of prom to elt2 connections, with indication of list of genes of prom in field fld as comma separated list of
# chr:gbeg:gend:str:geneid:genename, and indication in the last 8 fields of 8 comma separated lists of gene ids (the ones whose exons
# are overlapped by the elt2, the ones whose introns encompass the elt2, the ones whose tss are overlapped by the elt2, same for tss+-1kb,
# same for tss+-5kb, and then same for tts), and provides the following summary stats in a single row and separated by tabs
# - total number of connections in the bedpe file
# - the number and % of total connections for which at least one gene of the prom gene list is identical to one of the genes where the elt2 lies
#   (therefore at most 5kb from a gene)

# example
# dir=~/regenet/workspace/sdjebali/egprediction/from3D/predictions/capturehic/jung.ren.2019
# pgm=~/fragencode/tools/multi/Scripts/conn2samegnpromelt2.awk
# ct=HCmerge
# e=elt2A
# cd $dir/$ct/score2/pall
# time zcat $ct.pall.score2.gninfo.elt2info.bedpe.gz | awk -v e=$e '$24==e' | awk -v fld=21 -f $pgm 
# real	0m0.905s

# input
# chr1	1183838	1209657	chr1	1491937	1496017	chr1:1183838:1209657,chr1:1491937:1496017	2.18497234006452	.	.	po	22511	4859	6	297230	2.0085	2.9873	1.0369	1.9571	1	chr1:1189288:1209265:-:ENSG00000160087.16:UBE2J2,	0	NA	elt2A	282280	NA	ENSG00000160075.10,	NA	NA	NA	NA	NA	ENSG00000160075.10,
# 64514 (33 fields)
# chr1:1853395:1935276:-:ENSG00000142609.13:C1orf222,chr1:1950779:1962192:+:ENSG00000187730.6:GABRD,	ENSG00000142609.13,

# output
# totconn	samegnpromelt2.nb	samegnpromelt2.pcent
# 56716	11968	21.1016


BEGIN{
    OFS="\t"
}

{
    tot++;
    split($fld,a,",");
    k=1;
    ok=0;
    while(ok==0&&a[k]!="")   # going over the list of genes of the prom
    {
	split(a[k],a1,":");   # a1[5] should be the current gene id of the prom
	i=1;
	while(ok==0&&i<=8)   # going over each type of gx overlap between the elt2 and the genome (exon, intron, ..., tts+-5kb), located in field $(NF-8+i)
	{
	    split($(NF-8+i),b,",");
	    j=1;
	    while(ok==0&&b[j]!="")   # going over the list of genes of the field $(NF-8+i), the current gene id being in b[j]
	    {
		if(a1[5]==b[j])
		{
		    ok=1;
		}
		j++;
	    }
	    i++;
	}
	k++;
    }
    same+=ok;
}

END{
    print "totconn", "samegnpromelt2.nb", "samegnpromelt2.pcent";
    print tot, same, same/tot*100;
}

