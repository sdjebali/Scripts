# assign_ens_tid_gid_to_retrocopies.awk
# takes as input the output of bed intersect on retrocopy exon file with transcript_id as 3rd key in 9th field (14th column)
# and ens exon gff file. By default the transcript_id and gene_id of the ensembl exon are in columns 32 and 28
# but this could be changed with fldtr and fldgn
# note: could be made even more general allowing the retrocopy transcript id to be elsewhere than in 14th column

# usage
# cd ~/fragencode/workspace/sdjebali/geneswitch/retrogliss
# pgm=~/fragencode/tools/multi/Scripts/assign_ens_tid_gid_to_retrocopies.awk
# time awk -v fldtr=32 -v fldgn=28 -f $pgm intersect_rtc_exon_vs_ens_exon.tsv > rtcid.ens.trid.gnid.overlap.tsv

# main input = intersect_rtc_exon_vs_ens_exon.tsv
# 1	ensembl	exon	900528	900875	56.78	+	0	gene_id "RetroG_1"; gene_biotype "processed_pseudogene"; transcript_id "RetroT_1"; transcript_biotype "processed_pseudogene"; exon_number "1"; 	.	.	.	-1	-1	.	.	.	.	0
# 1	ensembl	exon	1430846	1432510	51.53	-	0	gene_id "RetroG_2"; gene_biotype "processed_pseudogene"; transcript_id "RetroT_2"; transcript_biotype "processed_pseudogene"; exon_number "1"; 	1	ensembl	exon	1430447	1432730	.	-	.	gene_id "ENSBTAG00000012594"; gene_version "7"; transcript_id "ENSBTAT00000110316"; transcript_version "1"; exon_number "1"; gene_name "MRPS6"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "MRPS6-202"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSBTAE00000602654"; exon_version "1"; tag "Ensembl_canonical";	1665
# 1126 (28 fields)
# 158 (49 fields)
# 1045 (51 fields)
# 336 (53 fields)
# 360 (55 fields)
# 30 (57 fields)

# main output = rtcid.ens.trid.gnid.overlap.tsv
# RetroT_80	ENSBTAT00000095539	ENSBTAG00000069459	197
# RetroT_40	ENSBTAT00000074364	ENSBTAG00000050699	435
# 1183 (4 fields)

BEGIN{
    if(fldtr=="")
    {
	fldtr=32;
    }
    if(fldgn=="")
    {
	fldgn=28;
    }
}


$NF!=0{
    split($14,a,"\"");
    split($fldtr,b,"\"");
    split($fldgn,c,"\"");
    gnid[b[2]]=c[2];
    rtok[a[2]]=1;
    enstrok[a[2],b[2]]++;
    if(enstrok[a[2],b[2]]==1)
    {
	nbenstrid[a[2]]++;
	enstrid[a[2],nbenstrid[a[2]]]=b[2];
    }
    lg[a[2],b[2]]+=$NF;
}

END{
    OFS="\t";
    for(r in rtok)
    {
	l=0;
	etrid="NA";
	for(i=1; i<=nbenstrid[r]; i++)
	{
	    if(lg[r,enstrid[r,i]]>l)
	    {
		l=lg[r,enstrid[r,i]];
		etrid=enstrid[r,i];
	    }
	}
	print r, etrid, gnid[etrid], l;
    }
}
