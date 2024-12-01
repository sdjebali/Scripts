# addcdscoord.and.chimjuncincds.awk
# takes as input a tsv file that has tagiso exon-exon junctions identified as chimeric wrt ens108
# with 11 columns which are
# - exon1 coord
# - exon2 coord
# - tagiso trid
# - tagiso gnid
# - exon1 ens108 gnid
# - exon2 ens108 gnid
# - chimeric junction coord (wrt ens108, be careful for - strand the 1st coord is higher than the 2nd)
# - a boolean if chimpipe also found this chimeric junction
# - a boolean whether the tr is coding acc to feelnc
# - same with transdecoder
# - same with final (coding by feelnc and transdecoder)
# and as fileRef the complete tagiso annotation with cds rows and produces the same file as the input tsv
# but with 2 additional fields
# - the coordinates of the cds of the tagiso tr when it is considered coding by transdecoder (NA otherwise)
# - whether the chimeric junction is included in the CDS (NA if no cds)

# Example
# indir=/work/project/fragencode/workspace/geneswitch/results/rnaseq
# sp=gallus_gallus
# cd $indir/$sp/tagiso/annotation/compare.to.chimpipe
# pgm=~/fragencode/tools/multi/Scripts/addcdscoord.and.chimjuncincds.awk
# awk -v fileRef=$indir/$sp/nf-feelnc.tagiso.23-03-03/results/annotation/novel.clean.full.gff -f $pgm novel.exon-exon.with.ens108gnoverex.singlegneachside.diff.incp.cod.feelnc.TD.final.tsv > novel.exon-exon.with.ens108gnoverex.singlegneachside.diff.incp.cod.feelnc.TD.final.cdscoord.chimjuncincl.tsv
  
# fileRef=$indir/$sp/nf-feelnc.tagiso.23-03-03/results/annotation/novel.clean.full.gff
# 1	transdecoder	CDS	5315	5512	.	-	0	gene_id "LOC_000000158014"; transcript_id "TM_000000000012"; feelnc_transcript_biotype "mRNA"; feelnc_gene_biotype "mRNA"; transdecoder_gene_biotype "mRNA"; final_gene_biotype "mRNA"; transdecoder_transcript_biotype "mRNA"; final_transcript_biotype "mRNA";
# 14764 (20 fields)
# ...
# 122339 (60 fields)

# main input = novel.exon-exon.with.ens108gnoverex.singlegneachside.diff.incp.cod.feelnc.TD.final.tsv
# 6:6457020:6458896:-	6:6474895:6475025:-	G28850.2	LOC_000000063038	ENSGALG00000063938	ENSGALG00000049233	6_6474895_-:6_6458896_-	0	00	0
# 5311 (11 fields)

# output
# 6:6457020:6458896:-	6:6474895:6475025:-	G28850.2	LOC_000000063038	ENSGALG00000063938	ENSGALG00000049233	6_6474895_-:6_6458896_-	0	00	0	NA	NA
# 1:60042164:60042244:+	1:60062554:60062842:+	TM_000000108024	LOC_000000021857	ENSGALG00000067512	ENSGALG00000012975	1_60042244_+:1_60062554_+	11	1	1	60069270-60099039	0
# 5311 (13 fields)

BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	if($3=="CDS")
	{
	    cod[$12]=1;
	    if(cdsbeg[$12]==""||$4<cdsbeg[$12])
	    {
		cdsbeg[$12]=$4;
	    }
	    if(cdsend[$12]==""||$5>cdsend[$12])
	    {
		cdsend[$12]=$5;
	    }
	}
    }
}

{
    # 6_6474895_-:6_6458896_-
    split($7,a,":");
    split(a[1],a11,"_");
    split(a[2],a12,"_");
    chbeg=((a11[2]<a12[2]) ? a11[2] : a12[2]);
    chend=((a11[2]<a12[2]) ? a12[2] : a11[2]);
    t="\""$3"\"\;";
    cb=cdsbeg[t];
    ce=cdsend[t];
    cdscoord=(cod[t]==1 ? (cb"-"ce) : "NA");
    print $0, cdscoord, (cod[t]==1 ? included(chbeg,chend,cb,ce) : "NA");
}

function included(x,y,z,w){
    return (x>=z&&y<=w ? 1 : 0);
}
