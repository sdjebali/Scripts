# novel.annot.to.exprtr.novelty.nbex.struct.nbsamp.4trclasses.awk
# same as novel.annot.to.exprtr.novelty.nbex.struct.nbsamp.awk except that
# - it allows 4 tr classes (1.tumor.igas, 2.tumor.over, 3.healthy.igas, 4.healthy.over) and not just 2 (Genc for exact Genc or LMS when not over Genc exons)
# - it takes the tr class from the main input file as column 3 (and not a simple boolean saying over gencode or not) and does not guess it 
# this script takes as input a 
# - fileRef is a complete annotation file from tagada that was run with gencode annot on human and with tr id in $12 for exon and transcript rows
# - main input file is a 4 column tsv file without header that has all tr present in this annot and as columns
#   * trid
#   * gene id
#   * its class out of 4 
# and returns as output a tsv file with header that only has transcripts from fileRef that are expressed (seen in at least one sample)
# and that are either Gencode tr or not exonically and strandedly overlapping gencode tr
# - trid
# - tr class out of 4
# - nbex
# - structure class
# - nb samples where it is seen according to longest_FL_supporters field

# example of usage on the genotoul cluster
##########################################
# srun --x11 --mem=16G --time=10:00:00 --pty bash
# cd ~/bridge/workspace/sdjebali/LMS/novel.annot/tumor.healthy
# annot=~/bridge/workspace/clatour/LMS_RNAseq/05_merged_annotation//LMS_annotMergedQC_batch003_20260417/LMS.combined.trclass.gtf
# pgm=~/fragencode/tools/multi/Scripts/novel.annot.to.exprtr.novelty.nbex.struct.nbsamp.4trclasses.awk
# time awk -v fileRef=$annot -f $pgm tr.id.gnid.class.tsv > tr.id.noveltyclass.nbex.structureclass.nbsamp.tsv
# real	0m53.994s 

# fileRef= $annot (~/bridge/workspace/clatour/LMS_RNAseq/05_merged_annotation//LMS_annotMergedQC_batch003_20260417/LMS.combined.trclass.gtf)
# GL000008.2	tagada	exon	83354	83545	0	+	.	gene_id "LOC_000000194708"; transcript_id "TM_000000000001"; ref_gene_id ".";
# GL000008.2	tagada	exon	83354	84014	0	+	.	gene_id "LOC_000000194708"; transcript_id "TM_000000000002"; ref_gene_id ".";
# 90817 (12 fields)
# 1448970 (14 fields)
# 3469181 (18 fields)
# 205789 (42 fields)
# 474171 (46 fields)  

# main input file = tr.id.gnid.class.tsv
# TM_000001021932	LOC_000000071825	1.tumor.igas
# TM_000001997780	LOC_000000079730	1.tumor.igas
# 160085 (3 fields)

# output file = tr.id.noveltyclass.nbex.structureclass.nbsamp.tsv
# trid	trclass	nbex	structure	nbsamples
# TM_000001021932	1.tumor.igas	1	2.Monoexonic	1
# 160086 (5 fields) 

# !!! be careful a gencode transcript can be covered by several partial lms transcripts and therefore we need to increment the nb of samples !!!
# !!! only once for each transcript and each sample seen in 
BEGIN{
    OFS="\t";
    print "trid", "trclass", "nbex", "structure", "nbsamples";
    while (getline < fileRef >0)
    {
	split($0,a,"\t");
	split(a[9],b,"; ");
	split($12,c,"\"");
	if(a[3]=="exon")
	{
	    nbex[c[2]]++;
	}
	else
	{
	    if(a[3]=="transcript")
	    {
		k=1;
		found=0;
		while(found==0&&b[k]!="")
		{
		    split(b[k],d," ");
		    if(d[1]=="longest_FL_supporters")
		    {
			found=1;
			split(d[2],e,"\"");
			split(e[2],f,",");
			l=1;
			while(f[l]!="")  # f[l] is like this LMS33R.filtered:STRG.5060.1 or ref:ENST00000613695.1
			{
			    split(f[l],g,".");
			    seen[c[2],g[1]]++;
			    if(g[1]!~/ref/&&seen[c[2],g[1]]==1)
			    {
				nbsamp[c[2]]++;
			    }
			    l++;
			}
		    }
		    k++;
		}
	    }
	}
    }
}

# nbsamp[$1]>=1 means it is seen in at least one sample but in theory not needed here since we only took the expressed tr
nbsamp[$1]>=1{
    print $1, $3, nbex[$1], (nbex[$1]==1 ? "2.Monoexonic" : "1.Spliced"), nbsamp[$1];
}
