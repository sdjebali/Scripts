# novel.annot.to.exprtr.novelty.nbex.struct.nbsamp.awk
# this script takes as input a 
# - fileRef is a complete annotation file from tagada that was run with gencode annot on human and with tr id in $12 for exon and transcript rows
# - main input file is a 3 column tsv file without header that has all tr present in this annot and as columns
#   * trid
#   * gene id
#   * whether it overlaps a gencode tr in an exonic and stranded way (boolean)
# and returns as output a tsv file with header that only has transcripts from fileRef that are expressed (seen in at least one sample)
# and that are either exactly Gencode tr or not exonically and strandedly overlapping gencode tr
# - trid
# - tr class
# - nbex
# - structure class
# - nb samples where it is seen according to longest_FL_supporters field

# example of usage on the genotoul cluster
##########################################
# srun --x11 --mem=16G --time=10:00:00 --pty bash
# cd ~/work/sarcomas/novel.annot/compare.new.with.gencv49
# annot=~/bridge/workspace/clatour/LMS_RNAseq/05_merged_annotation/LMS_annotMergedQC_batch001_20260313/LMS.combined.trclass.gtf
# pgm=~/fragencode/tools/multi/Scripts/novel.annot.to.exprtr.novelty.nbex.struct.nbsamp.awk
# time awk -v fileRef=$annot -f $pgm novel.trid.gnid.trovergencv49.tsv > lms.and.genc.exprtr.id.noveltyclass.nbex.structureclass.nbsamp.tsv
# real	2m39.326s  *** takes long

# fileRef= $annot (~/bridge/workspace/clatour/LMS_RNAseq/05_merged_annotation/LMS_annotMergedQC_batch001_20260313/LMS.combined.trclass.gtf) is like this
# GL000008.2	tagada	exon	83353	83545	0	+	.	gene_id "LOC_000000044264"; transcript_id "TM_000000000001"; ref_gene_id ".";
# GL000008.2	tagada	transcript	83353	89823	0	+	.	gene_id "LOC_000000044264"; transcript_id "TM_000000000001"; contains "LMS88R.filtered:STRG.38525.1,LMS92R.filtered:STRG.38499.1"; contains_count "2"; 3p_dists_to_3p "0,0"; 5p_dists_to_5p "0,26"; flrpm "33.820702"; longest "LMS88R.filtered:STRG.38525.1"; longest_FL_supporters "LMS92R.filtered:STRG.38499.1,LMS88R.filtered:STRG.38525.1"; longest_FL_supporters_count "2"; mature_RNA_length "4450"; meta_3p_dists_to_5p "1,1"; meta_5p_dists_to_5p "0,0.00584269662921348"; rpm "33.820702"; spliced "1"; ref_gene_id "."; trclass "Intergenic_or_antisense:novel";
# 169678 (12 fields)
# 7625023 (14 fields)
# 3300939 (18 fields)
# 1117434 (42 fields)
# 447234 (46 fields)

# main input file = novel.trid.gnid.trovergencv49.tsv
# TM_000003078957	LOC_000000006654	1
# ENST00000607970.3	LOC_000000003162	1
# 1564668 (3 fields)

# output file = lms.and.genc.exprtr.id.noveltyclass.nbex.structureclass.nbsamp.tsv is like this
# trid	trclass	nbex	structure	nbsamples
# ENST00000607970.3	2.Genc	4	1.Spliced	3
# 480136 (5 fields)


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

# $3==0 means it does not overlap a gencode exon on the same strand
# $1~/ENST/ means it exactly corresponds to a gencode tr
# nbsamp[$1]>=1 means it is seen in at least one lms sample
($3==0||$1~/ENST/)&&nbsamp[$1]>=1{
    print $1, ($3==0 ? "1.LMS" : ($1~/ENST/ ? "2.Genc" : "NA")), nbex[$1], (nbex[$1]==1 ? "2.Monoexonic" : "1.Spliced"), nbsamp[$1];
}
