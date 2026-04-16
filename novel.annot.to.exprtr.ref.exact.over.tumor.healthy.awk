# novel.annot.to.exprtr.ref.exact.over.tumor.healthy.awk
# similar to novel.annot.to.exprtr.novelty.nbex.struct.nbsamp.awk because it takes the same two inputs, which are
# - fileRef is a complete annotation file from tagada that was run with gencode annot on human and with tr id in $12 for exon and transcript rows
# - main input file is a 3 column tsv file without header that has all tr present in this annot and as columns
#   * trid
#   * gene id
#   * whether it overlaps a gencode tr in an exonic and stranded way (boolean)
# but here it returns for all expressed transcripts (meaning it takes out the ref tr that are not associated to any sample)
# - trid
# - whether it exactly corresponds to a ref tr (has ref in longest_FL_supporters field)
# - whether it overlaps exons of a ref tr on the same strand but does not exactly corresponds to a ref tr
# - whether it was seen in the lms tumor samples (if LMSXRY in the contains field)
# - whether it was seen in the lms healthy samples (if LMSXRNY in the contains field)

# example of usage on the genotoul cluster
##########################################
# srun --x11 --mem=16G --time=10:00:00 --pty bash
# cd ~/bridge/workspace/sdjebali/LMS/novel.annot/tumor.healthy
# annot=~/bridge/workspace/clatour/LMS_RNAseq/03_tagada_annotation/LMS_tagadaAnnot_batch003_20260403/annotation/novel.gtf
# pgm=~/fragencode/tools/multi/Scripts/novel.annot.to.exprtr.ref.exact.over.tumor.healthy.awk
# time awk -v fileRef=$annot -f $pgm compare.new.with.gencv49/novel.trid.gnid.trovergencv49.tsv > lms.and.genc.exprtr.id.ref.exact.over.tumor.healthy.tsv
# real	0m59.301s


# fileRef = $annot = ~/bridge/workspace/clatour/LMS_RNAseq/03_tagada_annotation/LMS_tagadaAnnot_batch003_20260403/annotation/novel.gtf is like this
# chr1	tagada	transcript	11121	14413	0	+	.	gene_id "LOC_000000033196"; transcript_id "ENST00000832824.1"; contains "ref:ENST00000832824.1"; contains_count "1"; 3p_dists_to_3p "0"; 5p_dists_to_5p "0"; flrpm "2.247738"; longest "ref:ENST00000832824.1"; longest_FL_supporters "ref:ENST00000832824.1"; longest_FL_supporters_count "1"; mature_RNA_length "1379"; meta_3p_dists_to_5p "1"; meta_5p_dists_to_5p "0"; rpm "2.247738"; spliced "1"; ref_gene_id "ENSG00000290825.2"; tmerge_tr_id "TM_000000001418"; transcript_type "lncRNA";
# GL000009.2	tagada	transcript	49627	97288	0	-	.	gene_id "LOC_000000123428"; transcript_id "TM_000000000041"; contains "LMS9R1.filtered:STRG.46590.1,LMS1R.filtered:STRG.39975.1,LMS8R1.filtered:STRG.44674.1,LMS7R1.filtered:STRG.39712.1,LMS1R.filtered:STRG.39974.1,LMS8R1.filtered:STRG.44675.1,LMS7R1.filtered:STRG.39713.1,LMS8R1.filtered:STRG.44676.1,LMS7R1.filtered:STRG.39714.1,LMS8RN.filtered:STRG.35182.1"; contains_count "10"; 3p_dists_to_3p "0,5391,39856,39860,522,45217,46284,46290,47251,46401"; 5p_dists_to_5p "23539,0,1628,1628,42592,1661,1051,1060,143,1054"; flrpm "2.247738"; longest "LMS9R1.filtered:STRG.46590.1"; longest_FL_supporters "LMS9R1.filtered:STRG.46590.1"; longest_FL_supporters_count "1"; mature_RNA_length "47662"; meta_3p_dists_to_5p "1,0.886891024296085,0.163778271998657,0.163694347698376,0.98904787881331,0.0512987285468507,0.0289119214468549,0.0287860349964332,0.00862322185388775,0.0264571356636314"; meta_5p_dists_to_5p "0.493873526079476,0,0.0341571902144266,0.0341571902144266,0.893625949393647,0.034849565691746,0.0220511098988712,0.0222399395745038,0.00300029373505098,0.0221140531240821"; rpm "22.47738"; spliced "0"; ref_gene_id ".";
# 90817 (12 fields)
# 1448970 (14 fields)
# 3469181 (18 fields)
# 205789 (40 fields)
# 474171 (44 fields) 

# main input file = compare.new.with.gencv49/novel.trid.gnid.trovergencv49.tsv is like this
# ENST00000966617.1	LOC_000000006627	1
# TM_000001021931	LOC_000000071825	0
# 679960 (3 fields)

# output file = lms.and.genc.exprtr.id.ref.exact.over.tumor.healthy.tsv
# trid	refexact	refovernotexact	tumor	healthy
# TM_000001021931	0	0	0	1
# 395437 (5 fields) 

BEGIN{
    OFS="\t";
    print "trid", "refexact", "refovernotexact", "tumor", "healthy";
    while (getline < fileRef >0)
    {
	split($0,a,"\t");
	split(a[9],b,"; ");
	split($12,c,"\"");
	if(a[3]=="transcript")
	{
	    k=1;
	    while(b[k]!="")
	    {
		split(b[k],d," ");
		# we look at the contains field to see if it detected in lms or healthy samples or both
		if(d[1]=="contains")
		{
		    split(d[2],e,"\"");
		    split(e[2],f,",");
		    l=1;
		    while(f[l]!="")  # f[l] is like this LMS9R1.filtered:STRG.46590.1 or LMS8RN.filtered:STRG.35182.1 or ref:ENST00000832828.1
		    {
			split(f[l],g,".");
			seen[c[2],g[1]]++;
			if(g[1]!~/ref/&&seen[c[2],g[1]]==1)
			{
			    nbsamp[c[2]]++;
			    if(g[1]~/N/)
			    {
				health[c[2]]=1;
			    }
			    else
			    {
				tumor[c[2]]=1;
			    }
			}
			l++;
		    }
		}
		else
		{
		    # we look at the longest_FL_supporters field to see if it is exactly a gencode tr, in which case it would have a ref in this field
		    if(d[1]=="longest_FL_supporters")
		    {
			if(d[2]~/ref/)
			{
			    refexact[c[2]]=1;
			}
		    }
		}
		k++;
	    }
	}
    }
}


# nbsamp[$1]>=1 means it is seen in at least one lms sample
# we want to print for those expr tr, the following info
# - trid
# - whether it exactly corresponds to a ref tr (has ref in longest_FL_supporters field)
# - whether it overlaps exons of a ref tr on the same strand but does not exactly corresponds to a ref tr ($3 from our main input file but not in 1st class)
# - whether it was seen in the lms tumor samples (if LMSXRY in contains field)
# - whether it was seen in the lms healthy samples (if LMSXRNY in contains field)
nbsamp[$1]>=1{
    print $1, (refexact[$1]==1 ? 1 : 0), (($3==1&&refexact[$1]!=1) ? 1 : 0), (tumor[$1]==1 ? 1 : 0), (health[$1]==1 ? 1 : 0);
}
