# add_novelty_status_to_tagada_file.awk
# takes as input a tagada novel gtf file with indication of gene biotype and adds a novelty status
# for exon, transcript and gene rows
# - for exons and transcripts it is based on the presence of transcript_biotype from ensembl in the row
# and can only be known or unknown
# - for genes it is based on what we found for the different transcripts of the gene and can be
# known if all tr are known, unknown if all tr are unknown and mix otherwise

# example
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/geneswitch/results/rnaseq/sus_scrofa/TAGADA.v2.1.0.Sscrofa11.1.102.22-06-07/annotation
# time awk -f $pgm novel.with.gnbt.gtf > novel.with.gnbt.gnnov.gtf

# input
# # 12	tagada	gene	31695404	31723780	0	-	.	gene_id "LOC_000000010029"; ref_gene_id "."; feelnc_gene_biotype "mRNA";
# 12	tagada	transcript	31695404	31723729	0	-	.	gene_id "LOC_000000010029"; transcript_id "TM_000000214961"; contains "lung_stage3.filtered:STRG.74333.1,liver_stage2.filtered:STRG.57124.1,ileum_stage2.filtered:STRG.42019.1,cerebellum_stage2.filtered:STRG.30398.1,cerebellum_stage1.filtered:STRG.39652.1,lung_stage2.filtered:STRG.61223.1,liver_stage1.filtered:STRG.38572.1,cerebellum_stage3.filtered:STRG.32836.2,skin_stage2.filtered:STRG.43040.1,muscle_stage3.filtered:STRG.35479.1,lung_stage1.filtered:STRG.51252.1,skin_stage3.filtered:STRG.52500.1,muscle_stage1.filtered:STRG.39008.2,skin_stage1.filtered:STRG.51505.1,kidney_stage3.filtered:STRG.43454.1,kidney_stage1.filtered:STRG.38877.1,ileum_stage3.filtered:STRG.37374.1,muscle_stage2.filtered:STRG.59944.1,kidney_stage2.filtered:STRG.52524.1,ileum_stage1.filtered:STRG.48093.1,liver_stage3.filtered:STRG.61641.2"; contains_count "21"; 3p_dists_to_3p "0,154,414,431,438,406,428,447,441,439,425,439,464,439,435,438,440,446,475,472,475"; 5p_dists_to_5p "0,87,18,1,1,37,16,1,15,18,34,22,0,27,39,37,39,39,16,37,38"; flrpm "27.953205"; longest "lung_stage3.filtered:STRG.74333.1"; longest_FL_supporters "kidney_stage1.filtered:STRG.38877.1,skin_stage1.filtered:STRG.51505.1,muscle_stage3.filtered:STRG.35479.1,liver_stage1.filtered:STRG.38572.1,muscle_stage2.filtered:STRG.59944.1,kidney_stage3.filtered:STRG.43454.1,ileum_stage1.filtered:STRG.48093.1,ileum_stage2.filtered:STRG.42019.1,skin_stage3.filtered:STRG.52500.1,kidney_stage2.filtered:STRG.52524.1,skin_stage2.filtered:STRG.43040.1,liver_stage2.filtered:STRG.57124.1,ileum_stage3.filtered:STRG.37374.1,muscle_stage1.filtered:STRG.39008.2,liver_stage3.filtered:STRG.61641.2,lung_stage3.filtered:STRG.74333.1,cerebellum_stage2.filtered:STRG.30398.1,cerebellum_stage1.filtered:STRG.39652.1,lung_stage1.filtered:STRG.51252.1,cerebellum_stage3.filtered:STRG.32836.2,lung_stage2.filtered:STRG.61223.1"; longest_FL_supporters_count "21"; mature_RNA_length "3138"; meta_3p_dists_to_5p "1,0.950924155513066,0.868068833652008,0.862651370299554,0.860420650095602,0.870618228170809,0.863607393244105,0.85755258126195,0.859464627151052,0.860101975780752,0.864563416188655,0.860101975780752,0.852135117909496,0.860101975780752,0.861376673040153,0.860420650095602,0.859783301465902,0.857871255576801,0.848629700446144,0.849585723390695,0.848629700446144"; meta_5p_dists_to_5p "0,0.0277246653919694,0.00573613766730402,0.000318674314850223,0.000318674314850223,0.0117909496494583,0.00509878903760357,0.000318674314850223,0.00478011472275335,0.00573613766730402,0.0108349267049076,0.00701083492670491,0,0.00860420650095602,0.0124282982791587,0.0117909496494583,0.0124282982791587,0.0124282982791587,0.00509878903760357,0.0117909496494583,0.0121096239643085"; rpm "27.953205"; spliced "1"; ref_gene_id "."; feelnc_gene_biotype "mRNA";
# 51171 (14 fields)
# 931554 (16 fields)
# 141602 (18 fields)
# 592967 (20 fields)
# 30530 (22 fields)
# 87739 (42 fields)
# 48935 (44 fields)
# 50488 (46 fields)
# 10234 (48 fields)

# output
# 12	tagada	transcript	31695404	31723729	0	-	.	gene_id "LOC_000000010029"; transcript_id "TM_000000214961"; contains "lung_stage3.filtered:STRG.74333.1,liver_stage2.filtered:STRG.57124.1,ileum_stage2.filtered:STRG.42019.1,cerebellum_stage2.filtered:STRG.30398.1,cerebellum_stage1.filtered:STRG.39652.1,lung_stage2.filtered:STRG.61223.1,liver_stage1.filtered:STRG.38572.1,cerebellum_stage3.filtered:STRG.32836.2,skin_stage2.filtered:STRG.43040.1,muscle_stage3.filtered:STRG.35479.1,lung_stage1.filtered:STRG.51252.1,skin_stage3.filtered:STRG.52500.1,muscle_stage1.filtered:STRG.39008.2,skin_stage1.filtered:STRG.51505.1,kidney_stage3.filtered:STRG.43454.1,kidney_stage1.filtered:STRG.38877.1,ileum_stage3.filtered:STRG.37374.1,muscle_stage2.filtered:STRG.59944.1,kidney_stage2.filtered:STRG.52524.1,ileum_stage1.filtered:STRG.48093.1,liver_stage3.filtered:STRG.61641.2"; contains_count "21"; 3p_dists_to_3p "0,154,414,431,438,406,428,447,441,439,425,439,464,439,435,438,440,446,475,472,475"; 5p_dists_to_5p "0,87,18,1,1,37,16,1,15,18,34,22,0,27,39,37,39,39,16,37,38"; flrpm "27.953205"; longest "lung_stage3.filtered:STRG.74333.1"; longest_FL_supporters "kidney_stage1.filtered:STRG.38877.1,skin_stage1.filtered:STRG.51505.1,muscle_stage3.filtered:STRG.35479.1,liver_stage1.filtered:STRG.38572.1,muscle_stage2.filtered:STRG.59944.1,kidney_stage3.filtered:STRG.43454.1,ileum_stage1.filtered:STRG.48093.1,ileum_stage2.filtered:STRG.42019.1,skin_stage3.filtered:STRG.52500.1,kidney_stage2.filtered:STRG.52524.1,skin_stage2.filtered:STRG.43040.1,liver_stage2.filtered:STRG.57124.1,ileum_stage3.filtered:STRG.37374.1,muscle_stage1.filtered:STRG.39008.2,liver_stage3.filtered:STRG.61641.2,lung_stage3.filtered:STRG.74333.1,cerebellum_stage2.filtered:STRG.30398.1,cerebellum_stage1.filtered:STRG.39652.1,lung_stage1.filtered:STRG.51252.1,cerebellum_stage3.filtered:STRG.32836.2,lung_stage2.filtered:STRG.61223.1"; longest_FL_supporters_count "21"; mature_RNA_length "3138"; meta_3p_dists_to_5p "1,0.950924155513066,0.868068833652008,0.862651370299554,0.860420650095602,0.870618228170809,0.863607393244105,0.85755258126195,0.859464627151052,0.860101975780752,0.864563416188655,0.860101975780752,0.852135117909496,0.860101975780752,0.861376673040153,0.860420650095602,0.859783301465902,0.857871255576801,0.848629700446144,0.849585723390695,0.848629700446144"; meta_5p_dists_to_5p "0,0.0277246653919694,0.00573613766730402,0.000318674314850223,0.000318674314850223,0.0117909496494583,0.00509878903760357,0.000318674314850223,0.00478011472275335,0.00573613766730402,0.0108349267049076,0.00701083492670491,0,0.00860420650095602,0.0124282982791587,0.0117909496494583,0.0124282982791587,0.0124282982791587,0.00509878903760357,0.0117909496494583,0.0121096239643085"; rpm "27.953205"; spliced "1"; ref_gene_id "."; feelnc_gene_biotype "mRNA"; novelty "unknown";
# 12	tagada	exon	31695404	31697761	0	-	.	gene_id "LOC_000000010029"; transcript_id "TM_000000214961"; ref_gene_id "."; feelnc_gene_biotype "mRNA"; novelty "unknown";
# 51171 (16 fields)
# 931554 (18 fields)
# 141602 (20 fields)
# 592967 (22 fields)
# 30530 (24 fields)
# 87739 (44 fields)
# 48935 (46 fields)
# 50488 (48 fields)
# 10234 (50 fields)

BEGIN{OFS="\t"}

{
    if($3=="exon"||$3=="transcript")
    {
	if($0~/transcript_biotype/)
	{
	    known[$12]=1;
	    print $0" novelty \"known\"\;";
	}
	else
	{
	    known[$12]=0;
	    print $0" novelty \"unknown\"\;";
	}
	if($3=="transcript")
	{
	    novlist[$10]=(novlist[$10])(known[$12])(",");
	}
    }
    else
    {
	if($3=="gene")
	{
	    row[$10]=$0;
	}
    }
}

END{
    for(g in row)
    {
	stat=((novlist[g]~/1/) ? ((novlist[g]!~/0/) ? "known" : "mix") : "unknown");
	print row[g]" novelty \""stat"\"\;";
    }
}
