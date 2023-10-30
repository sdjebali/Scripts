# clean_gtf_after2feelncruns.awk
# this scripts takes as input a gtf file for which feelnc has been run twice and
# - for exon and transcript rows when there are two feelnc tags then only keep the second one
# - for gene rows remove all the feelnc tags (feelnc_biotype) 

# srun --x11 --mem=8G --pty bash
# dir=~/fragencode/workspace/geneswitch/results/rnaseq
# sp=sus_scrofa
# cd $dir/$sp/nf-feelnc.tagiso.23-03-03/results/annotation
# pgm=~/fragencode/tools/multi/Scripts/clean_gtf_after2feelncruns.awk
# time awk -f $pgm novel.gtf > novel.clean.gtf
# real	0m41.542s  *** same number of rows as the input file and visually checked ok

# input novel.gtf
# 1	tagada	exon	1	961	0	+	.	gene_id "LOC_000000075193"; transcript_id "G1.1"; ref_gene_id "G1"; tmerge_tr_id "TM_000000000001"; transcript_biotype ""; feelnc_biotype "lncRNA"; ensembl_gene_id "ENSSSCG00000048769"; best_ref "ENSSSCG00000048769"; feelnc_biotype "lncRNA";
# 1	tagada	gene	1	3862	0	+	.	gene_id "LOC_000000075193"; ref_gene_id "G1"; ensembl_gene_id "ENSSSCG00000048769"; best_ref "ENSSSCG00000048769"; feelnc_biotype "lncRNA"; feelnc_biotype "mRNA"; feelnc_biotype "noORF"; feelnc_biotype "TUCp";
# 1	tagada	transcript	1	3862	0	+	.	gene_id "LOC_000000075193"; transcript_id "G1.1"; contains "skin_stage2.filtered:STRG.1.1,cerebellum_stage3.filtered:STRG.1.1,kidney_stage3.filtered:STRG.1.1,cerebellum_stage1.filtered:STRG.1.1,skin_stage1.filtered:STRG.1.1,kidney_stage2.filtered:STRG.1.1,ileum_stage2.filtered:STRG.1.1,kidney_stage1.filtered:STRG.1.1,lung_stage1.filtered:STRG.1.1,skin_stage3.filtered:STRG.1.1,liver_stage3.filtered:STRG.1.1,liver_stage2.filtered:STRG.1.1,liver_stage1.filtered:STRG.1.1,muscle_stage1.filtered:STRG.1.1,ileum_stage1.filtered:STRG.1.1,lung_stage2.filtered:STRG.1.1,muscle_stage2.filtered:STRG.1.1,ileum_stage3.filtered:STRG.1.1,lung_stage3.filtered:STRG.1.1,cerebellum_stage2.filtered:STRG.1.1,muscle_stage3.filtered:STRG.1.1,ref:G1.1"; contains_count "22"; 3p_dists_to_3p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; 5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; flrpm "58.872946"; longest "skin_stage2.filtered:STRG.1.1,cerebellum_stage3.filtered:STRG.1.1,kidney_stage3.filtered:STRG.1.1,cerebellum_stage1.filtered:STRG.1.1,skin_stage1.filtered:STRG.1.1,kidney_stage2.filtered:STRG.1.1,ileum_stage2.filtered:STRG.1.1,kidney_stage1.filtered:STRG.1.1,lung_stage1.filtered:STRG.1.1,skin_stage3.filtered:STRG.1.1,liver_stage3.filtered:STRG.1.1,liver_stage2.filtered:STRG.1.1,liver_stage1.filtered:STRG.1.1,muscle_stage1.filtered:STRG.1.1,ileum_stage1.filtered:STRG.1.1,lung_stage2.filtered:STRG.1.1,muscle_stage2.filtered:STRG.1.1,ileum_stage3.filtered:STRG.1.1,lung_stage3.filtered:STRG.1.1,cerebellum_stage2.filtered:STRG.1.1,muscle_stage3.filtered:STRG.1.1,ref:G1.1"; longest_FL_supporters "liver_stage1.filtered:STRG.1.1,liver_stage3.filtered:STRG.1.1,muscle_stage2.filtered:STRG.1.1,ileum_stage3.filtered:STRG.1.1,muscle_stage3.filtered:STRG.1.1,cerebellum_stage2.filtered:STRG.1.1,cerebellum_stage3.filtered:STRG.1.1,cerebellum_stage1.filtered:STRG.1.1,skin_stage3.filtered:STRG.1.1,lung_stage1.filtered:STRG.1.1,kidney_stage1.filtered:STRG.1.1,lung_stage2.filtered:STRG.1.1,muscle_stage1.filtered:STRG.1.1,ileum_stage1.filtered:STRG.1.1,liver_stage2.filtered:STRG.1.1,ref:G1.1,lung_stage3.filtered:STRG.1.1,skin_stage2.filtered:STRG.1.1,kidney_stage2.filtered:STRG.1.1,ileum_stage2.filtered:STRG.1.1,skin_stage1.filtered:STRG.1.1,kidney_stage3.filtered:STRG.1.1"; longest_FL_supporters_count "22"; mature_RNA_length "1801"; meta_3p_dists_to_5p "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1"; meta_5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; rpm "58.872946"; spliced "1"; ref_gene_id "G1"; tmerge_tr_id "TM_000000000001"; transcript_biotype ""; feelnc_biotype "lncRNA"; ensembl_gene_id "ENSSSCG00000048769"; best_ref "ENSSSCG00000048769"; feelnc_biotype "lncRNA";
# 5241 (18 fields)
# 902361 (20 fields)
# 71444 (22 fields)
# 1247357 (24 fields)
# 22682 (26 fields)
# 3767 (44 fields)
# 101812 (46 fields)
# 16930 (48 fields)
# 141100 (50 fields)
# 5025 (52 fields)


# output novel.clean.gtf
# 1	tagada	exon	1	961	0	+	.	gene_id "LOC_000000075193"; transcript_id "G1.1"; ref_gene_id "G1"; tmerge_tr_id "TM_000000000001"; transcript_biotype ""; ensembl_gene_id "ENSSSCG00000048769"; best_ref "ENSSSCG00000048769"; feelnc_biotype "lncRNA";
# 1	tagada	gene	1	3862	0	+	.	gene_id "LOC_000000075193"; ref_gene_id "G1"; ensembl_gene_id "ENSSSCG00000048769"; best_ref "ENSSSCG00000048769";
# 1	tagada	transcript	1	3862	0	+	.	gene_id "LOC_000000075193"; transcript_id "G1.1"; contains "skin_stage2.filtered:STRG.1.1,cerebellum_stage3.filtered:STRG.1.1,kidney_stage3.filtered:STRG.1.1,cerebellum_stage1.filtered:STRG.1.1,skin_stage1.filtered:STRG.1.1,kidney_stage2.filtered:STRG.1.1,ileum_stage2.filtered:STRG.1.1,kidney_stage1.filtered:STRG.1.1,lung_stage1.filtered:STRG.1.1,skin_stage3.filtered:STRG.1.1,liver_stage3.filtered:STRG.1.1,liver_stage2.filtered:STRG.1.1,liver_stage1.filtered:STRG.1.1,muscle_stage1.filtered:STRG.1.1,ileum_stage1.filtered:STRG.1.1,lung_stage2.filtered:STRG.1.1,muscle_stage2.filtered:STRG.1.1,ileum_stage3.filtered:STRG.1.1,lung_stage3.filtered:STRG.1.1,cerebellum_stage2.filtered:STRG.1.1,muscle_stage3.filtered:STRG.1.1,ref:G1.1"; contains_count "22"; 3p_dists_to_3p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; 5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; flrpm "58.872946"; longest "skin_stage2.filtered:STRG.1.1,cerebellum_stage3.filtered:STRG.1.1,kidney_stage3.filtered:STRG.1.1,cerebellum_stage1.filtered:STRG.1.1,skin_stage1.filtered:STRG.1.1,kidney_stage2.filtered:STRG.1.1,ileum_stage2.filtered:STRG.1.1,kidney_stage1.filtered:STRG.1.1,lung_stage1.filtered:STRG.1.1,skin_stage3.filtered:STRG.1.1,liver_stage3.filtered:STRG.1.1,liver_stage2.filtered:STRG.1.1,liver_stage1.filtered:STRG.1.1,muscle_stage1.filtered:STRG.1.1,ileum_stage1.filtered:STRG.1.1,lung_stage2.filtered:STRG.1.1,muscle_stage2.filtered:STRG.1.1,ileum_stage3.filtered:STRG.1.1,lung_stage3.filtered:STRG.1.1,cerebellum_stage2.filtered:STRG.1.1,muscle_stage3.filtered:STRG.1.1,ref:G1.1"; longest_FL_supporters "liver_stage1.filtered:STRG.1.1,liver_stage3.filtered:STRG.1.1,muscle_stage2.filtered:STRG.1.1,ileum_stage3.filtered:STRG.1.1,muscle_stage3.filtered:STRG.1.1,cerebellum_stage2.filtered:STRG.1.1,cerebellum_stage3.filtered:STRG.1.1,cerebellum_stage1.filtered:STRG.1.1,skin_stage3.filtered:STRG.1.1,lung_stage1.filtered:STRG.1.1,kidney_stage1.filtered:STRG.1.1,lung_stage2.filtered:STRG.1.1,muscle_stage1.filtered:STRG.1.1,ileum_stage1.filtered:STRG.1.1,liver_stage2.filtered:STRG.1.1,ref:G1.1,lung_stage3.filtered:STRG.1.1,skin_stage2.filtered:STRG.1.1,kidney_stage2.filtered:STRG.1.1,ileum_stage2.filtered:STRG.1.1,skin_stage1.filtered:STRG.1.1,kidney_stage3.filtered:STRG.1.1"; longest_FL_supporters_count "22"; mature_RNA_length "1801"; meta_3p_dists_to_5p "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1"; meta_5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; rpm "58.872946"; spliced "1"; ref_gene_id "G1"; tmerge_tr_id "TM_000000000001"; transcript_biotype ""; ensembl_gene_id "ENSSSCG00000048769"; best_ref "ENSSSCG00000048769"; feelnc_biotype "lncRNA";
# 29101 (14 fields)
# 17616 (16 fields)
# 92639 (18 fields)
# 850391 (20 fields)
# 32167 (22 fields)
# 1227171 (24 fields)
# 37912 (44 fields)
# 78154 (46 fields)
# 13411 (48 fields)
# 139157 (50 fields)



{
    # s is the string where we will put the new row to be written for each row
    s="";
    split($0,a,"\t");
    # first we wrtie the 8 first fields
    for(l=1; l<=8; l++)
    {
	s=(s)(a[l])("\t");
    }
    # then we deal with the more complicated 9th field in which we have to keep only the last feelnc_biotype value, 
    # but we dont know in advance how many there will be, so write it at the end
    split(a[9],b," ");
    k=1;
    while(b[k]!="")
    {
	if(b[k]!="feelnc_biotype")
	{
	    s=(s)(b[k])" "(b[k+1])" ";
	}
	else
	{
	    sfeelnc=(b[k])" "(b[k+1]);
	}
	k+=2;
    }
    if($3=="exon"||$3=="transcript")
    {
	print (s)(sfeelnc);
    }
    else
    {
	print substr(s,1,length(s)-1);
    }
}

