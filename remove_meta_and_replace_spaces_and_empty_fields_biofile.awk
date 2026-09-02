# remove_meta_and_replace_spaces_and_empty_fields_biofile.awk

# example
# basedir=/work/project/bridge/workspace/l3_loanemadrid_lms
# cd $basedir/results
# pgm=$basedir/scripts/remove_meta_and_replace_spaces_and_empty_fields_biofile.awk
# awk -f $pgm ../data/Annot_bio_LMS_ICGC.txt > Annot_bio_LMS_ICGC.wometa.ok.tsv

# input file = ../data/Annot_bio_LMS_ICGC.txt
# Sample	Proportion_RB1_WT	RB1_mut_typ_allele1	CIRCOS RB1	Proportion_TP53_WT	TP53_mut_typ_allele1	CIRCOS p53	ATRX_mut_typ_allele1	Proportion_PTEN_WT	PTEN CIRCOS	CINSARC_LOOCV	Gaelle_Ploidie	ssGSEA_NES_CINSARC	TMM	iTRAC_DNA	iRACIN_DNA	LMS_group
# LMS1	0	SV	chromo	0	CNA	chromo	SV	0	del hetero	C2	Tetra	3,945060208	ALT	High	Low	Other
# 3 (15 fields)
# 70 (16 fields)
# 5 (17 fields)
# 23 (18 fields)
# 6 (19 fields)
# 2 (20 fields)
# 1 (22 fields)   *** 110 rows  *** here and unlike in the clinical file where we had the primary tumor and the recidive, LMS85 is present only once

# output file = Annot_bio_LMS_ICGC.wometa.ok.tsv
# Sample	Proportion_RB1_WT	RB1_mut_typ_allele1	CIRCOS_RB1	Proportion_TP53_WT	TP53_mut_typ_allele1	CIRCOS_p53	ATRX_mut_typ_allele1	Proportion_PTEN_WT	PTEN_CIRCOS	CINSARC_LOOCV	Gaelle_Ploidie	ssGSEA_NES_CINSARC	TMM	iTRAC_DNA	iRACIN_DNA	LMS_group
# LMS1	0	SV	chromo	0	CNA	chromo	SV	0	del_hetero	C2	Tetra	3,945060208	ALT	High	Low	Other
# 108 (17 fields)  *** 2 rows removed because only LMS129 and LMS130 were there


BEGIN{
    OFS="\t"; 
    ko["LMS26"]=1; 
    ko["LMS75"]=1; 
    ko["LMS82"]=1; 
    ko["LMS129"]=1; 
    ko["LMS130"]=1;
} 

{
    if(NR==1||ko[$1]!=1)
    {
        gsub(/\ /,"_",$0); 
        n=split($0,a,"\t"); 
        s=""; 
        for(k=1; k<n; k++)
        {
            s=(s)(a[k]=="" ? "NA" : a[k])("\t");
        } 
        print (s)(a[k]=="" ? "NA" : a[k]);
    }
}   