# srnaseq.tablist.to.dcc.ena.tab.awk
# !!! script not totally general since need to put the species name in the description manually !!!
# !!! on March 16th 2022 I put a much smaller description for the sake of conversion from tsv to excel !!!
# !!! the long one for pig was !!!
# RNA sequencing of pig tissues for small RNA annotation and expression analysis. Tissue specific RNA-seq data was generated to support annotation of small non-coding genes (in particular miRNAs) and to measure tissue specific expression. This study is part of the FAANG project, promoting rapid prepublication of data to support the research community. These data are released under Fort Lauderdale principles, as confirmed in the Toronto Statement (Toronto International Data Release Workshop. Birney et al. 2009. Pre-publication data sharing. Nature 461:168-170). Any use of this dataset must abide by the FAANG data sharing principles. Data producers reserve the right to make the first publication of a global analysis of this data. If you are unsure if you are allowed to publish on this dataset, please contact the FAANG Data Coordination Centre and FAANG Consortium (email: sarah.djebali@inserm.fr, sylvain.foissac@inrae.fr, faang-dcc@ebi.ac.uk, and cc faang@iastate.edu) to enquire. The full guidelines can be found at http://www.faang.org/data-share-principle.
# !!! and it now is !!!
# RNA sequencing


# very similar to srnaseq.bedlist.to.dcc.ena.tab.awk so would be good to have a single awk file for both
# main differences are, apart from file extension (bed vs tab):
###############################################################
# - id in 1st column, which is basename_predictions for bed files and just the basename of the tab file for tab file
# - where to get the metadata info from the basename (from 9th to 13th for bed file, 5 first fields for tab file)
# - type of file, SEQUENCE_ANNOTATION for bed files, SEQUENCE_FLATFILE for tab files

# this script takes as input:
#############################
# - a tsv file without header with at least 2 columns of which the first one is the absolute path to an srnaseq tab file
#   corresponding to edgeR cpm normalized expression of mature or hairpin mirnas as provided by the nf-core smrnaseq pipeline
#   in edgeR/miRBase_* subdirs and split into as many files as samples (expect two tab files per sample, one for matures and one
#   for hairpins), and the second one is the md5sum of this file
# - as the fileRef1 variable a read file that has the correspondence between run number (ERR) and the features of metadata
#   file included in the tab file names such as pig_stage1_fetusday30_cerebellum_rep1 (sp, stage, tissue, repnb)
# - as the fileRef2 variable an ena file that has information about study, sample, experiment, run and
#   metadata of the sample (very specific to our example, not really generalizable to anything else)
# this script makes as output:
##############################
# - a tsv file that will be provided as part of the ena tab in the dcc submission file for the srnaseq analysis
#   of a set of samples. This file has a header and then as many rows as tab files in the input file
# Notes:
########
# 1. we will take the basename of the tab file without the .tab extension as an alias

# To improve:
#############
# 1. make it possible to have the columns of the ena file in any number and use the header instead to know what is what
# 2. have a single awk file to deal with bed and tab files

# example
#########
# cd /work/project/fragencode/workspace/geneswitch/data/metadata/srnaseq
# pgm=/work/project/fragencode/tools/multi/Scripts/srnaseq.tablist.to.dcc.ena.tab.awk
# time awk -v fileRef1=srnaseq.read1.sus_scrofa.tsv -v fileRef2=ena_export_pig_srnas_PRJEB42001.tsv -f $pgm sus_scrofa.srnaseq.tablist.md5sum.tsv > sus_scrofa.srnaseq.tabfile.ena.tab.tsv 
# real	0m0.014s

# fileRef1=srnaseq.read1.sus_scrofa.tsv
# file	species	tissue	stage	replicate	sex	development	study	run	specimen	organism	labExpId
# SSR-INRA-GS-WP1-smallRNAseq-Liver-FT-30dpf-POOL-1_R*_001.fastq.gz	pig	liver	stage1	rep1	male	fetusday30	PRJEB42041	ERR4991913	SAMEA7629259		pig_stage1_fetusday30_liver_rep1
# 84 (11 fields)
# 1 (12 fields)  *** here corres between features of metadata present in bam file name and run id
#                *** wrong study id here so trust the ena file better for that

# fileRef2=ena_export_pig_srnas_PRJEB42001.tsv
# study_accession	sample_accession	experiment_accession	run_accession	tax_id	scientific_name	library_name	experiment_title	study_title	study_alias	experiment_alias	run_alias	sample_alias	sample_title
# PRJEB42001	SAMEA7629259	ERX4801438	ERR4991913	9823	Sus scrofa	pig_small_RNAseq_01	Illumina HiSeq 4000 paired end sequencing; Pig transcriptome and gene expression atlas	GENE-SWitCH Pig transcriptome and gene expression atlas (smallRNA-seq)	SSC_INRA_GS_WP1_RNAseq_smallRNAseq	SSC_INRA_GS_WP1_RNASEQ_Liver_FT_30dpf_POOL_1_smallRNAseq	SSC_INRA_GS_WP1_smallRNAseq_1	SAMEA7629259	Pig_large_white_Liver_30dpf
# 1 (14 fields)
# 84 (33 fields)   *** here we have a 1 to 1 corres between the run id and the exp id (1 run per exp here)
#                  *** the study id from this file is more reliable than in the read file (where need to be corrected for pig)

# sus_scrofa.srnaseq.tablist.md5sum.tsv
# /work/project/fragencode/workspace/geneswitch/results/srnaseq/sus_scrofa/nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28/edgeR/miRBase_mature/pig_stage1_fetusday30_cerebellum_rep1_1.mature.cpm.tab	b3c09276943df25e9dc415ff7c1ed958
# 168 (2 fields)

# output sus_scrofa.srnaseq.tabfile.ena.tab.tsv 
# Alias	Title	Analysis Type	Description	Study	Samples	Samples	Experiments	Experiments	Runs	Runs	Related Analyses	File Names	File Types	Checksum Methods	Checksums	Analysis Center	Analysis Date	Unit
# pig_stage1_fetusday30_cerebellum_rep1_1.mature.cpm	GENE-SWitCH Pig transcriptome and gene expression atlas (smallRNA-seq)	SEQUENCE_FLATFILE	RNA sequencingPRJEB42001	SAMEA7629264		ERX4801454		ERR4991929			pig_stage1_fetusday30_cerebellum_rep1_1.mature.cpm.tab	tab	MD5	b3c09276943df25e9dc415ff7c1ed958	INRAE Centre Toulouse Occitanie	2021-11-29	YYYY-MM-DD
# 169 (26 fields) *** 2nd samples, 2nd experiments, 2nd runs, related analyses, fine names are all empty here


BEGIN{
    OFS="\t";
    # read the read file to get the correspondence between the tab file name and the run number (reverse procedure compared to making symlinks)
    while (getline < fileRef1 >0)
    {
	runid[$2"_"$4"_"$7"_"$3"_"$5]=$9;
    }
    
    # read the ena file to get the correspondence between the run id and the study id, the sample id and the expt id
    while (getline < fileRef2 >0)
    {
	studyid[$4]=$1;
	sampleid[$4]=$2;
	expid[$4]=$3;
    }
    
    # print the header (19 columns)
    print "Alias", "Title", "Analysis Type", "Description", "Study", "Samples", "Samples", "Experiments", "Experiments", "Runs", "Runs", "Related Analyses", "File Names", "File Types", "Checksum Methods", "Checksums", "Analysis Center", "Analysis Date", "Unit";
}

{
    # for each tab file from the list, put the basename of the file without the .tab extension as an alias, and from its metadata get the runid and the other ids
    # note the basename of the file looks like this
    # pig_stage1_fetusday30_cerebellum_rep1_1.mature.cpm.tab
    n=split($1,a,"/");
    split(a[n],b,".tab");
    split(b[1],c,"_");
    run=runid[c[1]"_"c[2]"_"c[3]"_"c[4]"_"c[5]];
    # write the 19 columnms (some of which are empty)
    print b[1], "GENE-SWitCH Chicken transcriptome and gene expression atlas (smallRNA-seq)", "SEQUENCE_ANNOTATION", "RNA sequencing", studyid[run], sampleid[run], "", expid[run], "", run, "", "", "srnaseq/tabfiles/"a[n], "tab", "MD5", $2, "INRAE Centre Toulouse Occitanie", "2021-11-29", "YYYY-MM-DD";
}
