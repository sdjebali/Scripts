# srnaseq.bamlist.to.dcc.ena.tab.awk


# this script takes as input:
#############################
# - a tsv file without header with at least 2 columns of which the first one is the absolute path to an srnaseq bam files
#   corresponding to alignments of srnaseq reads to either mirna matures or mirna hairpins output by the nf-core smrnaseq
#   pipeline in bowtie subdirs (expect twice as many bam files as the number of samples) and the second one is the md5sum
#   of this file
# - as the fileRef1 variable a read file that has the correspondence between run number (ERR) and the features of metadata
#   file included in the bam file names such as pig_stage1_fetusday30_cerebellum_rep1 (sp, stage, tissue, repnb)
# - as the fileRef2 variable an ena file that has information about study, sample, experiment, run and
#   metadata of the sample (very specific to our example, not really generalizable to anything else)
# this script makes as output:
##############################
# - a tsv file that will be provided as the ena tab in the dcc submission file for the srnaseq analysis
#   of a set of samples. This file has a header and then as many rows as bam files in the input file
# Notes:
########
# 1. we will take the basename of the bam file without .sorted.bam as an alias

# To improve:
# make it possible to have the columns of the ena file in any number and use the header instead to know what is what

# example
#########
# cd /work2/project/fragencode/workspace/geneswitch/data/metadata/srnase
# pgm=/work2/project/fragencode/tools/multi/Scripts/srnaseq.bamlist.to.dcc.ena.tab.awk
# time awk -v fileRef1=srnaseq.read1.sus_scrofa.tsv -v fileRef2=ena_export_pig_srnas_PRJEB42001.tsv -f $pgm sus_scrofa.srnaseq.bamlist.md5sum.tsv > sus_scrofa.srnaseq.pig.bamfile.ena.tab.tsv 


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

# sus_scrofa.srnaseq.bamlist.md5sum.tsv
# /work2/project/fragencode/workspace/geneswitch/results/srnaseq/sus_scrofa/nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28/bowtie/miRBase_mature/pig_stage1_fetusday30_cerebellum_rep1_1.mature.sorted.bam	2a1545e0a6e3950ff8cd28b50338aa31
# 168 (2 fields)


# output sus_scrofa.srnaseq.pig.bamfile.ena.tab.tsv 
# Alias	Title	Analysis Type	Description	Study	Samples	Samples	Experiments	Experiments	Runs	Runs	Related Analyses	File Names	File Types	Checksum Methods	Checksums	Analysis Center	Analysis Date	Unit
# pig_stage1_fetusday30_cerebellum_rep1_1.mature	GENE-SWitCH Pig transcriptome and gene expression atlas (smallRNA-seq)	PROCESSED_READS	RNA sequencing of pig tissues for small RNA annotation and expression analysis. Tissue specific RNA-seq data was generated to support annotation of small non-coding genes (in particular miRNAs) and to measure tissue specific expression. This study is part of the FAANG project, promoting rapid prepublication of data to support the research community. These data are released under Fort Lauderdale principles, as confirmed in the Toronto Statement (Toronto International Data Release Workshop. Birney et al. 2009. Pre-publication data sharing. Nature 461:168-170). Any use of this dataset must abide by the FAANG data sharing principles. Data producers reserve the right to make the first publication of a global analysis of this data. If you are unsure if you are allowed to publish on this dataset, please contact the FAANG Data Coordination Centre and FAANG Consortium (email: sarah.djebali@inserm.fr, sylvain.foissac@inrae.fr, faang-dcc@ebi.ac.uk, and cc faang@iastate.edu) to enquire. The full guidelines can be found at http://www.faang.org/data-share-principle.										bam	MD5	2a1545e0a6e3950ff8cd28b50338aa31	INRAE Centre Toulouse Occitanie	2021-11-29	YYYY-MM-DD
# 1 (26 fields)
# 168 (169 fields)  *** 2nd samples, 2nd experiments, 2nd runs, related analyses, fine names are all empty here


BEGIN{
    OFS="\t";
    # read the reads file to get the correspondence between the bam file name and the run number (reverse procedure compared to making symlinks)
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
    # for each bam file from the list, put the basename of the file without .sorted.bam as an alias, and from its metadata get the runid and the other ids
    n=split($1,a,"/");
    split(a[n],b,".sorted.bam");
    split(b[1],c,".mature");
    run=runid[c[1]];
    # write the 19 columnms (some of which are empty)
    print b[1], "GENE-SWitCH Pig transcriptome and gene expression atlas (smallRNA-seq)", "PROCESSED_READS", "RNA sequencing of pig tissues for small RNA annotation and expression analysis. Tissue specific RNA-seq data was generated to support annotation of small non-coding genes (in particular miRNAs) and to measure tissue specific expression. This study is part of the FAANG project, promoting rapid prepublication of data to support the research community. These data are released under Fort Lauderdale principles, as confirmed in the Toronto Statement (Toronto International Data Release Workshop. Birney et al. 2009. Pre-publication data sharing. Nature 461:168-170). Any use of this dataset must abide by the FAANG data sharing principles. Data producers reserve the right to make the first publication of a global analysis of this data. If you are unsure if you are allowed to publish on this dataset, please contact the FAANG Data Coordination Centre and FAANG Consortium (email: sarah.djebali@inserm.fr, sylvain.foissac@inrae.fr, faang-dcc@ebi.ac.uk, and cc faang@iastate.edu) to enquire. The full guidelines can be found at http://www.faang.org/data-share-principle.", studid[run], sampleid[run], "", expid[run], "", run, "", "", "", "bam", "MD5", $2, "INRAE Centre Toulouse Occitanie", "2021-11-29", "YYYY-MM-DD";
}
