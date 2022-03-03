# srnaseq.bedlist.to.dcc.faang.tab.awk
# very similar to srnaseq.tablist.to.dcc.faang.tab.awk, see improvement note below

# this script takes as input:
#############################
# - a tsv file without header that has at least 1 column with a list of absolute paths to srnaseq bed files corresponding to mirdeep2
#   mirna predictions output by the nf-core smrnaseq pipeline in mirdeep2/mirdeep subdir and further filtered on their prob to be a TP
#   (expect one per sample)
# - a genome variable that is the name of the genome assembly used (eg Sscrofa11.1)
# this script makes as output:
##############################
# - a tsv file that will be provided as part of the faang tab in the dcc submission file for the srnaseq analysis
#   of a set of samples. This file has a header and then as many rows as bed files in the input file
# Notes:
########
# 1. many things on those rows are hard-coded here such as the bioinformatics protocol and others
#    since the header and the first row should look like this
# Alias	Project	Secondary Project	Assay Type	Analysis Protocol	Analysis Code	Analysis Version	Reference Genome
# labExpId_processed_reads	FAANG	GENE-SWitCH	microRNA profiling by high throughput sequencing	https://data.faang.org/api/fire_api/analyses/INSERM-INRAE_SOP_srnaseq-processing_20211129.pdf	https://github.com/nf-core/smrnaseq	v1.1.0	Sscrofa11.1
# 2. we will take the basename of the bed file without the .bed extension as an alias, although in the case if mirdeep2 files
#    there is a strange begining of file ...
# Improvements
##############
# in the future we could make a single script for this one and srnaseq.bedlist.to.dcc.faang.tab.awk by taking out as a parameter the extension
# of the file (bed or tab) and by adding _predictions for the id of the 1st column in one case (bed files) and nothing in the other (tab files)

# example
#########
# cd ~/fragencode/workspace/geneswitch/data/metadata/srnaseq
# pgm=/work/project/fragencode/tools/multi/Scripts/srnaseq.bedlist.to.dcc.faang.tab.awk
# time awk -v genome=Sscrofa11.1 -f $pgm sus_scrofa.srnaseq.bedlist.md5sum.tsv > sus_scrofa.srnaseq.bedfile.faang.tab.tsv 
# real	0m0.005s

# sus_scrofa.srnaseq.bedlist.md5sum.tsv
# /work/project/fragencode/workspace/geneswitch/results/srnaseq/sus_scrofa/nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28/mirdeep2/mirdeep/result_28_06_2021_t_16_37_24_pig_stage3_newbornday1_liver_rep1_1.bed	2a7bcbbc3043e1863bec32898641cd44
# 84 (2 fields)

# output sus_scrofa.srnaseq.bedfile.faang.tab.tsv 
# Alias	Project	Secondary Project	Assay Type	Analysis Protocol	Analysis Code	Analysis Version	Reference Genome
# result_28_06_2021_t_16_37_24_pig_stage3_newbornday1_liver_rep1_1_predictions	FAANG	GENE-SWitCH	microRNA profiling by high throughput sequencing	https://data.faang.org/api/fire_api/analyses/INSERM-INRAE_SOP_srnaseq-processing_20211129.pdf	https://github.com/nf-core/smrnaseq	v1.1.0	Sscrofa11.1
# 84 (13 fields)
# 1 (14 fields)


BEGIN{
    OFS="\t";
    # print the header
    print "Alias", "Project", "Secondary Project", "Assay Type", "Analysis Protocol", "Analysis Code", "Analysis Version", "Reference Genome";
}

{
    n=split($1,a,"/");
    split(a[n],b,".bed");
    print b[1]"_predictions", "FAANG", "GENE-SWitCH", "microRNA profiling by high throughput sequencing", "https://data.faang.org/api/fire_api/analyses/INSERM-INRAE_SOP_srnaseq-processing_20211129.pdf", "https://github.com/nf-core/smrnaseq", "v1.1.0", genome;
}
