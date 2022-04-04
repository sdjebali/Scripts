# srnaseq.bamlist.to.dcc.faang.tab.awk
# !!! on March 31st 2022 add two new necesary colums !!! 
# !!! nextflow config url and	nextflow spreadsheet url !!! 
# !!!! but with empty values there !!!


# this script takes as input:
#############################
# - a tsv file without header that has at least 1 column with a list of absolute paths to srnaseq bam files corresponding to alignments
#   of srnaseq reads to either mirna matures or mirna hairpins as output by the nf-core smrnaseq pipeline
#   in bowtie subdirs and expected to end in .sorted.bam (expect twice as many bam files as the number of samples) 
# - a genome variable that is the name of the genome assembly used (eg Sscrofa11.1)
# this script makes as output:
##############################
# - a tsv file that will be provided as the faang tab in the dcc submission file for the srnaseq analysis
#   of a set of samples. This file has a header and then as many rows as bam files in the input file
# Notes:
########
# 1. many things on those rows are hard-coded here such as the bioinformatics protocol and others
#    since the header and the first row should look like this
# Alias	Project	Secondary Project	Assay Type	Analysis Protocol	Analysis Code	Analysis Code Version	Reference Genome
# labExpId_processed_reads	FAANG	GENE-SWitCH	microRNA profiling by high throughput sequencing	https://data.faang.org/api/fire_api/analyses/INSERM-INRAE_SOP_srnaseq-processing_20211129.pdf	https://github.com/nf-core/smrnaseq	v1.1.0	Sscrofa11.1
# 2. we will take the basename of the bam file without .bam as an alias 

# example
#########
# cd ~/fragencode/workspace/geneswitch/data/metadata/srnaseq
# pgm=/work/project/fragencode/tools/multi/Scripts/srnaseq.bamlist.to.dcc.faang.tab.awk
# time awk -v genome=Sscrofa11.1 -f $pgm sus_scrofa.srnaseq.bamlist.md5sum.tsv > sus_scrofa.srnaseq.bamfile.faang.tab.tsv 
# real	0m0.758s


# sus_scrofa.srnaseq.bamlist.md5sum.tsv
# /work2/project/fragencode/workspace/geneswitch/results/srnaseq/sus_scrofa/nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28/bowtie/miRBase_mature/pig_stage1_fetusday30_cerebellum_rep1_1.mature.sorted.bam	2a1545e0a6e3950ff8cd28b50338aa31
# /work2/project/fragencode/workspace/geneswitch/results/srnaseq/sus_scrofa/nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28/bowtie/miRBase_mature/pig_stage1_fetusday30_cerebellum_rep2_1.mature.sorted.bam	118ea871a68f1e47a9116de177ee2bf6
# 168 (2 fields)

# output sus_scrofa.srnaseq.bamfile.faang.tab.tsv 
# Alias	Project	Secondary Project	Assay Type	Analysis Protocol	Analysis Code	Analysis Code Version	Reference Genome
# pig_stage1_fetusday30_cerebellum_rep1_1.mature.sorted	FAANG	GENE-SWitCH	microRNA profiling by high throughput sequencing	https://data.faang.org/api/fire_api/analyses/INSERM-INRAE_SOP_srnaseq-processing_20211129.pdf	https://github.com/nf-core/smrnaseq	v1.1.0	Sscrofa11.1
# 168 (13 fields)
# 1 (14 fields)   *** plus the two last columns in theory 


BEGIN{
    OFS="\t";
    # print the header
    print "Alias", "Project", "Secondary Project", "Assay Type", "Analysis Protocol", "Analysis Code", "Analysis Code Version", "Reference Genome", "nextflow config url", "nextflow spreadsheet url";
}

{
    n=split($1,a,"/");
    split(a[n],b,".bam");
    print b[1], "FAANG", "GENE-SWitCH", "microRNA profiling by high throughput sequencing", "https://data.faang.org/api/fire_api/analyses/INSERM-INRAE_SOP_srnaseq-processing_20211129.pdf", "https://github.com/nf-core/smrnaseq", "v1.1.0", genome, "", "";
}
