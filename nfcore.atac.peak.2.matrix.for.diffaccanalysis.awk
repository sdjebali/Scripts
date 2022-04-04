# nfcore.atac.peak.2.matrix.for.diffaccanalysis.awk
# (inspired from the script atac.qc.techrep.2.biorep.awk)
# this script takes as input:
#############################
# - as main input a raw read count peak matrix from nf-core atacseq pipeline with bioreplicates as
#   columns and some additional information about the peaks before that
# - fileRef1 which is the csv file taken as input by the pipeline that has the correspondance between
#   fastq files of each technical replicate, their belonging to a given biorep and to a given combination
#   of tissue/stage/diet
# - fileRef2 which is a tsv file with header with metadata information connecting the fastq file name
#   with a lot of information such as tissue, stage, diet, animal id, animal name, gendre, mother id, mother breed and father id
# this script produces as output:
#################################
# - as main output a raw read count peak matrix which as almost the same information as the main input file except that
#   * it does not have a command line as first row but just a header
#   * it has n+1 columns where n is the number of bioreplicates. This means the only information about the peak
#     structure is its coordinates as chr:gbeg:gend:+ (eg 1:21547:24249:+)
#   * the biorep ids in the header are of the form animalname:tissue
# TODO: have a file with corres between biorep1 and biorep2 beforehand which would avoid the two first fileRef files here (and adding a simple tsv one)
# but also the big code reading them in the begin, but need to be tested (could also help in the script nfcore.atac.peak.2.matrix.for.diffaccanalysis.awk)

# Example of usage:
###################
# srun --mem=8G --pty bash
# cd /work/project/fragencode/workspace/sdjebali/geneswitch/wp5/atacseq/qc.measures
# pgm=~/fragencode/tools/multi/Scripts/nfcore.atac.peak.2.matrix.for.diffaccanalysis.awk
# indir=/work/project/geneswitch/results/atacseq/sus_scrofa/wp5/complete/nf-core.atacseq.v1.2.1.Sscrofa11.1.102.22-03-07
# meta=~/geneswitch/data/metadata/atacseq/wp5/fastq_tomotherinfo.ok.tsv
# peaks=$indir/bwa/mergedReplicate/macs/broadPeak/consensus/consensus_peaks.mRp.clN.featureCounts.txt
# time awk -v fileRef1=$indir/samplesheet.csv -v fileRef2=$meta -f $pgm $peaks > peakid.rawcount.eachbiorep.tsv
# real	0m24.835s


# fileRef1 input file
#####################
# $indir/samplesheet.csv
# group,replicate,fastq_1,fastq_2
# liver_fetus_control,1,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-AHJK7WDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-AHJK7WDSX2_L004_R2.fastq.gz
# liver_fetus_control,1,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-BHMYNVDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-BHMYNVDSX2_L004_R2.fastq.gz
# liver_fetus_control,2,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-2_GGCTTAAG-GGTCACGA-AHJK7WDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-2_GGCTTAAG-GGTCACGA-AHJK7WDSX2_L004_R2.fastq.gz
# liver_fetus_control,2,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-2_GGCTTAAG-GGTCACGA-BHMYNVDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-2_GGCTTAAG-GGTCACGA-BHMYNVDSX2_L004_R2.fastq.gz
# 641 (1 fields)

# fileRef2 input file
#####################
# $meta
# fastqfile	readtype	runnumber	tissue	age	motherid	diet	sampleid	samplename	sex	motherbreed	fatherid
# Skel-mus-Fetus-109-1_TTGGACTC-CTGCTTCC-BHMYNVDSX2_L004_R2.fastq.gz	R2	BHMYNVDSX2	skeletalmuscle	fetus	109	control	109.1	SSC_WU_WP5_FT_70dpf_109.1	male	TN70	E6211
# 1281 (12 fields)  *** $9":"$4":"$3 identify a techrep uniquely, while $9":"$4 identify a biorep uniquely

# main input file
##################
# $peaks
# # Program:featureCounts v2.0.1; Command:"featureCounts" "-F" "SAF" "-O" "--fracOverlap" "0.2" "-T" "12" "-p" "--donotsort" "-a" "consensus_peaks.mRp.clN.saf" "-o" "consensus_peaks.mRp.clN.featureCounts.txt" "skelmus_fetus_control_R28.mLb.clN.bam" "skelmus_fetus_lowfibre_R9.mLb.clN.bam" ... "skelmus_fetus_lowfibre_R7.mLb.clN.bam" "skelmus_piglet_lowfibre_R10.mLb.clN.bam" 
# Geneid	Chr	Start	End	Strand	Length	skelmus_fetus_control_R28.mLb.clN.bam	skelmus_fetus_lowfibre_R9.mLb.clN.bam	...	skelmus_fetus_lowfibre_R7.mLb.clN.bam	skelmus_piglet_lowfibre_R10.mLb.clN.bam
# Interval_1	1	21547	24249	+	2703	656	392	...	563	330
# 245472 (326 fields) *** almost same number as before with tech rep, which is good
# 1 (337 fields)


# main output
#############
# peakid.rawcount.eachbiorep.tsv
# peakid	SSC_WU_WP5_FT_70dpf_913.4:skeletalmuscle	SSC_WU_WP5_FT_70dpf_117.1:skeletalmuscle	...	SSC_WU_WP5_FT_70dpf_115.3:skeletalmuscle	SSC_WU_WP5_Piglet_5768:skeletalmuscle
# 1:21547:24249:+	656	392	...	563	330
# 245472 (321 fields) 


BEGIN{
    OFS="\t";

    # read the design file of the pipeline remembering for each biorep no, the first of the 4 fastq files corresponding
    # to this biorep to then be able to assoicate this biorep to the samplename:tissue biorep id we want in the main output file
    # group,replicate,fastq_1,fastq_2
    # liver_fetus_control,1,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-AHJK7WDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-AHJK7WDSX2_L004_R2.fastq.gz
    # liver_fetus_control,1,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-BHMYNVDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-BHMYNVDSX2_L004_R2.fastq.gz
    while (getline < fileRef1 >0)
    {
	n++;
	if(n>=2)
	{
	    split($1,a,",");
	    currbiorepid=a[1]"_R"a[2];
	    seen[currbiorepid]++;
	    if(seen[currbiorepid]==1)
	    {
		n1=split(a[3],b,"/");
		fastqfile[currbiorepid]=b[n1];
	    }
	}
    }

    # read the metadata file of the reads to associate each read1 file to its samplename:tissue
    # $9":"$4 (like SSC_WU_WP5_FT_70dpf_109.1:skeletalmuscle) is the unique identifier of a biorep that we want in the main output file header
    # fastqfile	readtype	runnumber	tissue	age	motherid	diet	sampleid	samplename	sex	motherbreed	fatherid
    # Skel-mus-Fetus-109-1_TTGGACTC-CTGCTTCC-BHMYNVDSX2_L004_R2.fastq.gz	R2	BHMYNVDSX2	skeletalmuscle	fetus	109	control	109.1	SSC_WU_WP5_FT_70dpf_109.1	male	TN70	E6211
    # 1281 (12 fields)
    while (getline < fileRef2 >0)
    {
	if($2=="R1")
	{
	    biorepid[$1]=$9":"$4
	}
    }
}


# Takes the 2nd row as the header and edit it to only have peakid and then the id of each biorep as samplename:tissue (where samplename is in fact animalname)
# Geneid	Chr	Start	End	Strand	Length	skelmus_fetus_control_R28.mLb.clN.bam	skelmus_fetus_lowfibre_R9.mLb.clN.bam	...	skelmus_fetus_lowfibre_R7.mLb.clN.bam	skelmus_piglet_lowfibre_R10.mLb.clN.bam
NR==2{
    s="peakid\t";
    for(i=7; i<=(NF-1); i++)
    {
	split($i,a,".mLb.clN.bam");
	s=(s)(biorepid[fastqfile[a[1]]])("\t");
    }
    split($i,a,".mLb.clN.bam");
    print (s)(biorepid[fastqfile[a[1]]]);
}

# From the 3rd row, write the peak coordinates and the raw count in each biorep so just remove columns 2 to 6 and make column 1 out of the coord in columns 2 to 5
# Interval_1	1	21547	24249	+	2703	656	392	...	563	330
NR>=3{
    s=$2":"$3":"$4":"$5"\t";
    for(i=7; i<=(NF-1); i++)
    {
	s=(s)($i)("\t")
    }
    print (s)($i);
}

 
