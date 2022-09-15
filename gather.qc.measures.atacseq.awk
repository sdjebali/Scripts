# gather.qc.measures.atacseq.awk
# this script gathers several kinds of qc metrix distributed in 5 different tsv files and derived from the output of
# the nf core atacseq pipeline run on data with 2 technical replicates for each bioreplicate:
# - a main input tsv file with header that has N rows representing N biological replicates and a header, with this info (9 columns)
#   * biorep
#   * clnread.tot
#   * clnreadmap.nb
#   * clnreadmap.pcent
#   * clnreadmt.nb
#   * clnreadmt.pcent
#   * frac.dupl.picard
#   * peak.nb
#   * frip
# - 6 fileRef files
#   * fileRef1 which is the csv file taken as input by the nf core atacseq pipeline that has the fastq files of each technical replicate, 
#     their belonging to a given biorep and to a given combination tissue/stage/diet
#   * fileRef2 which is a tsv file with header with metadata information connecting the fastq file name with a lot of information such as
#     tissue, stage, diet, animal id, animal name, gendre, mother id, mother breed and father id
#   * fileRef3=lid.nbinitreads.tsv that has for each technical replicate and read number the initial number of reads (4*N rows)
#   * fileRef4=lid.nbbpinit.aftertrim.nb.pcent.tsv that has for each techrep and read no the init nb of bp, the one after trimming and the % it represents (4*N rows)
#   * fileRef5=lid.gccontent.aftertrimming.summary.tsv that has for each techrep and read no whether it passes the fastqc GC content test (4*N rows)
#   * fileRef6=lid.dupllevel.aftertrimming.summary.tsv that has for each techrep and read no whether it passes the fastqc duplication level test (4*N rows)
# It will output a tsv file with header that has 34 fields in this order
########################################################################
#   * biorep id from the pipeline (such as liver_fetus_control_R10)
#   * biorep id of the raw count matrix (such as SSC_WU_WP5_FT_70dpf_913.4:skeletalmuscle)
#   * techrep1.read1.initnb
#   * techrep1.read2.initnb
#   * techrep2.read1.initnb
#   * techrep2.read2.initnb
#   * techrep1.read1.initbp.nb
#   * techrep1.read1.trimbp.nb
#   * techrep1.read1.trimbp.pcent
#   * techrep1.read2.initbp.nb
#   * techrep1.read2.trimbp.nb
#   * techrep1.read2.trimbp.pcent
#   * techrep2.read1.initbp.nb
#   * techrep2.read1.trimbp.nb
#   * techrep2.read1.trimbp.pcent
#   * techrep2.read2.initbp.nb
#   * techrep2.read2.trimbp.nb
#   * techrep2.read2.trimbp.pcent
#   * techrep1.read1.gcsummary
#   * techrep1.read2.gcsummary
#   * techrep2.read1.gcsummary
#   * techrep2.read2.gcsummary
#   * techrep1.read1.duplevsummary
#   * techrep1.read2.duplevsummary
#   * techrep2.read1.duplevsummary
#   * techrep2.read2.duplevsummary
#   * the same 9 first rows as the input file
# in a second step we can add some values, average others or select the worst in pass, warn and fail
# but here we just report the info about each tech rep of a biorep
# TODO: have a file with corres between biorep1 and biorep2 beforehand which would avoid the two first fileRef files here (and adding a simple tsv one)
# but also the big code reading them in the begin, but need to be tested (could also help in the script nfcore.atac.peak.2.matrix.for.diffaccanalysis.awk)
# TODO: allow a biorep to have more or less than 2 techrep

# example
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp5/atacseq/qc.measures
# pgm=~/fragencode/tools/multi/Scripts/gather.qc.measures.atacseq.awk
# indir=/work/project/geneswitch/results/atacseq/sus_scrofa/wp5/complete/nf-core.atacseq.v1.2.1.Sscrofa11.1.102.22-03-07
# meta=~/geneswitch/data/metadata/atacseq/wp5/fastq_tomotherinfo.ok.tsv
# time paste lid.totreads.mapped.nb.pcent.tsv lid.totmappedreads.mt.nb.pcent.tsv lid.fracduplicate.tsv lid.peaks.nb.frip.tsv | awk 'BEGIN{OFS="\t"; print "biorep", "clnread.tot", "clnreadmap.nb", "clnreadmap.pcent", "clnreadmt.nb", "clnreadmt.pcent", "frac.dupl.picard", "peak.nb", "frip"} {split($9,a,"_T1"); if($1==$5&&$1==a[1]&&$1==$11){print $1, $2, $3, $4, $7, $8, $10, $12, $13}}' | awk -v fileRef1=$indir/samplesheet.csv -v fileRef2=$meta -v fileRef3=lid.nbinitreads.tsv -v fileRef4=lid.nbbpinit.aftertrim.nb.pcent.tsv -v fileRef5=lid.gccontent.aftertrimming.summary.tsv -v fileRef6=lid.dupllevel.aftertrimming.summary.tsv -f $pgm > lid1.2.techrep1.2.read1.2.nb.bp.trim.nb.pcent.gc.duplev.techrep1.2.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.peaknb.frip.tsv
# real	0m0.060s

# fileRef1=$indir/samplesheet.csv
# group,replicate,fastq_1,fastq_2
# liver_fetus_control,1,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-AHJK7WDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-AHJK7WDSX2_L004_R2.fastq.gz
# liver_fetus_control,1,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-BHMYNVDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-1_CCAAGTCT-TCATCCTT-BHMYNVDSX2_L004_R2.fastq.gz
# liver_fetus_control,2,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-2_GGCTTAAG-GGTCACGA-AHJK7WDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-2_GGCTTAAG-GGTCACGA-AHJK7WDSX2_L004_R2.fastq.gz
# liver_fetus_control,2,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-2_GGCTTAAG-GGTCACGA-BHMYNVDSX2_L004_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-Fetus-109-2_GGCTTAAG-GGTCACGA-BHMYNVDSX2_L004_R2.fastq.gz
# 641 (1 fields)

# fileRef2=$meta
# fastqfile	readtype	runnumber	tissue	age	motherid	diet	sampleid	samplename	sex	motherbreed	fatherid
# Skel-mus-Fetus-109-1_TTGGACTC-CTGCTTCC-BHMYNVDSX2_L004_R2.fastq.gz	R2	BHMYNVDSX2	skeletalmuscle	fetus	109	control	109.1	SSC_WU_WP5_FT_70dpf_109.1	male	TN70	E6211
# 1281 (12 fields)

# fileRef3=lid.nbinitreads.tsv
# liver_fetus_control_R10_T1_1	68268904
# 1280 (2 fields)  *** always T1 and then T2 according to tech rep no and after that either 1 or 2 according to read no

# fileRef4=lid.nbbpinit.aftertrim.nb.pcent.tsv
# liver_fetus_control_R10_T1_1	10240335600	6449723020	62.9835
# 1280 (4 fields)  *** always T1 and then T2 according to tech rep no and after that either 1 or 2 according to read no

# fileRef5=lid.gccontent.aftertrimming.summary.tsv
# liver_fetus_control_R10_T1_1_val_1	PASS
# 1280 (2 fields)  *** only _T1_1_val_1 or _T1_2_val_2 or _T2_1_val_1 or _T2_2_val_2 according to tech rep no and then read no

# fileRef6=lid.dupllevel.aftertrimming.summary.tsv
# liver_fetus_control_R10_T1_1_val_1	PASS
# 1280 (2 fields)  *** only _T1_1_val_1 or _T1_2_val_2 or _T2_1_val_1 or _T2_2_val_2 according to tech rep no and then read no

# main input file
# biorep	clnread.tot	clnreadmap.nb	clnreadmap.pcent	clnreadmt.nb	clnreadmt.pcent	frac.dupl.picard	peak.nb	frip
# liver_fetus_control_R10	172461956	170078000	98.6177	1476758	0.85628	0.221747	76791	0.437052
# 321 (9 fields)

# main output file
# biorepid1	biorepid2	techrep1.read1.initnb	techrep1.read2.initnb	techrep2.read1.initnb	techrep2.read2.initnb	techrep1.read1.initbp.nb	techrep1.read1.trimbp.nb	techrep1.read1.trimbp.pcent	techrep1.read2.initbp.nb	techrep1.read2.trimbp.nb	techrep1.read2.trimbp.pcent	techrep2.read1.initbp.nb	techrep2.read1.trimbp.nb	techrep2.read1.trimbp.pcent	techrep2.read2.initbp.nb	techrep2.read2.trimbp.nb	techrep2.read2.trimbp.pcent	techrep1.read1.gcsummary	techrep1.read2.gcsummary	techrep2.read1.gcsummary	techrep2.read2.gcsummary	techrep1.read1.duplevsummary	techrep1.read2.duplevsummary	techrep2.read1.duplevsummary	techrep2.read2.duplevsummary	clnread.tot	clnreadmap.nb	clnreadmap.pcent	clnreadmt.nb	clnreadmt.pcent	frac.dupl.picard	peak.nb	frip
# liver_fetus_control_R10	SSC_WU_WP5_FT_70dpf_175.2:liver	68268904	68268904	65985728	65985728	10240335600	6449723020	62.9835	10240335600	6475616541	63.2364	9897859200	6266643954	63.3131	9897859200	6294061983	63.5901	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	172461956	170078000	98.6177	1476758	0.85628	0.221747	76791	0.437052
# 321 (34 fields)

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
    
    # reads the file of initial read number in each tech rep and read no
    while (getline < fileRef3 >0)
    {
	s="";
	s1="";
	split($1,a,"_");
	k=1;
	
	# until it hits T1 or T2 (indication of tech rep no) it continues to add info to s string corresponding to biorep id
	while(a[k]!="T1"&&a[k]!="T2")
	{
	    s=(s)(a[k])("_")
	    k++;
	}
	# if in techrep1 it writes to readinit1 (removing last _ to s into s1)
	if(a[k]=="T1")
	{
	    s1=substr(s,1,(length(s)-1));
	    readinit1[s1,a[k+1]]=$2;
	}
	else
	{
	    # if in techrep1 it writes to readinit2 (removing last _ to s into s1)
	    if(a[k]=="T2")
	    {
		s1=substr(s,1,(length(s)-1));
		readinit2[s1,a[k+1]]=$2;
	    }
	}
    }

    # reads the file of init bp, bp after trimming and percent of this
    while (getline < fileRef4 >0)
    {
	s="";
	s1="";
	split($1,a,"_");
	k=1;

	# until it hits T1 or T2 (indication of tech rep no) it continues to add info to s string corresponding to biorep id
	while(a[k]!="T1"&&a[k]!="T2")
	{
	    s=(s)(a[k])("_")
	    k++;
	}
	# if in techrep1 it writes to bpinit1, bptrim1 and pcenttrim1 (removing last _ to s into s1)
	if(a[k]=="T1")
	{
	    s1=substr(s,1,(length(s)-1));
	    bpinit1[s1,a[k+1]]=$2;
	    bptrim1[s1,a[k+1]]=$3;
	    pcenttrim1[s1,a[k+1]]=$4;
	}
	else
	{
	    # if in techrep1 it writes to bpinit2, bptrim2 and pcenttrim2 (removing last _ to s into s1)
	    if(a[k]=="T2")
	    {
		s1=substr(s,1,(length(s)-1));
		bpinit2[s1,a[k+1]]=$2;
		bptrim2[s1,a[k+1]]=$3;
		pcenttrim2[s1,a[k+1]]=$4;
	    }
	}
    }

    # reads the gc summary for each tech rep and read no
    # only _T1_1_val_1 or _T1_2_val_2 or _T2_1_val_1 or _T2_2_val_2 according to tech rep no and then read no
    while (getline < fileRef5 >0)
    {
	s="";
	s1="";
	split($1,a,"_");
	k=1;

	# until it hits T1 or T2 (indication of tech rep no) it continues to add info to s string corresponding to biorep id
	while(a[k]!="T1"&&a[k]!="T2")
	{
	    s=(s)(a[k])("_")
	    k++;
	}
	# if in techrep1 it writes to gc1 (removing last _ to s into s1)
	if(a[k]=="T1")
	{
	    s1=substr(s,1,(length(s)-1));
	    gc1[s1,a[k+1]]=$2;
	}
	else
	{
	    # if in techrep1 it writes to gc2 (removing last _ to s into s1)
	    if(a[k]=="T2")
	    {
		s1=substr(s,1,(length(s)-1));
		gc2[s1,a[k+1]]=$2;
	    }
	}
    }

    # reads the dupl level summary for each tech rep and read no
    # only _T1_1_val_1 or _T1_2_val_2 or _T2_1_val_1 or _T2_2_val_2 according to tech rep no and then read no
    while (getline < fileRef6 >0)
    {
	s="";
	s1="";
	split($1,a,"_");
	k=1;

	# until it hits T1 or T2 (indication of tech rep no) it continues to add info to s string corresponding to biorep id
	while(a[k]!="T1"&&a[k]!="T2")
	{
	    s=(s)(a[k])("_")
	    k++;
	}
	# if in techrep1 it writes to duplev1 (removing last _ to s into s1)
	if(a[k]=="T1")
	{
	    s1=substr(s,1,(length(s)-1));
	    duplev1[s1,a[k+1]]=$2;
	}
	else
	{
	    # if in techrep1 it writes to duplev2 (removing last _ to s into s1)
	    if(a[k]=="T2")
	    {
		s1=substr(s,1,(length(s)-1));
		duplev2[s1,a[k+1]]=$2;
	    }
	}
    }
}

NR==1{
    print "biorepid1", "biorepid2", "techrep1.read1.initnb", "techrep1.read2.initnb", "techrep2.read1.initnb", "techrep2.read2.initnb", "techrep1.read1.initbp.nb", "techrep1.read1.trimbp.nb", "techrep1.read1.trimbp.pcent", "techrep1.read2.initbp.nb", "techrep1.read2.trimbp.nb", "techrep1.read2.trimbp.pcent", "techrep2.read1.initbp.nb", "techrep2.read1.trimbp.nb", "techrep2.read1.trimbp.pcent", "techrep2.read2.initbp.nb", "techrep2.read2.trimbp.nb", "techrep2.read2.trimbp.pcent", "techrep1.read1.gcsummary", "techrep1.read2.gcsummary", "techrep2.read1.gcsummary", "techrep2.read2.gcsummary", "techrep1.read1.duplevsummary", "techrep1.read2.duplevsummary", "techrep2.read1.duplevsummary", "techrep2.read2.duplevsummary", $2, $3, $4, $5, $6, $7, $8, $9;
}

NR>=2{
    print $1, biorepid[fastqfile[$1]], readinit1[$1,1], readinit1[$1,2], readinit2[$1,1], readinit2[$1,2], bpinit1[$1,1], bptrim1[$1,1], pcenttrim1[$1,1], bpinit1[$1,2], bptrim1[$1,2], pcenttrim1[$1,2], bpinit2[$1,1], bptrim2[$1,1], pcenttrim2[$1,1], bpinit2[$1,2], bptrim2[$1,2], pcenttrim2[$1,2], gc1[$1,1], gc1[$1,2], gc2[$1,1], gc2[$1,2], duplev1[$1,1], duplev1[$1,2], duplev2[$1,1], duplev2[$1,2], $2, $3, $4, $5, $6, $7, $8, $9;
}
 
