# atac.qc.techrep.2.biorep.awk
# this script starts from a tsv file with header that has qc for atac seq technical replicates
# and outputs a tsv file with header that has the same qc but for biological replicates
# additionally it gives a unique and short biorep id and associates each biorep to all its techrep
# it is not easily generalisable since some metrics can be added, other averages, other like pass, warn or fail
# will be taken as the worst of them. We also assume that we only have two tech rep for each biorep here
# in particular in the end portion of the script

# example of usage
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp5/atacseq/qc.measures/firstrun.techrep.as.biorep
# meta=~/geneswitch/data/metadata/atacseq/wp5/fastq_tomotherinfo.ok.tsv
# pgm=~/fragencode/tools/multi/Scripts/atac.qc.techrep.2.biorep.awk
# time awk -v fileRef1=samplesheet.ok.csv -v fileRef2=$meta -f $pgm lid.read1.2.nb.bp.trim.nb.pcent.gc.duplev.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.peaknb.frip.tsv > biorepid.techreps.tiss.age.mid.mdiet.animal.id.name.gender.mbreed.fid.read1.2.nb.bp.trim.nb.pcent.gc.duplev.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.frip.tsv
# real	0m0.074s


# samplesheet.ok.csv
# group,replicate,fastq_1,fastq_2
# liver_fetus_lowfibre,1,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-941-1_CCGTGAAG-ATCCACTG-AHJK7WDSX2_L003_R1.fastq.gz,/work/project/geneswitch/data/reads/atacseq/sus_scrofa/wp5/complete/Liv-941-1_CCGTGAAG-ATCCACTG-AHJK7WDSX2_L003_R2.fastq.gz
# 633 (1 fields)

# $meta
# fastqfile	readtype	runnumber	tissue	age	motherid	diet	sampleid	samplename	sex	motherbreed	fatherid
# Skel-mus-Fetus-109-1_TTGGACTC-CTGCTTCC-BHMYNVDSX2_L004_R2.fastq.gz	R2	BHMYNVDSX2	skeletalmuscle	fetus	109	control	109.1	SSC_WU_WP5_FT_70dpf_109.1	male	TN70	E6211
# 1281 (12 fields)  *** $9":"$4":"$3 identify a techrep uniquely, while $9":"$4 identify a biorep uniquely

# lid.read1.2.nb.bp.trim.nb.pcent.gc.duplev.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.peaknb.frip.tsv
# techrep    read1.initnb    read2.initnb    read1.initbp.nb read1.trimbp.nb    read1.trimbp.pcent    read2.initbp.nb read2.trimbp.nb    read2.trimbp.pcent    read1.gcsummary read2.gcsummary    read1.duplevsummary    read2.duplevsummary clnread.tot    clnreadmap.nb    clnreadmap.pcent clnreadmt.nb    clnreadmt.pcent    frac.dupl.picard    peak.nb frip
# liver_fetus_control_R10    18718468    18718468    2807770200 1919447001    68.362    2807770200    1926371664    68.6086 PASS    PASS    PASS    PASS    24606658    24493104    99.5385 473292    1.92343    0.179214    35060    0.359275
# 633 (21 fields)

# output biorepid.techreps.tiss.age.mid.mdiet.animal.id.name.gender.mbreed.fid.read1.2.nb.bp.trim.nb.pcent.gc.duplev.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.frip.tsv
# biorep	techreps	tissue	age	motherid	motherdiet	animalid	animalname	gender	motherbreed	fatherid	read1.initnb	read2.initnb	read1.initbp.nb	read1.trimbp.nb	read1.trimbp.pcent	read2.initbp.nb	read2.trimbp.nb	read2.trimbp.pcent	read1.gcsummary	read2.gcsummary	read1.duplevsummary	read2.duplevsummary	clnread.tot	clnreadmap.nb	clnreadmap.pcent	clnreadmt.nb	clnreadmt.pcent	frac.dupl.picard	frip
# SSC_WU_WP5_Piglet_6447:skeletalmuscle	skelmus_piglet_highfibre_R38,skelmus_piglet_highfibre_R37,	skeletalmuscle	piglet	194	highfibre	6447	SSC_WU_WP5_Piglet_6447	female	TN60	E6021	72091288	72091288	10813693200	7031819374	65.027	10813693200	6992384383	64.6623	PASS	PASS	PASS	PASS	89828778	89318716	99.4322	2469796	2.74945	0.206196	0.226743
# 317 (30 fields)  *** many checks done in the readme of the output dir and seems fine

BEGIN{
    OFS="\t";
    # read the design file of the pipeline associating each tech rep id to its two fastq files, removing the absolute paths to the files
    while (getline < fileRef1 >0)
    {
	split($1,a,",");
	n1=split(a[3],b,"/");
	n2=split(a[4],c,"/");
	tid=a[1]"_R"a[2]
	trepid[b[n1]]=tid;
	trepid[c[n2]]=tid;
	fastq1[tid]=b[n1];
	fastq2[tid]=c[n2]
    }

    # read the metadata file of the reads to collect all possible metadata about the underlying samples
    # and also to associate each biorep to its two techreps
    # $9":"$4 (like SSC_WU_WP5_FT_70dpf_109.1:skeletalmuscle) is a unique identifier of a biorep, but
    # not very nice (mus_pig_F_TN70_high_E5918_865 nicer)
    # $9":"$4":"$3 (like SSC_WU_WP5_FT_70dpf_109.1:skeletalmuscle:BHMYNVDSX2) is a unique identifier of a techrep, but
    # not very nice (liver_fetus_control_R10 nicer)
    while (getline < fileRef2 >0)
    {
	# the info is the same for R1 and R2 so only read R1 rows
	if($2=="R1")
	{
	    # do not consider the read files from the metadata file that have not been processed by the pipeline
	    if(trepid[$1]!="")
	    {
		techrep[$9":"$4]=(techrep[$9":"$4])(trepid[$1])(",");
	    }
	    tiss[$9":"$4]=$4;
	    age[$9":"$4]=$5;
	    mid[$9":"$4]=$6;
	    diet[$9":"$4]=$7;
	    animid[$9":"$4]=$8;
	    animname[$9":"$4]=$9;
	    gender[$9":"$4]=$10;
	    mbreed[$9":"$4]=$11;
	    fid[$9":"$4]=$12;
	}
    }
}

{
    # read the input file to collect for each tech rep all its qc, here the tech rep is identified as liver_fetus_control_R10
    for(i=2; i<=NF; i++)
    {
	qc[$1,(i-1)]=$i;
    }
}

END{
    print "biorep", "techreps", "tissue", "age", "motherid", "motherdiet", "animalid", "animalname", "gender", "motherbreed", "fatherid", "read1.initnb", "read2.initnb", "read1.initbp.nb", "read1.trimbp.nb", "read1.trimbp.pcent", "read2.initbp.nb", "read2.trimbp.nb", "read2.trimbp.pcent", "read1.gcsummary", "read2.gcsummary", "read1.duplevsummary", "read2.duplevsummary", "clnread.tot", "clnreadmap.nb", "clnreadmap.pcent", "clnreadmt.nb", "clnreadmt.pcent", "frac.dupl.picard", "frip";
    for(biorep in techrep)
    {
	split(techrep[biorep],a,",");
	print biorep, techrep[biorep], tiss[biorep], age[biorep], mid[biorep], diet[biorep], animid[biorep], animname[biorep], gender[biorep], mbreed[biorep], fid[biorep], qc[a[1],1]+qc[a[2],1], qc[a[1],2]+qc[a[2],2], qc[a[1],3]+qc[a[2],3], qc[a[1],4]+qc[a[2],4], (qc[a[1],4]+qc[a[2],4])/(qc[a[1],3]+qc[a[2],3])*100, qc[a[1],6]+qc[a[2],6], qc[a[1],7]+qc[a[2],7], (qc[a[1],7]+qc[a[2],7])/(qc[a[1],6]+qc[a[2],6])*100, ((qc[a[1],9]=="FAIL"||qc[a[2],9]=="FAIL") ? "FAIL" : ((qc[a[1],9]=="WARN"||qc[a[2],9]=="WARN") ? "WARN" : "PASS")), ((qc[a[1],10]=="FAIL"||qc[a[2],10]=="FAIL") ? "FAIL" : ((qc[a[1],10]=="WARN"||qc[a[2],10]=="WARN") ? "WARN" : "PASS")), ((qc[a[1],11]=="FAIL"||qc[a[2],11]=="FAIL") ? "FAIL" : ((qc[a[1],11]=="WARN"||qc[a[2],11]=="WARN") ? "WARN" : "PASS")), ((qc[a[1],12]=="FAIL"||qc[a[2],12]=="FAIL") ? "FAIL" : ((qc[a[1],12]=="WARN"||qc[a[2],12]=="WARN") ? "WARN" : "PASS")), qc[a[1],13]+qc[a[2],13], qc[a[1],14]+qc[a[2],14], (qc[a[1],14]+qc[a[2],14])/(qc[a[1],13]+qc[a[2],13])*100, qc[a[1],16]+qc[a[2],16], (qc[a[1],16]+qc[a[2],16])/(qc[a[1],13]+qc[a[2],13])*100, (qc[a[1],18]+qc[a[2],18])/2, (qc[a[1],20]+qc[a[2],20])/2;
    }
} 
