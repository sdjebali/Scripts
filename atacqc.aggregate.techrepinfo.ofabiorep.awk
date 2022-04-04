# atacqc.aggregate.techrepinfo.ofabiorep.awk
# this script starts from a tsv file with header that has qc for atac seq biological replicates
# and outputs a tsv file with header that has the same qc but aggregating the info we have for
# the two tech rep of the biorep on the row of this biorep, this is the case for
# - nb init read (addition of the reads of the two tech rep for each read no)
# - nb init bp, nb trimmed bp, pcent trimmed bp (addition for the two first, recomputation for the 3rd for each read no)
# - gc summary (worst of the two for each read no)
# - duplev summary (worst of the two for each read no)
# it also adds many metadata info about the biorep using fileRef=$meta file (type fastq2motherinfo)
# Note: it is not easily generalisable since some metrics can be added, other averaged, other like pass, warn or fail
# will be taken as the worst of them. We also assume that we only have two tech rep for each biorep here in particular
# in the end portion of the script

# Example of usage
##################
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp5/atacseq/qc.measures
# meta=~/geneswitch/data/metadata/atacseq/wp5/fastq_tomotherinfo.ok.tsv
# pgm=~/fragencode/tools/multi/Scripts/atacqc.aggregate.techrepinfo.ofabiorep.awk
# time awk -v fileRef=$meta -f $pgm lid1.2.techrep1.2.read1.2.nb.bp.trim.nb.pcent.gc.duplev.techrep1.2.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.peaknb.frip.tsv > biorepid1.2.tiss.age.mid.mdiet.animal.id.name.gender.mbreed.fid.read1.2.nb.bp.trim.nb.pcent.gc.duplev.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.frip.tsv
# real	0m0.041s


# fileRef=$meta
# fastqfile	readtype	runnumber	tissue	age	motherid	diet	sampleid	samplename	sex	motherbreed	fatherid
# Skel-mus-Fetus-109-1_TTGGACTC-CTGCTTCC-BHMYNVDSX2_L004_R2.fastq.gz	R2	BHMYNVDSX2	skeletalmuscle	fetus	109	control	109.1	SSC_WU_WP5_FT_70dpf_109.1	male	TN70	E6211
# 1281 (12 fields)  *** $9":"$4 identify a biorep uniquely

# main input=lid1.2.techrep1.2.read1.2.nb.bp.trim.nb.pcent.gc.duplev.techrep1.2.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.peaknb.frip.tsv 
# biorepid1	biorepid2	techrep1.read1.initnb	techrep1.read2.initnb	techrep2.read1.initnb	techrep2.read2.initnb	techrep1.read1.initbp.nb	techrep1.read1.trimbp.nb	techrep1.read1.trimbp.pcent	techrep1.read2.initbp.nb	techrep1.read2.trimbp.nb	techrep1.read2.trimbp.pcent	techrep2.read1.initbp.nb	techrep2.read1.trimbp.nb	techrep2.read1.trimbp.pcent	techrep2.read2.initbp.nb	techrep2.read2.trimbp.nb	techrep2.read2.trimbp.pcent	techrep1.read1.gcsummary	techrep1.read2.gcsummary	techrep2.read1.gcsummary	techrep2.read2.gcsummary	techrep1.read1.duplevsummary	techrep1.read2.duplevsummary	techrep2.read1.duplevsummary	techrep2.read2.duplevsummary	clnread.tot	clnreadmap.nb	clnreadmap.pcent	clnreadmt.nb	clnreadmt.pcent	frac.dupl.picard	peak.nb	frip
# liver_fetus_control_R10	SSC_WU_WP5_FT_70dpf_175.2:liver	68268904	68268904	65985728	65985728	10240335600	6449723020	62.9835	10240335600	6475616541	63.2364	9897859200	6266643954	63.3131	9897859200	6294061983	63.5901	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	172461956	170078000	98.6177	1476758	0.85628	0.221747	76791	0.437052
# 321 (34 fields)

# main output=biorepid1.2.tiss.age.mid.mdiet.animal.id.name.gender.mbreed.fid.read1.2.nb.bp.trim.nb.pcent.gc.duplev.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.frip.tsv
# biorepid1	biorepid2	tissue	age	motherid	motherdiet	animalid	animalname	gender	motherbreed	fatherid	read1.initnb	read2.initnb	read1.initbp.nb	read1.trimbp.nb	read1.trimbp.pcent	read2.initbp.nb	read2.trimbp.nb	read2.trimbp.pcent	read1.gcsummary	read2.gcsummary	read1.duplevsummary	read2.duplevsummary	clnread.tot	clnreadmap.nb	clnreadmap.pcent	clnreadmt.nb	clnreadmt.pcent	frac.dupl.picard	peak.nb	frip
# liver_fetus_control_R10	SSC_WU_WP5_FT_70dpf_175.2:liver	liver	fetus	175	control	175.2	SSC_WU_WP5_FT_70dpf_175.2	male	TN60	E6272	134254632	134254632	20138194800	12716366974	63.1455	20138194800	12769678524	63.4102	PASS	PASS	PASS	PASS	172461956	170078000	98.6177	1476758	0.85628	0.221747	76791	0.437052
# 321 (31 fields)


BEGIN{
    OFS="\t";
    
    # read the metadata file of the reads to collect all possible metadata about the biorep we have (its 2nd id is of the type
    # SSC_WU_WP5_FT_70dpf_109.1:skeletalmuscle and corresponds to $9":"$4 in the metadata file
    # fastqfile	readtype	runnumber	tissue	age	motherid	diet	sampleid	samplename	sex	motherbreed	fatherid
    # Skel-mus-Fetus-109-1_TTGGACTC-CTGCTTCC-BHMYNVDSX2_L004_R2.fastq.gz	R2	BHMYNVDSX2	skeletalmuscle	fetus	109	control	109.1	SSC_WU_WP5_FT_70dpf_109.1	male	TN70	E6211
    # 1281 (12 fields)
    while (getline < fileRef >0)
    {
	# the info is the same for R1 and R2 so only read R1 rows to retrieve the info we want 
	if($2=="R1")
	{
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

# print the header
# all metadata info at the begining, then the aggregated values and finally the ones that we already have per biorep
NR==1{
    print $1, $2, "tissue", "age", "motherid", "motherdiet", "animalid", "animalname", "gender", "motherbreed", "fatherid", "read1.initnb", "read2.initnb", "read1.initbp.nb", "read1.trimbp.nb", "read1.trimbp.pcent", "read2.initbp.nb", "read2.trimbp.nb", "read2.trimbp.pcent", "read1.gcsummary", "read2.gcsummary", "read1.duplevsummary", "read2.duplevsummary", $27, $28, $29, $30, $31, $32, $33, $34;
}

# print the body
# all metadata info at the begining and then the aggregated values and finally the ones that we already have per biorep
NR>=2{
    print $1, $2, tiss[$2], age[$2], mid[$2], diet[$2], animid[$2], animname[$2], gender[$2], mbreed[$2], fid[$2], $3+$5, $4+$6, $7+$13, $8+$14, ($8+$14)/($7+$13)*100, $10+$16, $11+$17, ($11+$17)/($10+$16)*100, agg($19, $21), agg($20, $22), agg($23, $25), agg($24, $26), $27, $28, $29, $30, $31, $32, $33, $34;
}

# aggregation function for the fastqc summaries of type pass, warn, fail, we take the worst between the two techrep
function agg(x,y)
{
    return ((x=="FAIL"||y=="FAIL") ? "FAIL" : ((x=="WARN"||y=="WARN") ? "WARN" : "PASS"));
}

