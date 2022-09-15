# gather.qc.measures.atacseq.2.awk
# similar to gather.qc.measures.atacseq.awk but for wp2 (the first one was made for wp5) = accepts 4 tech rep instead of 2
# and does not use external metadata information (so uses 4 fileRef files instead of 6 which is much simpler)

# like gather.qc.measures.atacseq.awk this script gathers several kinds of qc metrix distributed in 5 different tsv files and derived from the output of
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
# - 4 fileRef files
#   * fileRef1=lid.nbinitreads.tsv that has for each technical replicate and read number the initial number of reads (8*N rows)
#   * fileRef2=lid.nbbpinit.aftertrim.nb.pcent.tsv that has for each techrep and read no the init nb of bp, the one after trimming and the % it represents (8*N rows)
#   * fileRef3=lid.gccontent.aftertrimming.summary.tsv that has for each techrep and read no whether it passes the fastqc GC content test (8*N rows)
#   * fileRef4=lid.dupllevel.aftertrimming.summary.tsv that has for each techrep and read no whether it passes the fastqc duplication level test (8*N rows)

# It will output a tsv file with header that has 57 fields in this order
########################################################################
#   * biorep id from the pipeline
#   * techrep1.read1.initnb
#   * techrep1.read2.initnb
#   * techrep2.read1.initnb
#   * techrep2.read2.initnb
#   * techrep3.read1.initnb
#   * techrep3.read2.initnb
#   * techrep4.read1.initnb
#   * techrep4.read2.initnb
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
#   * techrep3.read1.initbp.nb
#   * techrep3.read1.trimbp.nb
#   * techrep3.read1.trimbp.pcent
#   * techrep3.read2.initbp.nb
#   * techrep3.read2.trimbp.nb
#   * techrep3.read2.trimbp.pcent
#   * techrep4.read1.initbp.nb
#   * techrep4.read1.trimbp.nb
#   * techrep4.read1.trimbp.pcent
#   * techrep4.read2.initbp.nb
#   * techrep4.read2.trimbp.nb
#   * techrep4.read2.trimbp.pcent
#   * techrep1.read1.gcsummary
#   * techrep1.read2.gcsummary
#   * techrep2.read1.gcsummary
#   * techrep2.read2.gcsummary
#   * techrep3.read1.gcsummary
#   * techrep3.read2.gcsummary
#   * techrep4.read1.gcsummary
#   * techrep4.read2.gcsummary
#   * techrep1.read1.duplevsummary
#   * techrep1.read2.duplevsummary
#   * techrep2.read1.duplevsummary
#   * techrep2.read2.duplevsummary
#   * techrep3.read1.duplevsummary
#   * techrep3.read2.duplevsummary
#   * techrep4.read1.duplevsummary
#   * techrep4.read2.duplevsummary
#   * the same 9 first rows as the input file
# in a second step we can add some values, average others or select the worst in pass, warn and fail
# but here we just report the info about each tech rep of a biorep
# TODO: allow a biorep to have more or less than 4 techrep

# example
# cd ~/fragencode/workspace/sdjebali/geneswitch/analysis/atacseq/qc.measures
# pgm=~/fragencode/tools/multi/Scripts/gather.qc.measures.atacseq.2.awk
# time paste lid.totreads.mapped.nb.pcent.tsv lid.totmappedreads.mt.nb.pcent.tsv lid.fracduplicate.tsv lid.peaks.nb.frip.tsv | awk 'BEGIN{OFS="\t"; print "biorep", "clnread.tot", "clnreadmap.nb", "clnreadmap.pcent", "clnreadmt.nb", "clnreadmt.pcent", "frac.dupl.picard", "peak.nb", "frip"} {split($9,a,"_T1"); if($1==$5&&$1==a[1]&&$1==$11){print $1, $2, $3, $4, $7, $8, $10, $12, $13}}' | awk -v fileRef1=lid.nbinitreads.tsv -v fileRef2=lid.nbbpinit.aftertrim.nb.pcent.tsv -v fileRef3=lid.gccontent.aftertrimming.summary.tsv -v fileRef4=lid.dupllevel.aftertrimming.summary.tsv -f $pgm > lid1.2.techrep1.2.read1.2.nb.bp.trim.nb.pcent.gc.duplev.techrep1.2.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.peaknb.frip.tsv
# real	0m0.060s

# fileRef1=lid.nbinitreads.tsv
# pig_cerebellum_stage1_R1_T1_1	28428160
# 672 (2 fields)  *** always T1 and then T2, T3 and T4 according to tech rep no and after that either 1 or 2 according to read no

# fileRef2=lid.nbbpinit.aftertrim.nb.pcent.tsv
# pig_cerebellum_stage1_R1_T1_1	4264224000	4218584327	98.9297
# 672 (4 fields) *** always T1 and then T2, T3 and T4 according to tech rep no and after that either 1 or 2 according to read no

# fileRef3=lid.gccontent.aftertrimming.summary.tsv
# pig_cerebellum_stage1_R1_T1_1_val_1	PASS
# 672 (2 fields)  *** only _TX_1_val_1 or _TX_2_val_2 for X in 1 to 4 according to tech rep no and then read no

# fileRef4=lid.dupllevel.aftertrimming.summary.tsv
# pig_cerebellum_stage1_R1_T1_1_val_1	PASS
# 672 (2 fields)  *** only _TX_1_val_1 or _TX_2_val_2 for X in 1 to 4 according to tech rep no and then read no

# main input file
# biorep	clnread.tot	clnreadmap.nb	clnreadmap.pcent	clnreadmt.nb	clnreadmt.pcent	frac.dupl.picard	peak.nb	frip
# pig_cerebellum_stage1_R1	87593124	87580540	99.9856	521132	0.594946	0.186957	48654	0.218651
# 85 (9 fields)

# main output file
# biorepid	techrep1.read1.initnb	techrep1.read2.initnb	techrep2.read1.initnb	techrep2.read2.initnb	techrep3.read1.initnb	techrep3.read2.initnb	techrep4.read1.initnb	techrep4.read2.initnb	techrep1.read1.initbp.nb	techrep1.read1.trimbp.nb	techrep1.read1.trimbp.pcent	techrep1.read2.initbp.nb	techrep1.read2.trimbp.nb	techrep1.read2.trimbp.pcent	techrep2.read1.initbp.nb	techrep2.read1.trimbp.nb	techrep2.read1.trimbp.pcent	techrep2.read2.initbp.nb	techrep2.read2.trimbp.nb	techrep2.read2.trimbp.pcent	techrep3.read1.initbp.nb	techrep3.read1.trimbp.nb	techrep3.read1.trimbp.pcent	techrep3.read2.initbp.nb	techrep3.read2.trimbp.nb	techrep3.read2.trimbp.pcent	techrep4.read1.initbp.nb	techrep4.read1.trimbp.nb	techrep4.read1.trimbp.pcent	techrep4.read2.initbp.nb	techrep4.read2.trimbp.nb	techrep4.read2.trimbp.pcent	techrep1.read1.gcsummary	techrep1.read2.gcsummary	techrep2.read1.gcsummary	techrep2.read2.gcsummary	techrep3.read1.gcsummary	techrep3.read2.gcsummary	techrep4.read1.gcsummary	techrep4.read2.gcsummary	techrep1.read1.duplevsummary	techrep1.read2.duplevsummary	techrep2.read1.duplevsummary	techrep2.read2.duplevsummary	techrep3.read1.duplevsummary	techrep3.read2.duplevsummary	techrep4.read1.duplevsummary	techrep4.read2.duplevsummary	clnread.tot	clnreadmap.nbclnreadmap.pcent	clnreadmt.nb	clnreadmt.pcent	frac.dupl.picard	peak.nb	frip
# pig_cerebellum_stage1_R1	28428160	28428160	28000872	28000872	27885387	27885387	21747882	21747882	4264224000	4218584327	98.9297	4264224000	4201252796	98.5233	4200130800	4155302704	98.9327	4200130800	4137878323	98.5178	4182808050	4137268580	98.9113	4182808050	4119363407	98.4832	4200130800	3223169199	98.8041	3262182300	3207389535	98.3204	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	PASS	87593124	87580540	99.9856	521132	0.594946	0.186957	48654	0.218651
# 85 (57 fields)

BEGIN{
    OFS="\t";

    # reads the file of initial read number in each tech rep and read no
    while (getline < fileRef1 >0)
    {
	s="";
	s1="";
	split($1,a,"_");
	k=1;
	
	# until it hits TX (indication of tech rep no X) it continues to add info to s string corresponding to biorep id
	while(a[k]!="T1"&&a[k]!="T2"&&a[k]!="T3"&&a[k]!="T4")
	{
	    s=(s)(a[k])("_");
	    k++;
	}
	
	# if in techrepX it writes to readinit[X,...] (removing last _ to s into s1)
	if(a[k]=="T1"||a[k]=="T2"||a[k]=="T3"||a[k]=="T4")
	{
	    s1=substr(s,1,(length(s)-1));
	    split(a[k],b,"T");
	    readinit[b[2],s1,a[k+1]]=$2;
	}
    }

    # reads the file of init bp, bp after trimming and percent of this
    while (getline < fileRef2 >0)
    {
	s="";
	s1="";
	split($1,a,"_");
	k=1;

	# until it hits TX (indication of tech rep no X) it continues to add info to s string corresponding to biorep id
	while(a[k]!="T1"&&a[k]!="T2"&&a[k]!="T3"&&a[k]!="T4")
	{
	    s=(s)(a[k])("_");
	    k++;
	}
	
	# if in techrepX it writes to arrayname[X,...] (removing last _ to s into s1)
	if(a[k]=="T1"||a[k]=="T2"||a[k]=="T3"||a[k]=="T4")
	{
	    s1=substr(s,1,(length(s)-1));
	    split(a[k],b,"T");
	    bpinit[b[2],s1,a[k+1]]=$2;
	    bptrim[b[2],s1,a[k+1]]=$3;
	    pcenttrim[b[2],s1,a[k+1]]=$4;
	}
    }

    # reads the gc summary for each tech rep and read no
    # only _TX_1_val_1 or _TX_2_val_2 according to tech rep no X and then read no
    while (getline < fileRef3 >0)
    {
	s="";
	s1="";
	split($1,a,"_");
	k=1;

	# until it hits TX (indication of tech rep no X) it continues to add info to s string corresponding to biorep id
	while(a[k]!="T1"&&a[k]!="T2"&&a[k]!="T3"&&a[k]!="T4")
	{
	    s=(s)(a[k])("_");
	    k++;
	}
	
	# if in techrepX it writes to gc[X,...] (removing last _ to s into s1)
	if(a[k]=="T1"||a[k]=="T2"||a[k]=="T3"||a[k]=="T4")
	{
	    s1=substr(s,1,(length(s)-1));
	    split(a[k],b,"T");
	    gc[b[2],s1,a[k+1]]=$2;
	}
    }

    # reads the dupl level summary for each tech rep and read no
    # only _TX_1_val_1 or _TX_2_val_2 according to tech rep no X and then read no
    while (getline < fileRef4 >0)
    {
	s="";
	s1="";
	split($1,a,"_");
	k=1;

	# until it hits TX (indication of tech rep no X) it continues to add info to s string corresponding to biorep id
	while(a[k]!="T1"&&a[k]!="T2"&&a[k]!="T3"&&a[k]!="T4")
	{
	    s=(s)(a[k])("_");
	    k++;
	}
	
	# if in techrepX it writes to gc[X,...] (removing last _ to s into s1)
	if(a[k]=="T1"||a[k]=="T2"||a[k]=="T3"||a[k]=="T4")
	{
	    s1=substr(s,1,(length(s)-1));
	    split(a[k],b,"T");
	    duplev[b[2],s1,a[k+1]]=$2;
	}
    }
}

NR==1{
    print "biorepid", "techrep1.read1.initnb", "techrep1.read2.initnb", "techrep2.read1.initnb", "techrep2.read2.initnb", "techrep3.read1.initnb", "techrep3.read2.initnb", "techrep4.read1.initnb", "techrep4.read2.initnb", "techrep1.read1.initbp.nb", "techrep1.read1.trimbp.nb", "techrep1.read1.trimbp.pcent", "techrep1.read2.initbp.nb", "techrep1.read2.trimbp.nb", "techrep1.read2.trimbp.pcent", "techrep2.read1.initbp.nb", "techrep2.read1.trimbp.nb", "techrep2.read1.trimbp.pcent", "techrep2.read2.initbp.nb", "techrep2.read2.trimbp.nb", "techrep2.read2.trimbp.pcent", "techrep3.read1.initbp.nb", "techrep3.read1.trimbp.nb", "techrep3.read1.trimbp.pcent", "techrep3.read2.initbp.nb", "techrep3.read2.trimbp.nb", "techrep3.read2.trimbp.pcent", "techrep4.read1.initbp.nb", "techrep4.read1.trimbp.nb", "techrep4.read1.trimbp.pcent", "techrep4.read2.initbp.nb", "techrep4.read2.trimbp.nb", "techrep4.read2.trimbp.pcent", "techrep1.read1.gcsummary", "techrep1.read2.gcsummary", "techrep2.read1.gcsummary", "techrep2.read2.gcsummary", "techrep3.read1.gcsummary", "techrep3.read2.gcsummary", "techrep4.read1.gcsummary", "techrep4.read2.gcsummary", "techrep1.read1.duplevsummary", "techrep1.read2.duplevsummary", "techrep2.read1.duplevsummary", "techrep2.read2.duplevsummary", "techrep3.read1.duplevsummary", "techrep3.read2.duplevsummary", "techrep4.read1.duplevsummary", "techrep4.read2.duplevsummary", $2, $3, $4, $5, $6, $7, $8, $9;
}

NR>=2{
    print $1, readinit[1,$1,1], readinit[1,$1,2], readinit[2,$1,1], readinit[2,$1,2], readinit[3,$1,1], readinit[3,$1,2], readinit[4,$1,1], readinit[4,$1,2], bpinit[1,$1,1], bptrim[1,$1,1], pcenttrim[1,$1,1], bpinit[1,$1,2], bptrim[1,$1,2], pcenttrim[1,$1,2], bpinit[2,$1,1], bptrim[2,$1,1], pcenttrim[2,$1,1], bpinit[2,$1,2], bptrim[2,$1,2], pcenttrim[2,$1,2], bpinit[3,$1,1], bptrim[3,$1,1], pcenttrim[3,$1,1], bpinit[3,$1,2], bptrim[3,$1,2], pcenttrim[3,$1,2], bpinit[2,$1,1], bptrim[4,$1,1], pcenttrim[4,$1,1], bpinit[4,$1,2], bptrim[4,$1,2], pcenttrim[4,$1,2], gc[1,$1,1], gc[1,$1,2], gc[2,$1,1], gc[2,$1,2], gc[3,$1,1], gc[3,$1,2], gc[4,$1,1], gc[4,$1,2], duplev[1,$1,1], duplev[1,$1,2], duplev[2,$1,1], duplev[2,$1,2], duplev[3,$1,1], duplev[3,$1,2], duplev[4,$1,1], duplev[4,$1,2], $2, $3, $4, $5, $6, $7, $8, $9;
}
 
