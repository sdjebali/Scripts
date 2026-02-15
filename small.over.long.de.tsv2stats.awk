# small.over.long.de.tsv2stats.awk
# Takes a tsv file with header that has GS small genes overlapping GS long genes and involved in at least
# one differential analysis between stages (normal or ratio DE analysis) as rows, with indication of
# small gene id, small gene biotype, long gene id, long gene biotype in 4 1st fields and then in the last
# 14 fields, for each tissue (7) and stage pairs (2: s12 and s23), the indication of the result of the 1st
# DE analysis and then of the 2nd one separated by : (NA if not involved, ns if not significant, up if up
# and dw if down)
# And produces a stat tsv file that has on each row
# - the tissue and stage pair
# - whether we are considering only known, only unknown or all small genes
# - the number of small genes considered in the 1st de analysis
# - the number of small genes that are DE in the 1st de analysis
# - the % of those out of the ones considered
# - the same three numbers for the 2nd de analysis
# - the same three number for the ones in the 1st and 2nd de analysis (no consideration whether go in same direction)

# example
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# indir=~/fragencode/workspace/geneswitch/results
# pgm=~/fragencode/tools/multi/Scripts/small.over.long.de.tsv2stats.awk
# sp=gallus_gallus
# cd $indir/srnaseq/srna-rna-ratio/comparison.two.small.de.analyses/$sp
# time awk -f $pgm common.small.overlong.gnid.gnbt.deres12.eachtiss.stagepair.tsv > common.small.overlong.tisstagepair.gncat.considered.de.nb.pcent.1st.2nd.both.deanal.tsv
# real	0m0.137s


# input file  common.small.overlong.gnid.gnbt.deres12.eachtiss.stagepair.tsv
# smallgn.id	smallgn.bt	longgn.id	longgn.bt	ileum.stage1.stage2	ileum.stage2.stage3	cerebellum.stage1.stage2	cerebellum.stage2.stage3	kidney.stage1.stage2	kidney.stage2.stage3	liver.stage1.stage2	liver.stage2.stage3	lung.stage1.stage2	lung.stage2.stage3	muscle.stage1.stage2	muscle.stage2.stage3	skin.stage1.stage2	skin.stage2.stage3
# MIRNA_hairpin_0	miRNA	LOC_000000024161	mRNA	NA:NA	NA:NA	NA:NA	NA:NA	ns:ns	NA:NA	NA:NA	NA:NA	NA:NA	NA:NA	NA:NA	NA:NA	NA:NA	NA:NA
# 10175 (18 fields)

# output file common.small.overlong.tisstagepair.gncat.considered.de.nb.pcent.1st.2nd.both.deanal.tsv
# tisstagepair	smallgncat	1stDE.considered	1stDE.sig.nb	1stDE.sig.pcent	2ndDE.considered	2ndDE.sig.nb	2ndDE.sig.pcent	bothDE.considered	bothDE.sig.nb	bothDE.sig.pcent
# ileum.stage1.stage2	known	393	21	5.34351	393	20	5.08906	393	9	2.29008
# ileum.stage1.stage2	unknown	1471	117	7.95377	1471	2	0.135962	1471	1	0.067981
# ileum.stage1.stage2	all	1864	138	7.40343	1864	22	1.18026	1864	10	0.536481
# ileum.stage2.stage3	known	375	67	17.8667	375	137	36.5333	375	19	5.06667
# ileum.stage2.stage3	unknown	3144	899	28.5941	3144	36	1.14504	3144	17	0.540712
# ileum.stage2.stage3	all	3519	966	27.451	3519	173	4.91617	3519	36	1.02302
# 43 (11 fields)  *** ok format and nb rows since 1 for header and then 14 "sample"*3 small cat



BEGIN{
    # classify each small gene biotype into unknown and known (= if not unknown)
    unknown["long_RNA_overlapping"]=1;
    unknown["uncharacterized_small_RNAs"]=1;
}

# record the names of the tissue and stage pairs
NR==1{
    for(i=5; i<=NF; i++)
    {
	samp[i]=$i;
    }
}

# for each small look whether it is unknown or know and for each sample (in fact tissue and stage pair)
# increment the counters for total number of small considered in each analysis and in both and number
# of those that are DE so that we can compute the stats at the end
NR>=2{
    # computes the class of the small to have it in the double entry tables below
    c=(unknown[$2]==1 ? "u" : "k");
    for(i=5; i<=NF; i++)
    {
	s=samp[i];
	split($i,a,":");
	# if considered in the 1st analysis
	if(a[1]!="NA")
	{
	    tot1[c,s]++;
	    # if significant in the 1st analysis
	    if(a[1]!="ns")
	    {
		sig1[c,s]++;
	    }
	    # if also considered in the 2nd analysis
	    if(a[2]!="NA")
	    {
		tot3[c,s]++;
		# if significant in the 1st and 2nd analyses
		if(a[2]!="ns"&&a[1]!="ns")
		{
		    sig3[c,s]++;
		}
	    }
	}
	# if considered in the 2nd analysis
	if(a[2]!="NA")
	{
	    tot2[c,s]++;
	    # if significant in the 2nd analysis
	    if(a[2]!="ns")
	    {
		sig2[c,s]++;
	    }
	}
    }
}

END{
    OFS="\t";
    # print the header of the file
    print "tisstagepair", "smallgncat", "1stDE.considered", "1stDE.sig.nb", "1stDE.sig.pcent", "2ndDE.considered", "2ndDE.sig.nb", "2ndDE.sig.pcent", "bothDE.considered", "bothDE.sig.nb", "bothDE.sig.pcent";
    
    # for each "sample" print the stats for the know, the unknow and all small
    for(j=5; j<=(i-1); j++)
    {
	s=samp[j];
	totone=tot1["k",s]+tot1["u",s];
	sigone=sig1["k",s]+sig1["u",s];
	tottwo=tot2["k",s]+tot2["u",s];
	sigtwo=sig2["k",s]+sig2["u",s];
	totthree=tot3["k",s]+tot3["u",s];
	sigthree=sig3["k",s]+sig3["u",s];
	print s, "known", tot1["k",s], nn(sig1["k",s]), sig1["k",s]/tot1["k",s]*100, tot2["k",s], nn(sig2["k",s]), sig2["k",s]/tot2["k",s]*100, tot3["k",s], nn(sig3["k",s]), sig3["k",s]/tot3["k",s]*100;
	print s, "unknown", tot1["u",s], nn(sig1["u",s]), sig1["u",s]/tot1["u",s]*100, tot2["u",s], nn(sig2["u",s]), sig2["u",s]/tot2["u",s]*100, tot3["u",s], nn(sig3["u",s]), sig3["u",s]/tot3["u",s]*100;
	print s, "all", totone, sigone, sigone/totone*100, tottwo, sigtwo, sigtwo/tottwo*100, totthree, sigthree, sigthree/totthree*100; 	
    }
}

function nn(x){
    return (x=="" ? 0 : x)
}
