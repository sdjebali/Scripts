# add_codinginfo_about3ppartjunction.awk
# same as add_codinginfo_about5ppartjunction.awk but for the 3' part of the chimeric junction
# !!! TODO in the future !!!
# !!!! ideally we would want a single script that takes as input whether we are dealing with the 5' or the 3' part of the junction !!!
# !!! and that outputs a tsv file with 4 additional columns relative to the 5' or the 3' part of the junction !!!
# given:
########
# - a fileRef gff file containting the result of the overlap of the 3' part of chimeric junctions with cds exons
#   (in terms of nr lists of genes, transcripts, cds frames and cds coordinates)
# - a main input tsv file of junctions
# outputs a tsv file similar to the main input tsv file
#######################################################
# but with where the junctions that overlap cds of several genes and this gene is not present in $12 are eliminated and with 3 additional fields
# - threepcds which is a boolean saying whether the 3' part of the junction overlaps a cds  
# - threepcdsframe which is a boolean that if threepcds is 1 then says whether there is at least one cds in the same frame as the 3' part of the junction
# - threepokcdslist which is a list of such cds (and na if the list is empty)
# - threepallcdsok which is a boolean saying whether all the overlapped cds exons were in the same <<frame>> as the 3' part of the junction
#   by checking whether the distance between the cds acceptor site (AG) and the junction acceptor site is a multiple of 3
# !!! be careful the 3' part of the junction could be consistent with several cds exons, when all distances to AG (acceptor) of the junction are multiple of 3 !!!
# !!! without those cds exons being in the same frame with each other, so need to check later on, and also remove the cases where allcdsok is 0, see below !!!

# example
# cd ~/work/chimeras/sarcomas/crct/chimeras/analysis/merged.cases/1span1pe
# pgm=~/fragencode/tools/multi/Scripts/add_codinginfo_about3ppartjunction.awk
# time awk -v fileRef=tmp2.overcds.gff -f $pgm chimjunc.union.nopseudo.noercc.singlegneachside.pcgin5p.nbsamples.pcg.cds.orf.tsv > tmp

# input tmp2.overcds.gff
# chr1	.	.	1213034	1213093	.	-	.	junc: chr1_1223244_-:chr1_1213093_- gene: ENSG00000186827.10 nb_ov_cds_10: 1 list_cds_10: ENSG00000186827.10, nb_ov_cds_12: 1 list_cds_12: ENST00000379236.3, nb_ov_cds_22: 1 list_cds_22: 2, nb_ov_cds_24: 1 list_cds_24: chr1:1212992:1213093,
# chr1	.	.	1244448	1244497	.	-	.	junc: chr1_1267862_-:chr1_1244497_- gene: ENSG00000184163.3 nb_ov_cds_10: 1 list_cds_10: ENSG00000184163.3, nb_ov_cds_12: 1 list_cds_12: ENST00000330388.2, nb_ov_cds_22: 1 list_cds_22: 0, nb_ov_cds_24: 1 list_cds_24: chr1:1244381:1244497,
# 3361 (28 fields)

# input chimjunc.union.nopseudo.noercc.singlegneachside.pcgin5p.nbsamples.pcg.cds.orf.tsv 
# junc_id	beg	end	samechrstr	okgxorder	dist	type	ss1	ss2	gnlist1	gnlist2	gnname1	gnname2	gnbt1	gnbt2	LMS10R	LMS11R	LMS12RLMS13R	LMS14R	LMS15R	LMS16R	LMS18R1	LMS18R2	LMS18R3	LMS19R1	LMS19R1bis	LMS19R2	LMS19R3	LMS19R4	LMS19R5	LMS1R	LMS20R	LMS21R	LMS22R	LMS23R	LMS24R	LMS25RLMS26R	LMS27R	LMS28R	LMS29R	LMS2R	LMS30R	LMS31R	LMS32R	LMS33R	LMS34R	LMS35R	LMS36R	LMS37R	LMS38R	LMS39R	LMS3R	LMS40R	LMS41R	LMS42R	LMS43R	LMS44RLMS45R	LMS46R	LMS47R	LMS48R	LMS49R	LMS4R	LMS50R	LMS51R	LMS52R	LMS53R	LMS54R	LMS55R	LMS56R	LMS57R	LMS58R	LMS59R	LMS5R	LMS60R	LMS61R	LMS62R	LMS63RLMS64R	LMS65R	LMS66R	LMS67R	LMS68R	LMS69R	LMS6R	LMS7R1	LMS7R2	LMS7R3	LMS8R1	LMS8R2	LMS8R3	LMS9R1	LMS9R2	M961	M963	M964	M965	M966	M967	M968	M969	M970	M971	nbsamples	pcg	cds	orf	fiveppcg	fivepcds	fivepcdsframe	fivepokcdslist	fivepallcdsok	threeppcg
# chr16_21625115_+:chr16_21626947_+	21625041	21627021	1	1	1832	readthrough	GT	AG	ENSG00000197006.13	ENSG00000261596.2	METTL9	CTB-31N19.3	protein_coding	sense_intronic	0:0	0:0	0:0	0:0	0:0	0:0	0:0	5:3	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	11:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	8:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	6:0	0:0	0:0	0:0	0:0	21:3	0:0	8:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	6	0	0	0	1	1	1	chr16:21624931:21625115,	0	0
# 3362 (115 fields)
 
# output tmp
# junc_id	beg	end	samechrstr	okgxorder	dist	type	ss1	ss2	gnlist1	gnlist2	gnname1	gnname2	gnbt1	gnbt2	LMS10R	LMS11R	LMS12RLMS13R	LMS14R	LMS15R	LMS16R	LMS18R1	LMS18R2	LMS18R3	LMS19R1	LMS19R1bis	LMS19R2	LMS19R3	LMS19R4	LMS19R5	LMS1R	LMS20R	LMS21R	LMS22R	LMS23R	LMS24R	LMS25RLMS26R	LMS27R	LMS28R	LMS29R	LMS2R	LMS30R	LMS31R	LMS32R	LMS33R	LMS34R	LMS35R	LMS36R	LMS37R	LMS38R	LMS39R	LMS3R	LMS40R	LMS41R	LMS42R	LMS43R	LMS44RLMS45R	LMS46R	LMS47R	LMS48R	LMS49R	LMS4R	LMS50R	LMS51R	LMS52R	LMS53R	LMS54R	LMS55R	LMS56R	LMS57R	LMS58R	LMS59R	LMS5R	LMS60R	LMS61R	LMS62R	LMS63RLMS64R	LMS65R	LMS66R	LMS67R	LMS68R	LMS69R	LMS6R	LMS7R1	LMS7R2	LMS7R3	LMS8R1	LMS8R2	LMS8R3	LMS9R1	LMS9R2	M961	M963	M964	M965	M966	M967	M968	M969	M970	M971	nbsamples	pcg	cds	orf	fiveppcg	fivepcds	fivepcdsframe	fivepokcdslist	fivepallcdsok	threeppcg	threepcds	threepcdsframe	threepokcdslist	threepallcdsok
# chr16_21625115_+:chr16_21626947_+	21625041	21627021	1	1	1832	readthrough	GT	AG	ENSG00000197006.13	ENSG00000261596.2	METTL9	CTB-31N19.3	protein_coding	sense_intronic	0:0	0:0	0:0	0:0	0:0	0:0	0:0	5:3	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	11:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	8:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	6:0	0:0	0:0	0:0	0:0	21:3	0:0	8:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	6	0	0	0	1	1	1	chr16:21624931:21625115,	0	0	0	0	na	0
# 3344 (119 fields)  *** 4 more columns, fine  *** 18 junctions are lost due to overlapping cds of diff genes

BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	ko[$10]=0;
	if($16!=".")
	{
	    split($16,a,",");
	    k=1;
	    while(ko[$10]!=1&&a[k]!="")
	    {
		if(a[k]!=$12)
		{
		    ko[$10]=1;
		}
		k++;
	    }
	    # at this point $28 should not be empty
	    split($28,b,",");
	    l=1;
	    while(b[l]!="")
	    {
		split(b[l],c,":");
		if($7=="+")
		{
		    # distance between the beginings
		    accdist[$10]=c[2]-$4;
		}
		else
		{
		    # distance between the ends
		    accdist[$10]=c[3]-$5;
		}
		if(((accdist[$10])%3)==0)
		{
		    cdslist[$10]=(cdslist[$10])(b[l])(",");
		}
		else
		{
		    allcdsko[$10]=1;
		}
		l++;
	    }
	}
	if(ko[$10]==0)
	{
	    threepcds[$10]=($16!="." ? 1 : 0);
	    threepcdsframe[$10]=((($16!=".")&&(cdslist[$10]!="")) ? 1 : 0);
	}
    }
}

NR==1{
    print $0, "threepcds", "threepcdsframe", "threepokcdslist", "threepallcdsok";
}

NR>=2&&ko[$1]==0{
    print $0, threepcds[$1], threepcdsframe[$1], (cdslist[$1]!="" ? cdslist[$1] : "na"), ((threepcds[$1]&&threepcdsframe[$1]&&(allcdsko[$1]!=1)) ? 1 : 0);
}
