# add_codinginfo_about5ppartjunction.awk
# given
# - a fileRef gff file containting the result of the overlap of the 5' part of chimeric junctions with cds exons
#   (in terms of nr lists of genes, transcripts, cds frames and cds coordinates)
# - a main input tsv file of junctions
# outputs a tsv file similar to the main input tsv file but with where the junctions that overlap cds of several genes
# and this gene is not present in $12 are eliminated and with 3 additional fields
# - fivepcds which is a boolean saying whether the 5' part of the junction overlaps a cds  
# - fivepcdsframe which is a boolean that if fivepcds is 1 then says whether there is at least one cds in the same frame as the 5' part of the junction
# - fivepokcdslist which is a list of such cds (and na if the list is empty)
# - fivepallcdsok which is a boolean saying whether all the overlapped cds exons were in the same frame as the 5' part of the junction
#    by checking whether the distance between the cds donor site (GT) and the junction donor site is a multiple of 3
# !!! be careful the 5' part of the junction could be consistent with several cds exons, when all distances to GT of the junction are multiple of 3 !!!
# !!! without those cds exons being in the same frame with each other, so need to check later on, and also remove the cases where allcdsok is 0, see below !!!

# example
# cd ~/work/chimeras/sarcomas/crct/chimeras/analysis/merged.cases/gene.anything
# pgm=~/fragencode/tools/multi/Scripts/add_codinginfo_about5ppartjunction.awk
# time awk -v fileRef=tmp1.overcds.gff -f $pgm alljunctions_from2stag.tsv > tmp

# input tmp1.overcds.gff
# chr1	.	.	935856	935896	.	+	.	junc: chr1_935896_+:chr1_936154_+ gene: ENSG00000187634.11 nb_ov_cds_10: 1 list_cds_10: ENSG00000187634.11, nb_ov_cds_12: 9 list_cds_12: ENST00000341065.8,ENST00000342066.7,ENST00000420190.6,ENST00000616016.4,ENST00000616125.4,ENST00000617307.4,ENST00000618181.4,ENST00000618779.4,ENST00000622503.4, nb_ov_cds_22: 1 list_cds_22: 1, nb_ov_cds_24: 1 list_cds_24: chr1:935772:935896,
# chr1	.	.	944386	944446	.	+	.	junc: chr1_944446_+:chr2_130231847_- gene: ENSG00000187634.11 nb_ov_cds_10: 0 list_cds_10: . nb_ov_cds_12: 0 list_cds_12: . nb_ov_cds_22: 0 list_cds_22: . nb_ov_cds_24: 0 list_cds_24: .
# 40948 (28 fields)

# input alljunctions_from2stag.tsv
# junc_id	beg	end	samechrstr	okgxorder	dist	type	ss1	ss2	gnlist1	gnlist2	gnname1	gnname2	gnbt1	gnbt2	LMS10R	LMS11R	LMS12RLMS13R	LMS14R	LMS15R	LMS16R	LMS18R1	LMS18R2	LMS18R3	LMS19R1	LMS19R1bis	LMS19R2	LMS19R3	LMS19R4	LMS19R5	LMS1R	LMS20R	LMS21R	LMS22R	LMS23R	LMS24R	LMS25RLMS26R	LMS27R	LMS28R	LMS29R	LMS2R	LMS30R	LMS31R	LMS32R	LMS33R	LMS34R	LMS35R	LMS36R	LMS37R	LMS38R	LMS39R	LMS3R	LMS40R	LMS41R	LMS42R	LMS43R	LMS44RLMS45R	LMS46R	LMS47R	LMS48R	LMS49R	LMS4R	LMS50R	LMS51R	LMS52R	LMS53R	LMS54R	LMS55R	LMS56R	LMS57R	LMS58R	LMS59R	LMS5R	LMS60R	LMS61R	LMS62R	LMS63RLMS64R	LMS65R	LMS66R	LMS67R	LMS68R	LMS69R	LMS6R	LMS7R1	LMS7R2	LMS7R3	LMS8R1	LMS8R2	LMS8R3	LMS9R1	LMS9R2	M961	M963	M964	M965	M966	M967	M968	M969	M970	M971	nbsamp	class	fiveppcg
# chr19_41320470_+:chr14_72895852_+	41320411	72895867	0	na	na	UNANNOTATED	GT	AG	ENSG00000142039.3	unannotated	CCDC97	na	protein_coding	na	0:0	0:0	1:1	0:0	0:0	2:1	0:0	0:0	4:1	5:1	12:1	0:0	6:1	7:1	1:1	7:1	1:1	3:1	0:0	2:1	0:0	1:1	2:1	2:1	3:1	6:1	0:0	4:1	0:0	1:1	2:1	2:1	0:0	0:0	0:0	0:0	0:0	1:1	0:0	0:0	3:1	1:1	1:1	3:1	1:1	2:1	0:0	0:0	1:1	1:1	0:0	2:1	0:0	2:1	0:0	0:0	0:0	0:0	0:0	1:1	1:1	2:1	1:1	1:1	2:1	0:0	0:0	1:1	0:0	4:1	1:1	0:0	3:1	3:1	1:1	2:1	0:0	2:1	0:0	1:1	2:2	1:1	3:2	6:1	3:2	1:1	0:0	7:1	11:2	0:0	55	interchromosomal	1
# 40949 (108 fields)
 
# output tmp
# junc_id	beg	end	samechrstr	okgxorder	dist	type	ss1	ss2	gnlist1	gnlist2	gnname1	gnname2	gnbt1	gnbt2	LMS10R	LMS11R	LMS12RLMS13R	LMS14R	LMS15R	LMS16R	LMS18R1	LMS18R2	LMS18R3	LMS19R1	LMS19R1bis	LMS19R2	LMS19R3	LMS19R4	LMS19R5	LMS1R	LMS20R	LMS21R	LMS22R	LMS23R	LMS24R	LMS25RLMS26R	LMS27R	LMS28R	LMS29R	LMS2R	LMS30R	LMS31R	LMS32R	LMS33R	LMS34R	LMS35R	LMS36R	LMS37R	LMS38R	LMS39R	LMS3R	LMS40R	LMS41R	LMS42R	LMS43R	LMS44RLMS45R	LMS46R	LMS47R	LMS48R	LMS49R	LMS4R	LMS50R	LMS51R	LMS52R	LMS53R	LMS54R	LMS55R	LMS56R	LMS57R	LMS58R	LMS59R	LMS5R	LMS60R	LMS61R	LMS62R	LMS63RLMS64R	LMS65R	LMS66R	LMS67R	LMS68R	LMS69R	LMS6R	LMS7R1	LMS7R2	LMS7R3	LMS8R1	LMS8R2	LMS8R3	LMS9R1	LMS9R2	M961	M963	M964	M965	M966	M967	M968	M969	M970	M971	nbsamp	class	fiveppcg	fivepcds	fivepcdsframe	fivepokcdslist	fivepallcdsok
# chr19_41320470_+:chr14_72895852_+	41320411	72895867	0	na	na	UNANNOTATED	GT	AG	ENSG00000142039.3	unannotated	CCDC97	na	protein_coding	na	0:0	0:0	1:1	0:0	0:0	2:1	0:0	0:0	4:1	5:1	12:1	0:0	6:1	7:1	1:1	7:1	1:1	3:1	0:0	2:1	0:0	1:1	2:1	2:1	3:1	6:1	0:0	4:1	0:0	1:1	2:1	2:1	0:0	0:0	0:0	0:0	0:0	1:1	0:0	0:0	3:1	1:1	1:1	3:1	1:1	2:1	0:0	0:0	1:1	1:1	0:0	2:1	0:0	2:1	0:0	0:0	0:0	0:0	0:0	1:1	1:1	2:1	1:1	1:1	2:1	0:0	0:0	1:1	0:0	4:1	1:1	0:0	3:1	3:1	1:1	2:1	0:0	2:1	0:0	1:1	2:2	1:1	3:2	6:1	3:2	1:1	0:0	7:1	11:2	0:0	55	interchromosomal	1	1	1	chr19:41320341:41320470,	0
# 40803 (112 fields)


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
		    # distance between the ends
		    dondist[$10]=c[3]-$5;
		}
		else
		{
		    # distance between the beginings
		    dondist[$10]=c[2]-$4;
		}
		if(((dondist[$10])%3)==0)
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
	    fivepcds[$10]=($16!="." ? 1 : 0);
	    fivepcdsframe[$10]=((($16!=".")&&(cdslist[$10]!="")) ? 1 : 0);
	}
    }
}

NR==1{
    print $0, "fivepcds", "fivepcdsframe", "fivepokcdslist", "fivepallcdsok";
}

NR>=2&&ko[$1]==0{
    print $0, fivepcds[$1], fivepcdsframe[$1], (cdslist[$1]!="" ? cdslist[$1] : "na"), ((fivepcds[$1]&&fivepcdsframe[$1]&&(allcdsko[$1]!=1)) ? 1 : 0);
}
