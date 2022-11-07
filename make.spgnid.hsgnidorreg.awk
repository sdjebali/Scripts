# ~/fragencode/tools/multi/Scripts/make.spgnid.hsgnidorreg.awk
# the idea here that we have best hits when projecting livestock species tr to human as well as a tsv file with for each livestock transcript that
# we managed to project to human its gene id and the corresponding human transcript id and gene id already computed. we also have a complete gene
# annotation file for the livestock species
# and we want a 2 column tsv file with all initial livestock species genes that has
# - livestock gene id
# - human gene id or region to which it project (with redundancy in case of gene since for a livestock gene we can have several tr
#   and all of them could in theory project to different human genes) (not that if a livestock gene has several tr of which some
#   project to human tr and some others not to human tr we only keep the human tr hits
# we will then remove redundancy with a sort | uniq
# but here we want all possible livestock genes even those that did not project to the human genome (this is why we have the complete annot file)

# example
# srun --x11 --mem=8G --pty bash
# p=geneswitch
# sp=sus_scrofa
# gtf=/work/project/fragencode/workspace/geneswitch/results/rnaseq/sus_scrofa/TAGADA.v2.1.0.Sscrofa11.1.102.22-06-07/annotation/novel.gtf
# base=`basename ${gtf%.gtf}`
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp3/compgx/rnaseq/$p/$sp
# pgm=~/fragencode/tools/multi/Scripts/make.spgnid.hsgnidorreg.awk
# time zcat $p.$sp.$base\_to_hg38.besthit.psl.gz | awk -v fileRef1=tmp -v fileRef2=$gtf -f $pgm | sort | uniq > $p.$sp.$base.gnid.hsgnorreg.tsv


# input zcat $p.$sp.$base\_to_hg38.besthit.psl.gz
# 2319	0	0	0	11	69	40	21432	+	ENSSSCT00000000003	2393	5	2393	chr22	50818468	46307053	46330804	52	14,21,58,81,292,19,25,121,85,84,3,10,53,195,13,120,70,4,5,4,26,7,11,2,160,79,75,26,9,10,57,110,31,102,9,33,10,7,14,13,14,12,7,14,9,5,10,19,11,41,100,9,	5,19,40,98,188,480,499,524,645,736,820,823,834,887,1083,1096,1216,1286,1290,1295,1299,1325,1332,1351,1366,1526,1608,1683,1709,1718,1728,1785,1898,1929,2036,2062,2095,2105,2112,2126,2139,2153,2165,2172,2186,2195,2200,2210,2232,2243,2284,2384,	46307053,46307068,46308149,46308318,46308399,46308732,46308773,46308822,46312140,46312225,46316008,46316014,46316024,46316083,46316278,46316292,46323189,46323260,46323265,46326435,46326442,46326470,46326481,46326492,46326494,46328687,46328766,46328953,46328999,46329013,46329357,46329426,46329536,46330046,46330148,46330157,46330233,46330246,46330261,46330320,46330336,46330351,46330376,46330408,46330455,46330468,46330554,46330602,46330621,46330650,46330692,46330795,	2268
# 181493 (22 fields)  *** column 19 contains block sizes, column 14, 16, 17, 9 the coordinates of the hit

# fileRef1 = tmp
# trid	gnid	nbex	exlg	overtrlist	overgnid	nbtr	nbexlist	exlglist	internbexlist	interexlglist	besttr	itsgn
# TM_000000025082	LOC_000000030991	185	10526	ENST00000369208,ENST00000262901,ENST00000505753,ENST00000511871,	ENSG00000112246,ENSG00000112246,ENSG00000112246,ENSG00000112246	4	12,11,2,2	8428,3999,565,1609	121,42,2,13	7683,3946,412,874	ENST00000369208	ENSG00000112246
# 181494 (13 fields)

# fileRef2 = $gtf 
# 1	tagada	exon	1	2465	0	+	.	gene_id "LOC_000000050654"; transcript_id "ENSSSCT00000066540"; ref_gene_id "ENSSSCG00000048769"; tmerge_tr_id "TM_000000000001"; transcript_biotype "lncRNA"; feelnc_biotype "lncRNA";
# 1	tagada	transcript	1	3780	0	+	.	gene_id "LOC_000000050654"; transcript_id "ENSSSCT00000066540"; contains "ileum_stage3.filtered:STRG.1.1,skin_stage3.filtered:STRG.1.1,skin_stage2.filtered:STRG.1.1,cerebellum_stage3.filtered:STRG.1.1,muscle_stage2.filtered:STRG.1.1,kidney_stage2.filtered:STRG.1.1,muscle_stage1.filtered:STRG.1.1,lung_stage1.filtered:STRG.1.1,liver_stage3.filtered:STRG.1.1,kidney_stage1.filtered:STRG.1.1,lung_stage3.filtered:STRG.1.1,skin_stage1.filtered:STRG.1.1,ileum_stage2.filtered:STRG.1.1,kidney_stage3.filtered:STRG.1.1,muscle_stage3.filtered:STRG.1.1,ileum_stage1.filtered:STRG.1.1,cerebellum_stage2.filtered:STRG.1.1,lung_stage2.filtered:STRG.1.1,liver_stage2.filtered:STRG.1.1,liver_stage1.filtered:STRG.1.1,ref:ENSSSCT00000066540,cerebellum_stage1.filtered:STRG.1.1"; contains_count "22"; 3p_dists_to_3p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; 5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; flrpm "29.28431"; longest "ileum_stage3.filtered:STRG.1.1,skin_stage3.filtered:STRG.1.1,skin_stage2.filtered:STRG.1.1,cerebellum_stage3.filtered:STRG.1.1,muscle_stage2.filtered:STRG.1.1,kidney_stage2.filtered:STRG.1.1,muscle_stage1.filtered:STRG.1.1,lung_stage1.filtered:STRG.1.1,liver_stage3.filtered:STRG.1.1,kidney_stage1.filtered:STRG.1.1,lung_stage3.filtered:STRG.1.1,skin_stage1.filtered:STRG.1.1,ileum_stage2.filtered:STRG.1.1,kidney_stage3.filtered:STRG.1.1,muscle_stage3.filtered:STRG.1.1,ileum_stage1.filtered:STRG.1.1,cerebellum_stage2.filtered:STRG.1.1,lung_stage2.filtered:STRG.1.1,liver_stage2.filtered:STRG.1.1,liver_stage1.filtered:STRG.1.1,ref:ENSSSCT00000066540,cerebellum_stage1.filtered:STRG.1.1"; longest_FL_supporters "muscle_stage1.filtered:STRG.1.1,liver_stage3.filtered:STRG.1.1,lung_stage1.filtered:STRG.1.1,kidney_stage2.filtered:STRG.1.1,cerebellum_stage3.filtered:STRG.1.1,skin_stage2.filtered:STRG.1.1,muscle_stage2.filtered:STRG.1.1,ileum_stage3.filtered:STRG.1.1,skin_stage3.filtered:STRG.1.1,cerebellum_stage1.filtered:STRG.1.1,ref:ENSSSCT00000066540,cerebellum_stage2.filtered:STRG.1.1,liver_stage2.filtered:STRG.1.1,liver_stage1.filtered:STRG.1.1,lung_stage2.filtered:STRG.1.1,muscle_stage3.filtered:STRG.1.1,ileum_stage1.filtered:STRG.1.1,kidney_stage3.filtered:STRG.1.1,ileum_stage2.filtered:STRG.1.1,skin_stage1.filtered:STRG.1.1,kidney_stage1.filtered:STRG.1.1,lung_stage3.filtered:STRG.1.1"; longest_FL_supporters_count "22"; mature_RNA_length "3128"; meta_3p_dists_to_5p "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1"; meta_5p_dists_to_5p "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"; rpm "29.28431"; spliced "1"; ref_gene_id "ENSSSCG00000048769"; tmerge_tr_id "TM_000000000001"; transcript_biotype "lncRNA"; feelnc_biotype "lncRNA";
# 51171 (12 fields)
# 931554 (14 fields)
# 141602 (16 fields)
# 592967 (18 fields)
# 30530 (20 fields)
# 87739 (40 fields)
# 48935 (42 fields)
# 50488 (44 fields)
# 10234 (46 fields)

# output = $p.$sp.$base.gnid.hsgnorreg.tsv
# LOC_000000000000	ENSG00000121898
# LOC_000000000001	ENSG00000183077
# 54378 (2 fields)  *** some human genes, some human regions, sometimes NA, fine
# one check
# awk '{print $1}' $p.$sp.$base.gnid.hsgnorreg.tsv | sort | uniq | wc -l
# 51171  *** do correspond to the total nb of initial livestock species gene ids


BEGIN{
    OFS="\t";
    # trid	gnid	nbex	exlg	overtrlist	overgnid	nbtr	nbexlist	exlglist	internbexlist	interexlglist	besttr	itsgn
    # TM_000000025082	LOC_000000030991	185	10526	ENST00000369208,ENST00000262901,ENST00000505753,ENST00000511871,	ENSG00000112246,ENSG00000112246,ENSG00000112246,ENSG00000112246	4	12,11,2,2	8428,3999,565,1609	121,42,2,13	7683,3946,412,874	ENST00000369208	ENSG00000112246
    # 181494 (13 fields)
    while(getline < fileRef1 >0)
    {
	n++;
	if(n>=2&&$13!="NA")
	{
	    nbregtr[$1]++;
	    nbreg[$2]++;
	    reg[$2,nbreg[$2]]=$13;
	}
    }

    # 1	tagada	exon	1	2465	0	+	.	gene_id "LOC_000000050654"; transcript_id "ENSSSCT00000066540"; ref_gene_id "ENSSSCG00000048769"; tmerge_tr_id "TM_000000000001"; transcript_biotype "lncRNA"; feelnc_biotype "lncRNA";
    while(getline < fileRef2 >0)
    {
	if($3=="gene")
	{
	    split($10,a,"\"");
	    ok[a[2]]=1;
	}
	else
	{
	    if($3=="transcript")
	    {
		split($10,a,"\"");
		split($12,b,"\"");
		gnid[b[2]]=a[2];
	    }
	}
    }
}

# 2319	0	0	0	11	69	40	21432	+	ENSSSCT00000000003	2393	5	2393	chr22	50818468	46307053	46330804	52	14,21,58,81,292,19,25,121,85,84,3,10,53,195,13,120,70,4,5,4,26,7,11,2,160,79,75,26,9,10,57,110,31,102,9,33,10,7,14,13,14,12,7,14,9,5,10,19,11,41,100,9,	5,19,40,98,188,480,499,524,645,736,820,823,834,887,1083,1096,1216,1286,1290,1295,1299,1325,1332,1351,1366,1526,1608,1683,1709,1718,1728,1785,1898,1929,2036,2062,2095,2105,2112,2126,2139,2153,2165,2172,2186,2195,2200,2210,2232,2243,2284,2384,	46307053,46307068,46308149,46308318,46308399,46308732,46308773,46308822,46312140,46312225,46316008,46316014,46316024,46316083,46316278,46316292,46323189,46323260,46323265,46326435,46326442,46326470,46326481,46326492,46326494,46328687,46328766,46328953,46328999,46329013,46329357,46329426,46329536,46330046,46330148,46330157,46330233,46330246,46330261,46330320,46330336,46330351,46330376,46330408,46330455,46330468,46330554,46330602,46330621,46330650,46330692,46330795,	2268
# 181493 (22 fields)  *** column 19 contains block sizes, column 14, 16, 17, 9 the coordinates of the hit
(nbregtr[$10]=="")&&(nbreg[gnid[$10]]==""){
    split($19,a,",");
    k=1;
    while(a[k]!="")
    {
	cumulprojexlg[$10]+=a[k];
	k++;
    }
    if((nbreg[gnid[$10]]=="")||(cumulprojexlg[$10]>cumulprojexlggn[gnid[$10]]))
    {
	cumulprojexlggn[gnid[$10]]=cumulprojexlg[$10];
	nbreg[gnid[$10]]=1;
	reg[gnid[$10],1]=$14":"$16":"$17":"$9;
    }
}

END{
    for(g in nbreg)
    {
	for(i=1; i<=nbreg[g]; i++)
	{
	    print g, reg[g,i];
	}
    }

    for(g in ok)
    {
	if(nbreg[g]=="")
	{
	    print g, "NA";
	}
    }
}
