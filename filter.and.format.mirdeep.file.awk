# filter.and.format.mirdeep.file.awk
# the idea here is to take for a given small rnaseq sample that has gone in the nf core smrnaseq pipeline
# the two main files from mirdeep2 as input:
# - the *_trimmed_collapsed.bed file that has the coordinates of the small found by mirdeep2 (but in fact we only use it
#   to know the id of the sample)
# - the *_trimmed_collapsed.csv tsv file that has additional information about the prob of being a true positive
# and to output:
# - a bed file with the small with prob more than 75% of being a true positive

# example
# sp=sus_scrofa
# resdir=nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28
# dir=/work/project/fragencode/workspace/geneswitch/results/srnaseq/$sp/$resdir/mirdeep2/mirdeep
# pgm=/work/project/fragencode/tools/multi/Scripts/filter.and.format.mirdeep.file.awk
# f=/work/project/fragencode/workspace/geneswitch/results/srnaseq/sus_scrofa/nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28/mirdeep2/mirdeep/result_28_06_2021_t_16_37_24_pig_stage3_newbornday1_liver_rep1_1_trimmed_collapsed.bed
# no=1
# awk -v sp=$sp -v f=$f -f $pgm ${f%.bed}.csv | sort -k1,1 -k2,2n -k3,3n > ${f%.bed}.75pcentTP.bed


# and to make a single bed file for all samples
# cat $dir/file.nbsmall.tsv | while read f no
# do
#     awk -v sp=$sp -v f=$f -f $pgm ${f%.bed}.csv 
# done | sort -k1,1 -k2,2n -k3,3n > $dir/allsamples.75pcentTP.bed


# input $f
# /work/project/fragencode/workspace/geneswitch/results/srnaseq/sus_scrofa/nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28/mirdeep2/mirdeep/result_28_06_2021_t_16_37_24_pig_stage3_newbornday1_liver_rep1_1_trimmed_collapsed.bed
# X	8239615	8239672	novel:X_34786	1696.4	-	8239615	8239672	0,0,255
# 1	224443152	224443214	novel:1_2744	123.9	-	224443152	224443214	0,0,255
# 737 (9 fields)

# input ${f%.bed}.csv
# /work/project/fragencode/workspace/geneswitch/results/srnaseq/sus_scrofa/nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28/mirdeep2/mirdeep/result_28_06_2021_t_16_37_24_pig_stage3_newbornday1_liver_rep1_1_trimmed_collapsed.csv
# miRDeep2 score	novel miRNAs reported by miRDeep2	novel miRNAs, estimated false positives	novel miRNAs, estimated true positives	known miRNAs in species	known miRNAs in data	known miRNAs detected by miRDeep2	estimated signal-to-noise	excision gearing
# 10	15	3 +/- 2	12 +/- 2 (81 +/- 12%)	48885	10110	6597 (65%)	35	7
# 20177 (0 fields)
# 1 (5 fields)
# 2 (6 fields)
# 21 (17 fields)
# 737 (19 fields)
# 1 (27 fields)
# 1 (34 fields)
# 2 (55 fields)

# output, ${f%.bed}.75pcentTP.bed
# /work/project/fragencode/workspace/geneswitch/results/srnaseq/sus_scrofa/nf-core.smrnaseq.1.1.0.Sscrofa11.1.102.21-06-28/mirdeep2/mirdeep/result_28_06_2021_t_16_37_24_pig_stage3_newbornday1_liver_rep1_1_trimmed_collapsed.75pcentTP.bed
# 1	51383792	51383852	stage3_newbornday1_liver_rep1:known:oar-miR-30c_MIMAT0030058_Ovis_aries_miR-30c:1_1845	69	-
# 1	51412850	51412913	stage3_newbornday1_liver_rep1:known:cpi-miR-30a-5p_MIMAT0037674_Chrysemys_picta_miR-30a-5p:1_1851	85	-
# 395 (6 fields)



BEGIN{
    # here the aim is to get id of the experiment in lid by looking at the name of the input file $f
    # of the type mirdeep2/mirdeep/result_28_06_2021_t_16_37_24_pig_stage3_newbornday1_liver_rep1_1_trimmed_collapsed.bed
    OFS="\t";
    nbnov1=0;
    nbkn1=0;
    n=split(f,a,"/");
    if(sp=="sus_scrofa")
    {
	split(a[n],b,"_pig_");
    }
    else
    {
	if(sp=="gallus_gallus")
	{
	    split(a[n],b,"_chicken_");
	}
    }
    split(b[2],c,"_1_trimmed_collapsed.bed");
    lid=c[1];
}


# here we start reading the mirdeep file of the type
# mirdeep2/mirdeep/result_28_06_2021_t_16_37_24_pig_stage3_newbornday1_liver_rep1_1_trimmed_collapsed.csv
# and if we enter the section novel, then we assign 1 to nov 
# novel miRNAs predicted by miRDeep2
# provisional id	miRDeep2 score	estimated probability that the miRNA candidate is a true positive	rfam alert	total read count	mature read count	loop read count	star read count	significant randfold p-value	miRBase miRNA	example miRBase miRNA with the same seed	UCSC browser	NCBI blastn	consensus mature sequence	consensus star sequence	consensus precursor sequence	precursor coordinate
$1=="novel"{
    nov=1;
}

# in case we are in the section novel and not in the header
# then we look at the probability of being a true positive in the 3rd field
# and if above 75 then we write the coordinates in bed format with information of this prob and the novelty
# X_34786	1696.4	81 +/- 12%	-	3325	3321	0	4	yes	-	-	-	-	aaccagacucugagagcaggac	ucguugcccucucugucuuggucu	ucguugcccucucugucuuggucuggauuuucccaaaccagacucugagagcaggac	X:8239615..8239672:-
$1!="novel"&&$1!="mature"&&$2!="id"&&nov==1&&NF>=1{
    split($0,a,"\t");
    split(a[3],a1," ");
    split(a1[3],a2,"%");
    if(a1[1]>=75)
    {
	split(a[1],b,"_");
	split(a[17],c,":");
	split(c[2],c1,".");
	print c[1], c1[1], c1[3], lid":novel:NA:"a[1], a1[1]-a2[1], c[3];
    }
}


# if we enter the section mature, then we assign 0 to nov and 1 to kn
# mature miRBase miRNAs detected by miRDeep2
# tag id	miRDeep2 score	estimated probability that the miRNA is a true positive	rfam alert	total read count	mature read count	loop read count	star read count	significant randfold p-value	mature miRBase miRNA	example miRBase miRNA with the same seed	UCSC browser	NCBI blastn	consensus mature sequence	consensus star sequence	consensus precursor sequence	precursor coordinate
$1=="mature"{
    nov=0;
    kn=1;
}

# when we are in the section mature but not in the header and not in the last undetected section
# with know small, we do the same as for the novel ones for the filtering but here we write known
# and with the closest known mirna found in all species and reported by mirdeep2
# 9_21350	-20.7	0.00 +/- 0.00	-	45864	45864	0	0	no	cfa-miR-1271_MIMAT0006685_Canis_familiaris_miR-1271	-	-	-	cuuggcaccuaguaagcacuca	acaccugucuggucacaggua	acaccugucuggucacagguaucccauguaauucccagcaucuagaagaauacuuggcaccuaguaagcacuca	9:122957223..122957297:-
$1!="mature"&&$1!="#miRBase"&&$1!="miRBase"&&$2!="id"&&kn==1&&NF>=1{
    split($0,a,"\t");
    split(a[3],a1," ");
    split(a1[3],a2,"%");
    if(a1[1]>=75)
    {
	split(a[1],b,"_");
	split(a[17],c,":");
	split(c[2],c1,".");
	print c[1], c1[1], c1[3], lid":known:"a[10]":"a[1], a1[1]-a2[1], c[3];
    }
}
