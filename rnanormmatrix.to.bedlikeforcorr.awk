# rnanormmatrix.to.bedlikeforcorr.awk

# example
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp5/multi/atac.rna.correlation/allpeaks
# annot=~/fragencode/data/species/sus_scrofa/Sscrofa11.1.102/sus_scrofa.genes.in.gx.order.tsv
# samples=~/fragencode/workspace/sdjebali/geneswitch/wp5/multi/atac.rna.correlation/tsspeaks/common.samples.scid.myid.atacnormid.josid.rnanormid.tsv
# pgm=~/fragencode/tools/multi/Scripts/rnanormmatrix.to.bedlikeforcorr.awk
# rna=~/fragencode/workspace/sdjebali/geneswitch/wp5/rnaseq/tagada.results/normalization/res/normfeatures_normc_TMM.tsv
# time awk -v fileRef1=$annot -v fileRef2=$samples -f $pgm $rna > gn.tmmnorm.eachcommonsample.bedlike.tsv

# fileRef1 = $annot
# 1	ensembl	1:3782	+	ENSSSCG00000048769
# 1	ensembl	23782:40033	+	ENSSSCG00000027257
# 31908 (5 fields)

# fileRef2 = $samples
# m_f_F_70_C_E6192_913_913.4	SSC_WU_WP5_FT_70dpf_913.4:skeletalmuscle	skelmusfetuscontrol_R28	Fetus913_4_H10_S272_L004	skelmusfetuscontrol_28
# m_f_M_70_L_E6046_117_117.1	SSC_WU_WP5_FT_70dpf_117.1:skeletalmuscle	skelmusfetuslowfibre_R9	Fetus117_1_A12_S281_L004	skelmusfetuslowfibre_9
# 283 (5 fields)

# main input $rna
# "skelmusfetuscontrol_1"	"liverfetuscontrol_1"	...	"skelmuspigletlowfibre_20"	"liverpigletlowfibre_20"
# "ENSSSCG00000000002"	5.62773874468582	18.2041290058372	...	1.48877558582336	0.714765727974258
# 1 (312 fields)
# 26872 (313 fields)


# main output = gn.tmmnorm.eachcommonsample.bedlike.tsv

BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	# when reading the gene file, remember the coord of the tss of each gene
	split($3,a,":");
	chr[$5]=$1;
	if($4=="+")
	{
	    gbeg[$5]=a[1];
	    gend[$5]=a[1]+1;
	}
	else
	{
	    if($4=="-")
	    {
		gbeg[$5]=a[2]-1;
		gend[$5]=a[2];
	    }
	}
    }

    while (getline < fileRef2 >0)
    {
	# when reading the sample file, remember the samples to put in column and in which order
	i++;
	id[i]=$1;
	newid[$5]=$1;
	# print i, $5, $1 > "tmp1";
	s=(s)($1)("\t");
    }
}

# "skelmusfetuscontrol_1"	"liverfetuscontrol_1"	...	"skelmuspigletlowfibre_20"	"liverpigletlowfibre_20"
# "ENSSSCG00000000002"	5.62773874468582	18.2041290058372	...	1.48877558582336	0.714765727974258
# 1 (312 fields)
# 26872 (313 fields)
NR==1{
    # when in the header of the norm matrix, remember the column in which we can find all samples with ok ids
    # and print the samples in correct order as remembered when reading sample file above
    # be careful the header of the matrix has n-1 columns compared to the body
    gsub(/\"/,"",$0);
    for(j=1; j<=NF; j++)
    {
	# print $j, newid[$j] > "tmp2";
	idx[newid[$j]]=(j+1);
    }
    print "#chr", "gbeg", "gend", "gnid", s;
}

NR>=2{
    # when in the body, remember the gene tss coordinate and gene id and then go over the i ok samples of the sample file
    # and write the values in this order
    gsub(/\"/,"",$0);
    s=(chr[$1])"\t"(gbeg[$1])"\t"(gend[$1])"\t"$1"\t";
    for(k=1; k<=i; k++)
    {
	# if(NR==2)
	#    print id[k], idx[id[k]] > "tmp3";
	s=(s)($(idx[id[k]]))("\t");
    }
    print s;
}
