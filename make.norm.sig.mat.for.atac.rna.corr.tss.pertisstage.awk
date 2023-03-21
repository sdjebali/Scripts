# make.norm.sig.mat.for.atac.rna.corr.tss.pertisstage.awk
# derived from the script make.norm.sig.mat.for.atac.rna.corr.tss.awk that was done for complete matrices and not per tissue stage
# and where we had quite complex headers compared to here for tissue, stage
# !!! contrary to make.norm.sig.mat.for.atac.rna.corr.tss.awk takes into account that genes in the norm matrix can be named ensid|name or just ensid !!!
# !!! however it also supposes that the atac matrix has chr:beg:end:+ while the main input file has chr:beg:end !!!
# given:
########
# - a tsv file of common samples between atacseq and rnaseq w/o header with 2 columns that has the wanted ids first
#   that are also the ones found in the header of the norm atac-seq matrix, and then Jos ids that are present in the 
#   header of the normalized rna matrix (fileRef1)
# - a tsv file with header and n-1 columns in the header representing the atacseq matrix (normalized signal) with the
#   peak ids and the sample ids with double quotes but the normalized values w/o double quotes (fileRef2)
# - the same for the rnaseq matrix (fileRef3)
# - an input file that is a tsv file w/o header with 2 columns, the peak ids and the gene ids for the peaks at the tss
#   of genes and with a 1 to 1 correspondance between the peaks and the genes
# outputs two files that are tsv files with headers and that are named always the same:
#######################################################################################
# - tsspeaks1to1.loessnorm.eachcommonsample.tsv
# - gnoftsspeaks1to1.tmmnorm.eachcommonsample.tsv
# those files are the subsets of fileRef2 and fileRef3 but that only have the peaks and the genes of the input file
# and with only the experiments in common between atac and rna as indicated in fileRef1 but with the ids in column1
# of fileRef1. The two files will be in the order given in the input file. This is to be able to apply the compute
# corr python script by CM with two files as input (in same order)
# Note: here we discard the (peak, gene) pairs for which there is no normalized values in one of the two input matrices (atac or rna)


# example
#########
# !!! finally run with sbatch to be sure of the ressources !!!
# srun --x11 --mem=20G --pty bash
# atacdir=~/geneswitch/workspace/schalabi/biostat/gs_wp5_wo115pig_ancien.diet.label/output/res
# rnadir=~/geneswitch/results/rnaseq/sus_scrofa/teams/DGE_p1/res
# outdir=~/fragencode/workspace/sdjebali/geneswitch/wp5/multi/atac.rna.correlation/tsspeaks
# pgm=~/fragencode/tools/multi/Scripts/make.norm.sig.mat.for.atac.rna.corr.tss.pertisstage.awk
# t=liv
# s=fet
# cd $outdir/$t\_$s
# atac=$atacdir/$t\_$s\_remove-40-outlier_motherboth_fatherboth_norm_counts_LOESS.tsv
# rna=$rnadir/$t\_$s\_norm_counts_TMM.tsv
# awk -v fileRef1=$t.$s.atac.samples.ok.scid.josid.tsv -v fileRef2=$atac -v fileRef3=$rna -f $pgm ../promtss.peak.id.geneid.v108.1to1.tsv
# real	2m23.500s 

# $t.$s.atac.samples.ok.scid.josid.tsv (fileRef1)
# l_f_M_70_C_E6192_913_913.2	l_f_Male_70_C_E6192_M913_Fetus913_2
# l_f_M_60_C_Ykn_204_204.2	l_f_Male_60_C_Ykn_M204_Fetus204_2
# 76 (2 fields)

# $atac (fileRef2)
# "l_f_M_70_C_E6192_913_913.2"	"l_f_M_60_C_Ykn_204_204.2"	...	"l_f_F_70_L_E6046_117_117.3"	"l_f_F_70_H_E6272_856_856.4"
# "1:21547:24249:+"	54.4445468558339	61.1264399326516	...	60.3780108932388	64.4199357952793
# 1 (76 fields)
# 245471 (77 fields)

# $rna (fileRef3)
# "l_f_Male_70_C_E6211_M109_Fetus109_1"	"l_f_Male_70_C_E6211_M109_Fetus109_2"	...	"l_f_FeM_70_H_E6211_M994_Fetus994_3"	"l_f_FeM_70_H_E6211_M994_Fetus994_4"
# "ENSSSCG00000000002"	18.0885001992994	16.3777516613038	...	15.0411576807724	16.9205012041786
# 1 (84 fields)
# 29698 (85 fields)

# ../promtss.peak.id.geneid.v108.1to1.tsv (main input)
# 15:45189477:45191872	ENSSSCG00000015780
# 13:139608955:139610111	ENSSSCG00000063377
# 10614 (2 fields)

# output1 tsspeaks1to1.loessnorm.eachcommonsample.tsv
# peakid	l_f_M_70_C_E6192_913_913.2	l_f_M_60_C_Ykn_204_204.2	...	l_f_F_70_L_E6046_117_117.3	l_f_F_70_H_E6272_856_856.4
# 15:45189477:45191872	14.116985041959	9.98509681708511	...	11.8491220951286	9.23694126540372
# 9986 (77 fields)

# output2 gnoftsspeaks1to1.tmmnorm.eachcommonsample.tsv
# gnid	l_f_M_70_C_E6192_913_913.2	l_f_M_60_C_Ykn_204_204.2	...	l_f_F_70_L_E6046_117_117.3	l_f_F_70_H_E6272_856_856.4
# ENSSSCG00000015780	20.48791191842	12.8463372829921	...	8.75950343672055	11.8559962217951	15.0465102105975	18.894230641017	18.7029342914079	18.391362926893
# 9986 (77 fields) 

BEGIN{
    OFS="\t";
    # reads the common sample file and record in which order we want the samples in the final files and how they are called
    # records the correspondance between the id we want and the other two ids
    # l_f_M_70_C_E6192_913_913.2	l_f_Male_70_C_E6192_M913_Fetus913_2
    while (getline < fileRef1 >0)
    {
	i++;
	samp[i]=$1;
	corrmyid1[$1]=$1;
	corrmyid2[$2]=$1;
    }

    # record the list of samples to be in the headers of the two wanted files
    for(k=1; k<i; k++)
    {
	header=(header)(samp[k])("\t");
    }
    header=(header)(samp[k]);
    # write the headers of the two files
    print "peakid\t"(header) > "tsspeaks1to1.loessnorm.eachcommonsample.tsv";
    print "gnid\t"(header) > "gnoftsspeaks1to1.tmmnorm.eachcommonsample.tsv";

    # record the atacseq matrix = first in which order appear the wanted ids in the header and then the normalized value of each peak
    # in each sample named by those ids, after removing the double quotes in the sample ids and in the peak id
    # "l_f_M_70_C_E6192_913_913.2"	"l_f_M_60_C_Ykn_204_204.2"	...	"l_f_F_70_L_E6046_117_117.3"	"l_f_F_70_H_E6272_856_856.4"
    # "1:21547:24249:+"	54.4445468558339	61.1264399326516	...	60.3780108932388	64.4199357952793
    while (getline < fileRef2 >0)
    {
	n1++;
	if(n1==1)
	{
	    gsub(/\"/,"", $0);
	    for(k=1; k<=NF; k++)
	    {
		myid1[k+1]=corrmyid1[$k];
	    }
	}
	else
	{
	    gsub(/\"/,"",$1);
	    ok1[$1]=1;
	    for(k=2; k<=NF; k++)
	    {
		val1[$1,myid1[k]]=$k;
	    }
	}
    }

    # do the same with the rna matrix
    # "l_f_Male_70_C_E6211_M109_Fetus109_1"	"l_f_Male_70_C_E6211_M109_Fetus109_2"	...	"l_f_FeM_70_H_E6211_M994_Fetus994_3"	"l_f_FeM_70_H_E6211_M994_Fetus994_4"
    # "ENSSSCG00000000002"	18.0885001992994	16.3777516613038	...	15.0411576807724	16.9205012041786
    # !!! here some genes can be named with ensid or with ensid|name !!!
    while (getline < fileRef3 >0)
    {
	n2++;
	if(n2==1)
	{
	    gsub(/\"/,"", $0);
	    for(k=1; k<=NF; k++)
	    {
		myid2[k+1]=corrmyid2[$k];
	    }
	}
	else
	{
	    gsub(/\"/,"",$1);
	    split($1,a,"|");
	    ok2[a[1]]=1;
	    for(k=2; k<=NF; k++)
	    {
		val2[a[1],myid2[k]]=$k;
	    }
	}
    }
}

# for each tss peak, gene pair of the input file, write in the final peak and gene matrices and in the same order
# the submatrix corresponding to the common samples and in the order of my sample seen in fileRef1
# !!! only write the peak and the gene if both of them have normalized values in the input file !!!
# 15:45189477:45191872	ENSSSCG00000015780
# !!! the atac matrix had a + in the peak id but we dont want any and we dont have it in input file so use peakid for the corres here !!!
{
    peakid=$1":+";
    if((ok1[peakid]==1)&&(ok2[$2]==1))
    {
	s1=$1"\t";
	s2=$2"\t";
	for(k=1; k<i; k++)
	{
	    s1=(s1)(val1[peakid,samp[k]])("\t");
	    s2=(s2)(val2[$2,samp[k]])("\t");
	}
	s1=(s1)(val1[peakid,samp[k]]);
	s2=(s2)(val2[$2,samp[k]]);
	print s1 > "tsspeaks1to1.loessnorm.eachcommonsample.tsv";
	print s2 > "gnoftsspeaks1to1.tmmnorm.eachcommonsample.tsv";
    }
}
