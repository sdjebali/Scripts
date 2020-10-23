# vcf_to_dose.awk
# This is a scripts that
# - takes as input from a vcf file of post-qc snps, that has REF and ALT alleles in the 4th and 5th field and at least DS information (0-2)
#   in the genotype field (format indicated in the 9th FORMAT field, eg GT:ADS:DS:GP) which is a reference allelic dosage and also
#   a boolean information about whether the ALT allele is minor in the INFO field
# and
# - outputs a dosage file for plink (variant-major allelic dosage) which means a tsv file with header that has the same number of columns
#   as the initial vcf file.
#   In this file the header is almost the same as the last row of the header of the input vcf file except that:
#    * ID is replaced by SNP in the 3rd field,
#    * REF and ALT are replaced by A1 and A1 in the 4th and 5th fields
#    * each of the individual's id is preceeded by its number and a space.
#   Then the body of this tsv file is the same as the vcf file except that:
#    * the 3rd field is CHROM:POS
#    * the 4th and 5th fields are the minor and major alleles and not the REF and ALT alleles
#    * the genotype fields from the 10th field only have the variant-major allelic dosage information. This dosage information goes from 0 to 2
#      and is the same as the one in the vcf file when the variant ALT allele is the minor allele and is 2-DS otherwise

# example:
##########
# datadir=~/data/parkinson/2017/HRC_Imputations/SANGER
# resdir=~/gwas/parkinson_france
# pgm=~/tools/multi/Scripts/Awk/vcf_to_dose.awk
# i=22
# cd $resdir/$i
# vcf=$datadir/$i/$i.filtered.vcf.gz
# time bcftools view $vcf | awk -f $pgm | gzip > $i.dose.plink.gz
# real	9m29.720s

# input:
########
# ##...
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT  SPD-999_3	SPD-997_3	...
# 22	16050822	.	G	A	.	PASS	RefPanelAF=0.137496;AN=6046;AC=685;INFO=0.687915;INFOin0.5-1=1;ACin60-5986=1;notstrambig=1;minor=1;infreq=0	GT:ADS:DS:GP	0|0:0,0:0:1,0,0	1|0:0.65,0.15:0.8:0.2975,0.605,0.0975	0|0:0.05,0.05:0.1:0.9025,0.095,0.0025	... 

# output:
#########
# CHROM	POS	SNP	A1	A2	QUAL	FILTER	INFO	FORMAT	1 SPD-999_3	2 SPD-997_3	...
# 22	16050822	22:16050822	A	G	.	PASS	RefPanelAF=0.137496;AN=6046;AC=685;INFO=0.687915;INFOin0.5-1=1;ACin60-5986=1;notstrambig=1;minor=1;infreq=0	DS	0	0.8	0.1 ...

BEGIN{
    OFS="\t";
}

# writes the header of the tsv file from the last row of the input vcf file header
$1=="#CHROM"{
    split($1,a,"#");
    str=a[2]"\t"$2"\tSNP\tA1\tA2\t"$6"\t"$7"\t"$8"\t"$9"\t";
    for(i=10; i<=NF-1; i++)
    {
	str=(str)((i-10+1)" "$i)("\t");
    }
    print (str)((i-10+1)" "$i);
}

# writes the body of the tsv file from the body of the input vcf file
$1!~/\#/{
    # change the 3rd field so that it is CHROM:POS
    $3=$1":"$2;
    # check whether the ALT allele is minor (usually the case), in which we case A1 and A2 will become ALT and REF but the dosage will stay the same
    # otherwise A1 and A2 will become REF and ALT and the dosage will be 2-dosage in the file
    altisminor=1;
    split($8,a,";");
    k=1;
    while(a[k]!="")
    {
	split(a[k],a1,"=");
	if(a1[1]=="minor")
	{
	    if(a1[2]==0)
	    {
		altisminor=0;
	    }
	}
	k++;
    }
    
    # Check where in the format field the dosage information is lying for the current variant/SNP and remember this as kdose
    kdose=0;
    split($9,a,":")
    k=1;
    while(a[k]!="")
    {
	if(a[k]=="DS")
	{
	    kdose=k;
	}
	k++;
    }
    
    # Makes the string of the first 9 fields
    split($9,a,":");
    str=$1"\t"$2"\t"$3"\t"((altisminor==1) ? $5 : $4)"\t"((altisminor==1) ? $4 : $5)"\t"$6"\t"$7"\t"$8"\t"a[kdose]"\t";
    # Adds the dosage information, keeping it as it is if altisminor and doing 2-it otherwise (keep the last one separate not to have a tab at the end)
    for(i=10; i<=(NF-1); i++)
    {
	split($i,a,":");
	str=(str)((altisminor==1) ? a[kdose] : (2-a[kdose]))("\t");
    }
    split($i,a,":");
    print (str)((altisminor==1) ? a[kdose] : (2-a[kdose]));
}
