# add.info.to.vcf.awk
# takes
# - as fileRef a vcf file of snps to retain (unzipped and with or without header) with at least 8 fields, meaning it should include the INFO field
# - as input a vcf file with phenotype info meaning at least 10 fields, a 9th field for FORMAT and many phenotypes from field 10 then
# and makes a join between the two meaning that it will output
# - the input file filtered for snps present in fileRef (based on $1":"$2":"$4":"$5) but adding the INFO field from the fileRef file in the middle of it

# note1: we suppose that the header of the input file is already complete with FORMAT 
# note2: does not work with gz file as fileRef
# note3: does not add INFO fields in header

# example
#########
# datadir=~/data/parkinson/2017/HRC_Imputations/SANGER
# mmdir=~mmartinez/PARKINSON/2017/HRC_Imputations/SANGER
# pgm=~/tools/multi/Scripts/Awk/add.info.to.vcf.awk
# i=22
# cd $datadir/$i
# time bcftools view $mmdir/$i.vcf.gz | awk -v fileRef=$i.basic.INFOin0.5-1.ACin60-5086.notstrambig.uniqchrompos.ok.minor.infreq.info.vcf -f $pgm | bcftools view -Oz -o $i.filtered.vcf.gz


BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	if($1!~/\#/)
	{
	    id=$1":"$2":"$4":"$5;
	    info[id]=$8;
	}
    }
}

$1~/\#/{
    print
}

$1!~/\#/{
    id=$1":"$2":"$4":"$5;
    if(info[id]!="")
    {
	$8=info[id];
	print $0;
    }
}


