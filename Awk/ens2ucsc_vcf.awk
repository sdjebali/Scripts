# ens2ucsc_vcf.awk
# This is to change chromosome names from ensembl convention (integers) to ucsc convention (chr<integer>)

# example
# cd ~/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
# pgm=~/fragencode/tools/multi/Scripts/Awk/ens2ucsc_vcf.awk
# time bcftools view 21.vcf.gz | awk -f $pgm | bcftools view -O z -o 21.basic.vcf.gz  
# 12 minutes

BEGIN{OFS="\t"}
{
    if($1=="#CHROM")
    {
	print $1, $2, $3, $4, $5, $6, $7, $8, $9;
    }
    else
    {
	if($0!~/^#/)
	{
	    print "chr"$1, $2, $3, $4, $5, $6, $7, $8, $9;
	}
	else
	{
	    if(match($0,/(##contig=<ID=)(.*)/,m))
	    {
		print m[1]"chr"m[2];
	    }
	    else
	    {
		print $0;
	    }
	}
    }
}
