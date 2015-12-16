
# make_gff_ok.awk

# this script takes a gff version 2 file (or gtf) with at least gene_id and transcript_id information
# and outputs the same file but with only gene_id and transcript_id information, as the two first
# (key,value) pairs

# usage
# awk -f make_gff_ok.awk input.gff > output.gff

$1!~/#/{
    ten="NA";
    twelve="NA";
    split($0,a,"\t");
    split(a[9],b,"; ");
    k=1;
    while(b[k]!="")
    {
	split(b[k],c," ");
	if(c[1]=="gene_id")
	{
	    ten=c[2]";";
	}
	else
	{
	    if(c[1]=="transcript_id")
	    {
		twelve=c[2]";";
	    }
	}
	k++
    }
    print a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]"\t"a[7]"\t"a[8]"\tgene_id "ten" transcript_id "twelve;
}