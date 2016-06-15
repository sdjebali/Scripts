
# make_gff_ok.awk

# this script takes a gff version 2 file (or gtf) with at least gene_id and transcript_id information
# and outputs the same file but with only gene_id and transcript_id information, as the two first
# (key,value) pairs
# on June 15th 2016 added the other (key,value) pairs after the gene_id and transcript_id and double quotes around NA

# usage
# awk -f make_gff_ok.awk input.gff > output.gff

# example
# input
# 1	ensembl	gene	1735	16308	.	+	.	gene_id "ENSGALG00000009771"; gene_version "4"; gene_source "ensembl"; gene_biotype "protein_coding";

# output
# 1	ensembl	gene	1735	16308	.	+	.	gene_id "ENSGALG00000009771"; transcript_id "NA"; gene_version "4"; gene_source "ensembl"; gene_biotype "protein_coding";


$1!~/#/{
    s="";
    ten="\"NA\";";
    twelve="\"NA\";";
    split($0,a,"\t");
    split(a[9],b,"; ");
    k=1;
    while(b[k]!="")
    {
	split(b[k],c," ");
	if(c[1]=="gene_id")
	{
	    ten=((substr(c[2],length(c[2]))==";") ? c[2] : c[2]";");
	    gsub(/;;/,";",ten);
	}
	else
	{
	    if(c[1]=="transcript_id")
	    {
		twelve=((substr(c[2],length(c[2]))==";") ? c[2] : c[2]";");
		gsub(/;;/,";",twelve);
	    }
	    else
	    {
		nth=((substr(c[2],length(c[2]))==";") ? (c[1]" "c[2]) : (c[1]" "c[2]";"));
		gsub(/;;/,";",nth);
		s=(s)(nth)("\t");
	    }
	}
	k++
    }
    
    
    print a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]"\t"a[7]"\t"a[8]"\tgene_id "ten" transcript_id "(twelve)" "(s);
}
