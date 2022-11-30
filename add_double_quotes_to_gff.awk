# add_double_quotes_to_gff.awk

# example
# cd ~/fragencode/results/rnaseq/gallus_gallus/enriched.annotation.rennes
# pgm=~/fragencode/tools/multi/Scripts/add_double_quotes_to_gff.awk
# time awk -f $pgm 1_LNCextendedEns101.gtf > 1_LNCextendedEns101.ok.gtf

# input
# 1	ensembl	gene	5273	10061	.	-	.	gene_id ENSGALG00000054818; gene_version 1; gene_source ensembl; gene_biotype protein_coding; gene_remappingGg5toGg6Info EnsV100replacing_mixed;
# 1	ensembl	transcript	5273	10061	.	-	.	gene_id ENSGALG00000054818; gene_version 1; transcript_id ENSGALT00000098984; transcript_version 1; gene_source ensembl; gene_biotype protein_coding; transcript_source ensembl; transcript_biotype protein_coding; gene_remappingGg5toGg6Info EnsV100replacing_mixed;
# 10768 (18 fields)
# 13588 (20 fields)
# 171720 (26 fields)
# 11066 (28 fields)
# 63573 (30 fields)
# 160997 (32 fields)
# 7 (34 fields)
# 570832 (36 fields)
# 38 (38 fields)

# output
# 1	ensembl	gene	5273	10061	.	-	.	gene_id "ENSGALG00000054818"; gene_version "1"; gene_source "ensembl"; gene_biotype "protein_coding"; gene_remappingGg5toGg6Info "EnsV100replacing_mixed"; 
# 1	ensembl	transcript	5273	10061	.	-	.	gene_id "ENSGALG00000054818"; gene_version "1"; transcript_id "ENSGALT00000098984"; transcript_version "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; gene_remappingGg5toGg6Info "EnsV100replacing_mixed"; 
# 10768 (18 fields)
# 13588 (20 fields)
# 171720 (26 fields)
# 11066 (28 fields)
# 63573 (30 fields)
# 160997 (32 fields)
# 7 (34 fields)
# 570832 (36 fields)
# 38 (38 fields) *** real	0m34.340s  *** seems fine



{
    split($0,a,"\t");
    s="";
    for(i=1; i<=8; i++)
    {
	s=(s)(a[i])("\t");
    }

    split(a[9],b," ");
    k=1;
    while(b[k]!="")
    {
	if(k%2==1)
	{
	    s=(s)(b[k])(" ");
	}
	else
	{
	    split(b[k],c,";");
	    s=(s)("\""c[1]"\";")(" ");
	}
	k++;
    }
    print s;
}
