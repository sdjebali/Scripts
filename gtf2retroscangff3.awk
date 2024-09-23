# gtf2retroscangff3.awk

# example
# pgm=/work/project/fragencode/tools/multi/Scripts/gtf2retroscangff3.awk
# awk -f $pgm tmp > tmp2

# input tmp
# 30	tagada	exon	475875	477918	0	-	.	gene_id "LOC_000000009373"; transcript_id "G22697.1"; ref_gene_id "G22697"; tmerge_tr_id "TM_000001540442"; transcript_biotype ""; ensembl_gene_id "."; feelnc_biotype "lncRNA"; feelnc_gene_biotype "lncRNA";
# 30	tagada	exon	475875	477918	0	-	.	gene_id "LOC_000000009373"; transcript_id "G22697.2"; ref_gene_id "G22697"; tmerge_tr_id "TM_000001540443"; transcript_biotype ""; ensembl_gene_id "."; feelnc_biotype "lncRNA"; feelnc_gene_biotype "lncRNA";
# 6 (12 fields)
# 1 (16 fields)
# 3 (20 fields)
# 13 (24 fields)
# 1 (46 fields)
# 5 (50 fields)

# output tmp2
# 30	tagada	gene	475875	479087	0	-	.	ID=LOC_000000009373;Note=protein_coding_gene;Name=LOC_000000009373
# 30	tagada	mRNA	475875	479087	0	-	.	ID=G22697.4;Parent=LOC_000000009373;Name=G22697.4;Index=6
# 34 (9 fields)


BEGIN{
    OFS="\t";
}

$3=="gene"{
    gnrow[$10]=$0;
}

$3=="transcript"{
    nbtr[$10]++;
    idx[$12]=nbtr[$10];
    trrow[$12]=$0;
    gnid[$12]=$10;
}

$3=="exon"{
    nex[$12]++;
    exrow[$12,nex[$12]]=$0;
}


$3=="CDS"{
    ncds[$12]++;
    cdsrow[$12,ncds[$12]]=$0;
}


END{
    for(t in ncds)
    {
	g=gnid[t];
	split(gnrow[g],a,"\t");
	split(a[9],b," ");
	split(g,c,"\"");
	if(ok[g]=="")
	{
	    print a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], "ID="c[2]";Note=protein_coding_gene;Name="c[2];
	    ok[g]=1;
	}

	split(trrow[t],a,"\t");
	split(a[9],b," ");
	split(t,d,"\"");
	print a[1], a[2], "mRNA", a[4], a[5], a[6], a[7], a[8], "ID="d[2]";Parent="c[2]";Name="d[2]";Index="idx[t];

	for(i=1; i<=ncds[t]; i++)
	{
	    split(cdsrow[t,i],a,"\t");
	    split(a[9],b," ");
	    print a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], "Parent="d[2]","d[2]"-Protein;";
	}

	for(i=1; i<=nex[t]; i++)
	{
	    split(exrow[t,i],a,"\t");
	    split(a[9],b," ");
	    print a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], "Parent="d[2];
	}
    }
} 
