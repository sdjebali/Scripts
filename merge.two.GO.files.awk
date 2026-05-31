# merge.two.GO.files.awk
# given a list of enriched GO terms obtained by the topGO package with their specific broad function (Metabolism here
# that could only be glucose, insulin, lipid of other) and the list of genes tested in each GO term, this script provides
# a list of testes genes present in the enriched GO terms with the list of broad functions of all the enriched GO terms
# they belong

# this script therefore takes as input two files
# - a fileRef tsv file that has a list of enriched GO terms with their GO id in column 1 and their broad function in last column
#   (kind of file provided by Cervin's topGO wrapper but with broad function at the end
# - a main input tsv file that has all the possible GO terms as rows and then a column for the number of tested genes in the term
#   and a column for the comma separared list of those (as gene ids)
# and it outputs a 2 column tsv file that has all genes that are in the enriched GO terms and with all the braod functions
# of all the enriched GO terms they are in

# example of usage
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd /work/project/bridge/workspace/sdjebali/origamir/rnaseq/norm.and.diff/liver/GO.analysis
# pgm=~/fragencode/tools/multi/Scripts/merge.two.GO.files.awk
# time awk -v fileRef=GO.BP.res.with.broad.function.tsv -f $pgm mm_go_enrichment_sig_genes_per_term.tsv > gene.in.enriched.GO.BP.id.function.tsv

# fileRef=GO.BP.res.with.broad.function.tsv 
# GO.ID	Term	Annotated	Significant	Expected	pval	Enrichment	Metabolism
# GO:0048821	erythrocyte_development	47	13	1.64	2.90E-09	7.93	other
# 97 (8 fields)

# main input file = mm_go_enrichment_sig_genes_per_term.tsv
# GO_ID	nb_sig_genes	sig_genes
# GO:0000002	1	ENSMUSG00000041064
# GO:0000018	5	ENSMUSG00000002603,ENSMUSG00000030867,ENSMUSG00000034206,ENSMUSG00000035365,ENSMUSG00000046295
# 13168 (2 fields)  
# 7043 (3 fields)  *** some rows have no gene in the 3rd column in which case column 2 is 0

# output file = gene.in.enriched.GO.BP.id.function.tsv
# ENSMUSG00000031400	glucose,other
# ENSMUSG00000030342	other
# 240 (2 fields)


BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	n=split($0,a,"\t");
	okgo[$1]=1;
	fun[$1]=a[n];
    }
}

# only when the go term is enriched and when it has a non null number of tested genes in it
okgo[$1]==1&&$2>0{
    f=fun[$1];
    split($3,a,",");
    k=1;
    while(a[k]!="")
    {
	nb[a[k]]=1;
	seen[a[k],f]=1;
	k++;
    }
}


END{
    for(g in nb)
    {
	s=(seen[g,"glucose"]==1 ? "glucose," : "")(seen[g,"insulin"]==1 ? "insulin," : "")(seen[g,"lipid"]==1 ? "lipid," : "")(seen[g,"other"]==1 ? "other" : "");
	print g, s;
    }
}
