# add_contrib_to_GO_file.awk

# Given a file of enriched GO terms for BP, MF or CC obtained by Cervin's wrapper of the 
# topGO package, followed by some operations on the main output file such as 
# - a different header 
# - erasing last column
# - with underscores instead of spaces
# ie with the following information:
# GO.ID Term    Annotated       Significant     Expected        pval
# that was obtained on a list of genes that were the expressed targets of the
# X best small genes of spls comp no i, adds one column with the gene names of 
# the GO term that were also tested 
# For this, the script uses 3 different files used before reading the main input file
# - 3 col tsv file wo header with gene id, gene biotype and gene name, this file 
#   just serves as a correspondance file and not to know the expressed genes
# - 2nd tsv file obtained by Cervin's wrapper of topGO that has for all go terms, 
#   not only the signif, the term id and the number and comma separated gene id 
#   list of the terms tested there

# usage
# cd $basedir/results/go.analysis/small/males
# i=1
# go=BP
# time awk -v fileRef1=../../../../data/block.spls.quant.allgn.id.bt.name.tsv
# -v fileRef2=hs_go_enrichment_male_sig_genes_per_term.tsv 
# -f ../../../../scripts/go.analysis/small/add_contrib_to_GO_file.awk 
# hs_go_enrichment_male_$go.tsv > hs_go_enrichment_male_$go\_contribtocomp_gnnames.tsv
# 

# fileRef1=../../../../data/block.spls.quant.allgn.id.bt.name.tsv is like this
# ENSG00000000003	protein_coding	TSPAN6
# 13682 (3 fields)

# fileRef2=hs_go_enrichment_male_sig_genes_per_term.tsv is like this
# GO_ID nb_sig_genes    sig_genes
# GO:0000002    1       ENSG00000143799
# 11489 (2 fields)
# 8613 (3 fields)  *** when 2 columns then 0 in 2nd column and means no tested gene in the term 

# main input file = hs_go_enrichment_male_$go.tsv  is like this
# GO.ID Term    Annotated       Significant     Expected        pval
# GO:0035019    somatic_stem_cell_population_maintenance        51      11      3.71    0.00099
# 21 (6 fields)

# main output file = hs_go_enrichment_male_$go\_contribtocomp_gnnames.tsv is like that
# GO.ID Term    Annotated       Significant     Expected        pval    testedgn.namelist
# GO:0035019    somatic_stem_cell_population_maintenance        51      11      3.71    0.00099 REST,NKAP,MED28,SOX4,MYC,YAP1,ZFP36L2,MED21,SKI,VPS72,STAT3,
# 21 (7 fields)

BEGIN{
    OFS="\t"; 
    # ENSG00000000003	protein_coding	TSPAN6
    while (getline < fileRef1 >0)
    {
        name[$1]=$3;
    } 

    # GO:0000002    1       ENSG00000143799
    while (getline < fileRef2 >0)
    {
        split($0,a,"\t");
        if(a[2]>0)
        {
            split(a[3],b,","); 
            k=1; 
            while(b[k]!="")
            {
                seen[a[1],name[b[k]]]++;
                if(seen[a[1],name[b[k]]]==1)
                {
                    nbgn[a[1]]++; 
                    gn[a[1],nbgn[a[1]]]=name[b[k]];
                }
                k++;
            }
        }
    }
} 

{
    if(NR==1)
    {
        print $0, "testedgn.namelist";
    } 
    else
    {
        split($0,a,"\t");
        st=$0"\t"; 
        sgnname=""; 
        for(k=1; k<=nbgn[a[1]]; k++)
        {
            # if(ok[name[gn[a[1],k]]]==1)
            # {
            sgnname=(sgnname)(gn[a[1],k])(",");
            # }
        } 
        print (st)(nn(sgnname));
    }
} 

function nn(x){
    return (x=="" ? "NA" : x);
}
