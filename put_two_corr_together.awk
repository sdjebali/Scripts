# put_two_corr_together.awk
# !!! takes far too much ram on matrices of 30k * 30k genes but work on smaller ones !!!
# this scripts takes two squared matrices of (mirna+long) genes with pearson and spearman corr (on the log of the norm expr + 0.001)
# and makes a single tsv file with header of all gene pairs with the 16 following pieces of information
# - chr, start, end, ID, gnclass (mirna or long) and number of samples with non null expr for the first gene
# - chr, start, end, ID, gnclass (mirna or long) and number of samples with non null expr for the second gene
# - cis or trans
# - distance if cis otherwise NA
# - pearson and spearman corr
# The pearson corr matrix file is the main input while the spearman corr matrix file is given as fileRef3
# there are 3 other files provided as fileRef0, fileRef1 and fileReg2, which are
# - the gene annotation in bed format for the gene extents
# - the list of gene ids of mirnas
# - the gene normalized expression matrix (with header and n-1 columns in header compared to body), from which the corr were computed

# example
# cd /work/project/fragencode/workspace/geneswitch/results/multi/rnaseq.srnaseq/encode/quantif.files/small.longpolya
# annot=/work/project/fragencode/data/species/homo_sapiens/GRCh38.gencv29/gencode.v29.annotation.genes.bed
# mirnas=/work/project/fragencode/workspace/geneswitch/results/multi/rnaseq.srnaseq/encode/quantif.files/small_RNA-seq/expr.mirna.id.txt
# pgm=/work/project/fragencode/tools/multi/Scripts/put_two_corr_together.awk
# time zcat normfeatures_normc_TMM_cor_matrix.tsv.gz | awk -v fileRef0=$annot -v fileRef1=$mirnas -v fileRef2=normfeatures_normc_TMM.tsv -v fileRef3=normfeatures_normc_TMM_cor_spearman_matrix.tsv -f $pgm > gnpairs.withtwocorrs.tsv
# but in fact on big 30k * 30k matrices it requires too much ram to work (more than 256G) so not continued
# here result is provided when we use as tmp1 and tmp2 the head of the two big matrices and only the 10 first columns (9 for the header)


# fileRef0=$annot 
# chr1	11868	14409	ENSG00000223972.5	0	+
# chr1	14403	29570	ENSG00000227232.5	0	-
# 58721 (6 fields)

# fileRef1=$mirnas
# ENSG00000278791.1
# ENSG00000198976.1
# 712 (1 fields)

# fileRef2=normfeatures_normc_TMM.tsv
# A549	HT1080	IMR.90	MCF.7	NCI.H460	neural_cell	SK.MEL.5	SK.N.DZ	SK.N.SH
# ENSG00000000003.14	37.7754036551637	11.051934268169	15.3522723228464	24.3213576270408	29.5651606898964	29.2514870712725	44.2043849306361	0.715292274180165	32.1228538102926
# 1 (9 fields)
# 31659 (10 fields)

# fileRef3=normfeatures_normc_TMM_cor_spearman_matrix.tsv 
# if we use a much smaller file obtained like this
# awk 'NR==1{print $1, $2, $3, $4, $5, $6, $7, $8, $9} NR>=2{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' normfeatures_normc_TMM_cor_spearman_matrix.tsv | head > tmp2
# "ENSG00000000003.14" "ENSG00000000419.12" "ENSG00000000457.13" "ENSG00000000460.16" "ENSG00000000938.12" "ENSG00000000971.15" "ENSG00000001036.13" "ENSG00000001084.12" "ENSG00000001167.14"
# "ENSG00000000003.14" 1 -0.3 -0.233333333333333 -0.333333333333333 0.365148371670111 0.216666666666667 0.133333333333333 0.416666666666667 -0.333333333333333
# 1 (9 fields)
# 9 (10 fields)


# main input zcat normfeatures_normc_TMM_cor_matrix.tsv.gz 
# if we use a much smaller file obtained like this
# zcat normfeatures_normc_TMM_cor_matrix.tsv.gz | awk 'NR==1{print $1, $2, $3, $4, $5, $6, $7, $8, $9} NR>=2{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' | head > tmp1
# "ENSG00000000003.14" "ENSG00000000419.12" "ENSG00000000457.13" "ENSG00000000460.16" "ENSG00000000938.12" "ENSG00000000971.15" "ENSG00000001036.13" "ENSG00000001084.12" "ENSG00000001167.14"
# "ENSG00000000003.14" 1 -0.565451126211402 -0.104923512189837 -0.330220869061649 0.309197701134861 0.0816259805034719 0.798917270861509 0.2322140108244 -0.403388397921093
# 1 (9 fields)
# 9 (10 fields)

# main output gnpairs.withtwocorrs.tsv would be like this if we only had 9*9 matrices
# gn1.chr	gn1.beg	gn1.end	gn1.id	gn1.class	gn1.nbnnexpr	gn2.chr	gn2.beg	gn2.end	gn2.id	gn2.class	gn2.nbnnexpr	cis.trans	dist	pearson.corr	spearman.corr
# chrX	100627108	100639991	ENSG00000000003.14	long	9	chrX	100627108	100639991	ENSG00000000003.14	long	9	cis	01	1
# 82 (16 fields)


BEGIN{
    OFS="\t";
    print "gn1.chr", "gn1.beg", "gn1.end", "gn1.id", "gn1.class", "gn1.nbnnexpr", "gn2.chr", "gn2.beg", "gn2.end", "gn2.id", "gn2.class", "gn2.nbnnexpr", "cis.trans", "dist", "pearson.corr", "spearman.corr";
    while (getline < fileRef0 >0)
    {
	chr[$4]=$1;
	gbeg[$4]=$2;
	gend[$4]=$3;
    }
    while (getline < fileRef1 >0)
    {
	mirna[$1]=1;
    }
    while (getline < fileRef2 >0)
    {
	n++;
	if(n>=2)
	{
	    for(i=2; i<=NF; i++)
	    {
		if($i>0)
		{
		    nbnn[$1]++;
		}
	    }
	}
    }
    n=0;
    while (getline < fileRef3 >0)
    {
	gsub(/\"/,"",$0);
	n++;
	if(n==1)
	{
	    for(i=1; i<=NF; i++)
	    {
		gnid[i+1]=$i;
	    }
	}
	else
	{
	    for(i=2; i<=NF; i++)
	    {
		spear[$1,gnid[i]]=$i;
	    }
	}
    }
}

{
    gsub(/\"/,"",$0);
    if(NR==1)
    {
	for(i=1; i<=NF; i++)
	{
	    gnid[i+1]=$i;
	}
    }
    else
    {
	for(i=2; i<=NF; i++)
	{
	    g=gnid[i];
	    c=(chr[$1]==chr[g] ? "cis" : "trans");
	    d=(c=="cis" ? ((gend[$1]<gbeg[g] ? (gbeg[g]-gend[$1]) : (gend[g]<gbeg[$1] ? (gbeg[$1]-gend[g]) : 0))) : "NA");
	    print chr[$1], gbeg[$1], gend[$1], $1, (mirna[$1]==1 ? "mirna" : "long"), nbnn[$1], chr[g], gbeg[g], gend[g], g, (mirna[g]==1 ? "mirna" : "long"), nbnn[g], c, d, $i, spear[$1,g];
	}
    }
}
