# add_info_to_twocorrfile.awk  (initial attempt in put_two_corr_together.awk with two 30k * 30k corr matrices but not efficient)
# this scripts takes as input a tsv file of gene pairs (from long+mirna) with pearson and spearman corr (on log10(normexpr+0.001)
# and 4 fileRef files and outputs a 18 column file with the pairs of genes and the two correlations but also many other pieces of info
# this output tsv file has a header and all the gene pairs present as input and the following colums:
# - chr, start, end, ID, gnclass (mirna or long) and number of samples with non null expr for the first gene
# - chr, start, end, ID, gnclass (mirna or long) and number of samples with non null expr for the second gene
# - cis or trans
# - distance if cis otherwise NA
# - pearson and spearman corr
# - if the pair involves a mirna then the genie3 score, otherwise NA
# - if the pair involves a mirna then a boolean saying if it is in the top 1M genie3 scores, otherwise NA

# the 4 files provides as fileRef0, fileRef1, fileRef2 and fileRef3 are
# - the gene annotation in bed format for the gene extents
# - the list of gene ids of mirnas
# - the gene normalized expression matrix (with header and n-1 columns in header compared to body), from which the corr were computed
# - the matrix of genie3 scores
# this script also takes as input a lower bound threshold for the top 1M genie3 scores

# example
# !!! more than 32G but less than 320G but need to assess by seff !!!
# srun --x11 --mem=320G --pty bash 
# cd /work/project/fragencode/workspace/geneswitch/results/multi/rnaseq.srnaseq/encode/quantif.files/small.longpolya
# annot=/work/project/fragencode/data/species/homo_sapiens/GRCh38.gencv29/gencode.v29.annotation.genes.bed
# mirnas=/work/project/fragencode/workspace/geneswitch/results/multi/rnaseq.srnaseq/encode/quantif.files/small_RNA-seq/expr.mirna.id.txt
# pgm=/work/project/fragencode/tools/multi/Scripts/add_info_to_twocorrfile.awk
# time zcat normfeatures_normc_TMM_twocor_gnpairs.tsv.gz | awk -v fileRef0=$annot -v fileRef1=$mirnas -v fileRef2=normfeatures_normc_TMM.tsv -v fileRef3=Genie3/GENIE3_RF_seed123.tsv -v t=0.0072398478807677 -f $pgm | gzip > normfeatures_normc_TMM_twocor_gnpairs.allinfo.tsv.gz 

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

# fileRef3=Genie3/GENIE3_RF_seed123.tsv 
# 1 (31659 fields)
# 712 (31660 fields)

# main input = zcat normfeatures_normc_TMM_twocor_gnpairs.tsv.gz
# gene1	gene2	pearson.corr	spearman.corr
# ENSG00000000003.14	ENSG00000000419.12	-0.565451126211402	-0.3
# 501130312 rows

# output = normfeatures_normc_TMM_twocor_gnpairs.allinfo.tsv.gz 
# gn1.chr	gn1.beg	gn1.end	gn1.id	gn1.class	gn1.nbnnexpr	gn2.chr	gn2.beg	gn2.end	gn2.id	gn2.class	gn2.nbnnexpr	cis.trans	dist	pearson.corr	spearman.corr	genie3.score	top1Mgenie3
# chrX	100627108	100639991	ENSG00000000003.14	long	9	chr20	50934866	50958555	ENSG00000000419.12	long	9	transNA	-0.565451126211402	-0.3	NA	NA
# 501130312 rows


BEGIN{
    OFS="\t";
    # read the gene bed file to get the coordinates of each gene
    while (getline < fileRef0 >0)
    {
	chr[$4]=$1;
	gbeg[$4]=$2;
	gend[$4]=$3;
    }
    # read the mirna id file to know which id corresponds to a mirna
    while (getline < fileRef1 >0)
    {
	mirna[$1]=1;
    }
    # read the tmm expression matrix to know in how many experiments each gene has a non null expr value
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
    # read the genie3 score matrix to associate each mirna-gene pair a genie3 score
    while (getline < fileRef3 >0)
    {
	gsub(/\"/,"",$0);
	n++;
	if(n==1)
	{
	    for(i=1; i<=NF; i++)
	    {
		gn[i+1]=$i;
	    }
	}
	else
	{
	    for(i=2; i<=NF; i++)
	    {
		genie3sc[$1,gn[i]]=$i;
		genie3sc[gn[i],$1]=$i;
	    }
	}
    }
}

NR==1{
    print "gn1.chr", "gn1.beg", "gn1.end", "gn1.id", "gn1.class", "gn1.nbnnexpr", "gn2.chr", "gn2.beg", "gn2.end", "gn2.id", "gn2.class", "gn2.nbnnexpr", "cis.trans", "dist", "pearson.corr", "spearman.corr", "genie3.score", "top1Mgenie3";
}

NR>=2{
    c=(chr[$1]==chr[$2] ? "cis" : "trans");
    d=(c=="cis" ? ((gend[$1]<gbeg[$2] ? (gbeg[$2]-gend[$1]) : (gend[$2]<gbeg[$1] ? (gbeg[$1]-gend[$2]) : 0))) : "NA");
    print chr[$1], gbeg[$1], gend[$1], $1, (mirna[$1]==1 ? "mirna" : "long"), nbnn[$1], chr[$2], gbeg[$2], gend[$2], $2, (mirna[$2]==1 ? "mirna" : "long"), nbnn[$2], c, d, $3, $4, (genie3sc[$1,$2]!="" ? genie3sc[$1,$2] : "NA"), (genie3sc[$1,$2]!="" ? (genie3sc[$1,$2]>=t ? 1 : 0) : "NA");
}
