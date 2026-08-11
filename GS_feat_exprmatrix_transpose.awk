# GS_feat_exprmatrix_transpose.awk
# Given:
# - a feature matrix (small or long genes) with normalised expression in 84 GS experiments (rows=feat, col=sample from 84, with header
#   that have the expt like stage1_cerebellum_rep1),
# - a 1 column txt file with the experiments in the order we want them in the output file (as rows) and with
#   tissue, stage number and biorep number (like cerebellum_1_1)
# transpose the matrix and at the same time put the expts in the wanted order and with smaller provided names
# so we should get a file with
# - as many rows as experiments + 1 (with the example below 85)
# - as many columns as nb small feat + 5 (with the example below 367)

# Example of usage
##################
# basedir=~/fragencode/workspace/geneswitch/results/srnaseq
# pgm=~/fragencode/tools/multi/Scripts/GS_feat_exprmatrix_transpose.awk
# sp=gallus_gallus
# cd $basedir/$sp/mixomics
# time awk -v fileRef=expt.txt -f $pgm maturemirna_tmm_84expts_1tissstage_4biorepnonzero.tsv > maturemirna_tmm_84expts_1tissstage_4biorepnonzero_transposed.tsv
# real	0m0.043s

# fileRef=expt.txt
# cerebellum_1_1
# 84 (1 fields) 

# input file = maturemirna_tmm_84expts_1tissstage_4biorepnonzero.tsv
# gnid	stage1_cerebellum_rep1	stage1_cerebellum_rep2	...	stage3_skin_rep3	stage3_skin_rep4
# MIRNA_hairpin_355	9.1309204885219	4.56990884532808	...	5.33139993759284	8.3973305403451
# 363 (85 fields)

# output file = maturemirna_tmm_84expts_1tissstage_4biorepnonzero_transposed.tsv
# labExpId	tissue	dvtstage	rep	subjid	MIRNA_hairpin_355	MIRNA_hairpin_356	...	MIRNA_hairpin_264	MIRNA_hairpin_265
# cerebellum_1_1	cerebellum	s1	r1	s1_r1	9.1309204885219	0	...	6.13921957760688	348.481216869675
# 85 (367 fields)

BEGIN{
    OFS="\t";
    # read the expt to get the wanted order and also the small n ames
    # cerebellum_1_1
    while (getline < fileRef >0)
    {
	i++;
	expt[i]=$1;
	split($1,a,"_");
	oldid="stage"a[2]"_"a[1]"_rep"a[3];
	corr[oldid]=$1;
	corr[$1]=oldid;
    }
}

# when reading the header of the main input file, remember the new names of the expts in the order they appear
# gnid	stage1_cerebellum_rep1	stage1_cerebellum_rep2	...	stage3_skin_rep3	stage3_skin_rep4
NR==1{
    for(i=2; i<=NF; i++)
    {
	exptnewid[i]=corr[$i];
    }
}

# when reading the body of the main input file, remember the ids of each small as well as its expr in each expts
# using expt new names
# MIRNA_hairpin_355	9.1309204885219	4.56990884532808	...	5.33139993759284	8.3973305403451
NR>=2{
    nbsmall++;
    smallid[nbsmall]=$1;
    for(i=2; i<=NF; i++)
    {
	expr[$1,exptnewid[i]]=$i;
    }
}

# after having read the main input file, make the main output file
# first the header with the small gene ids and then for each expt on a different row, the expr of the each small gene in this expt
END{
    # 1st row, the header
    s="labExpId\ttissue\tdvtstage\trep\tsubjid\t";
    for(k=1; k<nbsmall; k++)
    {
	s=(s)(smallid[k])("\t");
    }
    print (s)(smallid[k]);
    # next rows, each expt on a different row
    for(l=1; l<=84; l++)
    {
	e=expt[l];
	split(e,a,"_");
	s=e"\t"a[1]"\ts"a[2]"\tr"a[3]"\ts"a[2]"_r"a[3]"\t";
	for(k=1; k<nbsmall; k++)
	{
	    s=(s)(expr[smallid[k],e])("\t");
	}
	print (s)(expr[smallid[k],e]);
    }
}
