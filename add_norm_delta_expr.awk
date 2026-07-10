# add_norm_delta_expr.awk

# From:
# - a transposed matrix of long gene expression for samples representing
#   patients at cid2 or cid3 with an id like 1-005-2-2 (cid at the end)
#   and starting with labExpId	cidno	gender	age	center	subjid	diet	
#   and then the long gene expression values, 
# - a 1 column file with patients ids (subjid with _ and not -) to keep
# This script produces a matrix with patients as rows and as columns for each 
# long gene, its normalized expression difference between cid2 and cid3

# example
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/work/mixomics/Viguerie.Moro.Obesity/quantseq
# pgm=~/fragencode/tools/multi/Scripts/add_norm_delta_expr.awk
# time awk -v fileRef=normfeatures_normc_loess_keepgnless53zeros_transposed_withlid.cid.gender.age.center.subjid.diet.tsv 
# -f $pgm patients_with_quantseqexpr_cid23.txt 
# > patients_with_quantseqexpr_cid23_normgnexprdiff23.tsv
# patient	ENSG00000000003	ENSG00000000005	ENSG00000000419	...	ENSG00000283078	ENSG00000283632	ENSG00000283633
# 2_001_2	0.309029	-0.407913	0.31892	...	1104.17	2.11297	-0.800016
# 144 (13683 fields)

# inputs
########
# fileRef=normfeatures_normc_loess_keepgnless53zeros_transposed_withlid.cid.gender.age.center.subjid.diet.tsv 
# labExpId	cidno	gender	age	center	subjid	diet	ENSG00000000003	ENSG00000000005	ENSG00000000419	ENSG00000000457	...	ENSG00000283632	ENSG00000283633	ENSG00000283638
# 1-005-2-2	2	2	4	1	1-005-2	4	3.710517e+01	1.597701e+02	1.909600e+01	5.068725e+00	...	4.275128e+00	2.036942e+00	0.000000e+00
# 531 (13689 fields) 

# main input = patients_with_quantseqexpr_cid23.txt 
# 2_001_2
# 2_002_2
# 143 (1 fields) 

# output = patients_with_quantseqexpr_cid23_normgnexprdiff23.tsv
################################################################
# patient	ENSG00000000003	ENSG00000000005	ENSG00000000419	...	ENSG00000283078	ENSG00000283632	ENSG00000283633
# 2_001_2	0.309029	-0.407913	0.31892	...	1104.17	2.11297	-0.800016
# 144 (13683 fields)

BEGIN{
    OFS="\t"; 
    m=1; 
    while (getline < fileRef >0)
    {
        n++; 
        if(n==1)
        {
            s="patient\t"; 
            tot=NF; 
            for(i=8; i<NF; i++)
            {
                s=(s)($i)("\t");
            } 
            print (s)($i);
        }
        else
        {
            gsub(/-/,"_",$6); 
            for(i=8; i<=NF; i++)
            {
                expr[$6,i,$2]=$i; 
                if($i!=0&&$i<m)
                {
                    m=$i;
                }
            }
        }
    } 
    ps=int(log(m/10)/log(10))-1;
} 

{
    s=$1"\t"; 
    for(i=8; i<tot; i++)
    {
        s=(s)(div((expr[$1,i,3]-expr[$1,i,2]),expr[$1,i,2],10^ps))("\t");
    } 
    print (s)(div((expr[$1,i,3]-expr[$1,i,2]),expr[$1,i,2],10^ps));
} 

function div(x,y,z){
    return x/(y+z);
}