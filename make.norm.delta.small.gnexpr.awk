# make.norm.delta.small.gnexpr.awk

# From
# - fileRef tsv file with header of the expression of all small genes (as columns) 
#   for all samples (including cid1, as rows, with `-` delimiters) with additional 
#   indication of cidno, gender, center, subjid and diet for each sample after its id in col1
# - a main tsv file without header ids (with `_` delimiters) of the patients that have 
#   small gene expression at cid2 or cid3

# example
# basedir=~/bridge/workspace/sdjebali/papers/diogenes_paper_2026
# cd $basedir/results/small
# pgm=$basedir/scripts/small/make.norm.delta.small.gnexpr.awk
# time awk -v 
# fileRef=DiOGenes_microRNA_Arrays_data_allcids_ok_transposed_withlid.cid.gender.center.subjid.diet.tsv
# -f $pgm patients_with_mirnaexpr_cid23.txt 
# > patients_with_mirnaexpr_cid23_normgnexprdiff23.tsv
# real	0m6.429s

# fileRef=DiOGenes_microRNA_Arrays_data_allcids_ok_transposed_withlid.cid.gender.center.subjid.diet.tsv
# labExpId	cidno	gender	center	subjid	diet	hsa-let-7a-5p.MIMAT0000062_st	hsa-let-7a-3p.MIMAT0004481_st	hsa-let-7a-2-3p.MIMAT0010195_st	...	gi:555853.gi555853_copy7_st	gi:555853.gi555853_copy8_st	gi:555853.gi555853_copy9_st
# 1-003-1-1	1	1	1	1-003-1	3	11.13922114	5.192119973	5.335105051	...	11.07660404	10.71574777	10.86141329
# 1051 (6615 fields)

# main input file = patients_with_mirnaexpr_cid23.txt
# 2_001_2
# 2_002_2
# 163 (1 fields)

# main output file = patients_with_mirnaexpr_cid23_normgnexprdiff23.tsv
# patient	hsa-let-7a-5p.MIMAT0000062_st	hsa-let-7a-3p.MIMAT0004481_st	hsa-let-7a-2-3p.MIMAT0010195_st	...	gi:555853.gi555853_copy7_st	gi:555853.gi555853_copy8_st	gi:555853.gi555853_copy9_st
# 2_001_2	0.119966	0.00173605	-0.0171799	...	-0.0097605	-0.0524913	-0.0904286
# 164 (6610 fields) 

BEGIN{
    OFS="\t"; 
    while (getline < fileRef >0)
    {
        n++; 
        if(n==1)
        {
            s="patient\t"; 
            tot=NF; 
            for(i=7; i<NF; i++)
            {
                s=(s)($i)("\t");
            } 
            print (s)($i);
        }
        else {
        {
            gsub(/-/,"_",$5); 
            for(i=7; i<=NF; i++)
            {
                expr[$5,i,$2]=$i;
            }
        }
    }
} 

{
    s=$1"\t"; 
    for(i=7; i<tot; i++)
    {
        s=(s)(nn((expr[$1,i,3]-expr[$1,i,2]),expr[$1,i,2]))("\t");
    } 
    print (s)(nn((expr[$1,i,3]-expr[$1,i,2]),expr[$1,i,2]));
} 

function nn(x,y){
    return (y!=0 ? x/y : "NA");
}