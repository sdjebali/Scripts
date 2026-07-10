# fix.patientid.and.add.basic.info.awk

# From:
# - a fileRef1 2 column tsv file without header that has the old quantseq id (with _) and the new normalised matrix quantseq id
# - a fileRef2 5 column tsv file with header that has basic anthropomorphic information about the initial 622 patients (with -)
# - an input file of transposed normalised gene expression with samples as rows and with new ids
# Returns:
# - an output file which is very similar to the input file but that has quantseq sample id like in clinical data (with - as labExpId)
#   and with not only gene expression but also after the labExpId and before the genes, the cidno, the gender, the age, the center
#   the subjid and the diet

# example
# basedir=~/bridge/workspace/sdjebali/papers/diogenes_paper_2026
# cd $basedir/results/gathering
# pgm=$basedir/scripts/long/fix.patientid.and.add.basic.info.awk
# time zcat $basedir/results/long/normfeatures_normc_loess_keepgnless53zeros_transposed.tsv,gz | 
# awk -v fileRef1=$basedir/data/quantseq.wo27patients.cid23.lid.init.forde.tsv
# -v fileRef2=$basedir/data/622.subjects.metadata.tsv -f $pgm 
# > normfeatures_normc_loess_keepgnless53zeros_transposed_withlid.cid.gender.age.center.subjid.diet.tsv

# fileRef1=$basedir/data/quantseq.wo27patients.cid23.lid.init.forde.tsv
# 10052_2	cid2_rep1
# 10071_2	cid2_rep2
# 530 (2 fields)

# fileRef2=$basedir/data/622.subjects.metadata.tsv
# labExpId	gender	age	diet	center
# 1-001-2	2	4	3	1
# 623 (5 fields)

# main input file = zcat $basedir/results/long/normfeatures_normc_loess_keepgnless53zeros_transposed.tsv,gz
# gnid	ENSG00000000003	ENSG00000000005  ...  ENSG00000283632	ENSG00000283633
# cid2_rep1	3.710517e+01	1.597701e+02  ...  4.275128e+00	2.036942e+00
# 531 (13683 fields)

# main output file = normfeatures_normc_loess_keepgnless53zeros_transposed_withlid.cid.gender.age.center.subjid.diet.tsv
# labExpId	cidno	gender	age	center	subjid	diet	ENSG00000000003	ENSG00000000005	ENSG00000000419	ENSG00000000457	...	ENSG00000283078	ENSG00000283632	ENSG00000283633
# 1-005-2-2	2	2	4	1	1-005-2	4	3.710517e+01	1.597701e+02	1.909600e+01	5.068725e+00	...	1.051689e+00	4.275128e+00	2.036942e+00
# 531 (13689 fields)

BEGIN{
    OFS="\t"; 
    while (getline < fileRef1 >0)
    {
        lid[$2]=$1;
    } 
    while (getline < fileRef2 >0)
    {
        diet[$1]=$4; 
        age[$1]=$3;
    }
} 

NR==1{
    s="labExpId\tcidno\tgender\tage\tcenter\tsubjid\tdiet\t"; 
    for(i=2; i<NF; i++)
    {
        s=(s)($i)("\t");
    } 
    print (s)($i);
} 

NR>=2{
    split(lid[$1],a,"_"); 
    split(a[1],b,""); 
    ind=b[1]"-"b[2]b[3]b[4]"-"b[5]; 
    s=(ind"-"a[2]"\t"a[2]"\t"b[5]"\t"(age[ind]!="" ? age[ind] : "NA")"\t"b[1]"\t"ind"\t"(diet[ind]!="" ? diet[ind] : "6")"\t"); 
    for(i=2; i<NF; i++)
    {
        s=(s)($i)("\t");
    } 
    print (s)($i);
} 