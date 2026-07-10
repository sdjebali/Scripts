# fix.patientid.and.add.basic.info.forsmallexpr.awk

# From:
# - a fileRef 5 column tsv file with header that has basic anthropomorphic information about the initial 622 patients (with -)
# - an input file of transposed normalised small gene expression with samples as rows and with X first and the . as separators
# Returns:
# - an output file which is very similar to the input file but that has quantseq sample id like in clinical data (with - as labExpId)
#   and with not only gene expression but also after the labExpId and before the genes, the cidno, the gender, the age, the center
#   the subjid and the diet

# example
# cd $basedir/results/small
# pgm=$basedir/scripts/small/fix.patientid.and.add.basic.info.forsmallexpr.awk
# time awk -v fileRef=$basedir/data/622.subjects.metadata.tsv 
# -f $pgm DiOGenes_microRNA_Arrays_data_allcids_ok_transposed.tsv 
# > DiOGenes_microRNA_Arrays_data_allcids_ok_transposed_withlid.cid.gender.center.subjid.diet.tsv

# fileRef=622.subjects.metadata.tsv
# labExpId	gender	age	diet	center
# 1-001-2	2	4	3	1
# 623 (5 fields)

# main input file = DiOGenes_microRNA_Arrays_data_allcids_ok_transposed.tsv
# X1.003.1.1 11.13922114 5.192119973 11.07660404 10.71574777 10.86141329
# ...
# X9.092.1.3 10.50977795 5.150112885 11.22158347 11.02303167 10.86141329
# ...
# 1 (6609 fields)
# 1050 (6610 fields) *** after the X is the center, then we have the indiv, then the gender and finally the center

# main output file = DiOGenes_microRNA_Arrays_data_allcids_ok_transposed_withlid.cid.gender.center
# labExpId	cidno	gender	center	subjid	diet	hsa-let-7a-5p.MIMAT0000062_st	hsa-let-7a-3p.MIMAT0004481_st	hsa-let-7a-2-3p.MIMAT0010195_st	hsa-let-7b-5p.MIMAT0000063_st	...	gi:555853.gi555853_copy7_st	gi:555853.gi555853_copy8_st	gi:555853.gi555853_copy9_st
# 1-003-1-1	1	1	1	1-003-1	3	11.13922114	5.192119973	5.335105051	12.80316121	...	11.07660404	10.71574777	10.86141329
# 1051 (6615 fields)  


BEGIN{
    OFS="\t"; 
    while (getline < fileRef >0)
    {
        diet[$1]=$4;
    }
} 

NR==1{
    s="labExpId\tcidno\tgender\tcenter\tsubjid\tdiet\t"; 
    for(i=1; i<NF; i++)
    {
        s=(s)($i)("\t");
    }
    print (s)($i);
} 

NR>=2{
    split($1,a,"X"); 
    split(a[2],b,"."); 
    ind=b[1]"-"b[2]"-"b[3]; 
    s=(ind"-"b[4]"\t"b[4]"\t"b[3]"\t"b[1]"\t"ind"\t"(diet[ind]!="" ? diet[ind] : "6")"\t"); 
    for(i=2; i<NF; i++)
    {
        s=(s)($i)("\t");
    } 
    print (s)($i);
}