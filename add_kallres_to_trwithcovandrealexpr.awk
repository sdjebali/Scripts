# add_kallisto_results_to_transcript_file.awk
# This is a program to add Kallisto results (effective length, estimated counts and tpm)
# to the transcript file with real TPM, effective length, expected counts (simple join but
# checking that the transcript length is consistent)

# example of usage 
# dir=~/work/duplicons/rnaseq.simul
# pgm=~/fragencode/tools/multi/Scripts/add_kallres_to_trwithcovandrealexpr.awk
# m=se75
# cd $dir/fifth.varcovacrossgn.se75.pe75/kallisto.$m
# awk -v fileRef=abundance.tsv -f $pgm 
# $dir/data/fifth.varcovacrossgn.se75.pe75/Mus_musculus.GRCm38.102.tr.cov.length.$m.efflg.expreads.tpm.tsv 
# > Mus_musculus.GRCm38.102.tr.cov.length.$m.efflg.expreads.tpm.efflength.estcount.kalltpm.tsv

# fileRef=abundance.tsv
# target_id	length	eff_length	est_counts	tpm
# ENSMUST00000177564.1	16	2.66459	0	0
# 119415 (5 fields)

# main input file = $dir/data/fifth.varcovacrossgn.se75.pe75/Mus_musculus.GRCm38.102.tr.cov.length.$m.efflg.expreads.tpm.tsv 
# #ID	COV	tr.length	tr.efflength	expected.reads	TPM
# ENSMUST00000024880.10	0.8778014572639276	2708	2634	30.8284	0.796369
# 117641 (6 fields)

# output file = Mus_musculus.GRCm38.102.tr.cov.length.$m.efflg.expreads.tpm.efflength.estcount.kalltpm.tsv
# #ID	COV	tr.length	tr.efflength	expected.reads	TPM	kall.efflg	kall.estcounts	kall.TPM
# ENSMUST00000024880.10	0.8778014572639276	2708	2634	30.8284	0.796369	2509	29	0.664867
# 117641 (9 fields)

BEGIN{
    OFS="\t"; 
    while (getline < fileRef >0)
    {
        lg[$1]=$2; 
        elg[$1]=$3; 
        ec[$1]=$4; 
        tpm[$1]=$5;
    }
} 

NR==1{
    print $0, "kall.efflg", "kall.estcounts", "kall.TPM";
} 

NR>=2&&$3==lg[$1]{
    print $0, elg[$1], ec[$1], tpm[$1];
}