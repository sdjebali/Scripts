# complete_gnfile_with_kallres.awk
# this script completes the gene tsv file with no header that has information about gene id, gene bt, 
# whether the gene is in our paralog set, age class, list of transcripts and number of transcripts
# with the expected TPM and the Kallisto TPM by summing the TPM of all the transcripts of the gene 


# example of usage
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# dir=~/work/duplicons/rnaseq.simul 
# pgm=~/fragencode/tools/multi/Scripts/complete_gnfile_with_kallres.awk
# m=se75
# cd $dir/fifth.varcovacrossgn.se75.pe75/kallisto.$m  
# awk -v fileRef=Mus_musculus.GRCm38.102.tr.cov.length.$m.efflg.expreads.tpm.efflength.estcount.kalltpm.tsv -f $pgm $dir/data/Mus_musculus.GRCm38.cdna.all.gn.id.bt.inourparalogs.agecl.tr.list.nb.tsv > Mus_musculus.GRCm38.cdna.all.gn.id.bt.inourparalogs.agecl.tr.list.nb.tpm.exp.kall.tsv
 
# fileRef=Mus_musculus.GRCm38.102.tr.cov.length.$m.efflg.expreads.tpm.efflength.estcount.kalltpm.tsv is formatted like this 
# #ID	COV	tr.length	tr.efflength	expected.reads	TPM	kall.efflg	kall.estcounts	kall.TPM
# ENSMUST00000024880.10	0.8778014572639276	2708	2634	30.8284	0.796369	2509	29	0.664867
# 117641 (9 fields)

# main input file = $dir/data/Mus_musculus.GRCm38.cdna.all.gn.id.bt.inourparalogs.agecl.tr.list.nb.tsv is formatted like this
# ENSMUSG00000112209.1	processed_pseudogene	0	NA	ENSMUST00000218038.1,	1
# 36727 (6 fields)

# output file = Mus_musculus.GRCm38.cdna.all.gn.id.bt.inourparalogs.agecl.tr.list.nb.tpm.exp.kall.tsv is formatted like this
# gnid	gnbt	inparalogs	ageclass	trlist	nbtr	expected.TPM	kallisto.TPM
# ENSMUSG00000112209.1	processed_pseudogene	0	NA	ENSMUST00000218038.1,	1	5.12471	4.79354
# 36728 (8 fields)


BEGIN{
    OFS="\t"; 
    while (getline < fileRef >0)
    {
        tpmexp[$1]=$6; 
        tpmkall[$1]=$9;
    } 
    print "gnid", "gnbt", "inparalogs", "ageclass", "trlist", "nbtr", "expected.TPM", "kallisto.TPM";
} 

{
    ts=0; 
    tk=0; 
    split($5,a,","); 
    k=1; 
    while(a[k]!="")
    {
        ts+=tpmexp[a[k]]; 
        tk+=tpmkall[a[k]]; 
        k++;
    } 
    print $0, ts, tk;
}

