# make_expected_trexpr_file_from_cdna_and_art_files.awk

# takes as input two files
# - a fa.fai file of all the transcripts with their lengths, that it reads beforehand
# - a tsv file with header of transcript coverages given to art 
# and also
# - the read length for single end mode transcript generation by art that has or will be run (here 75 on the example)
# produces as output a tsv file with header that has the same content as the main
# input file except that it also has, for SE read generation by art, information about 
# - transcript length
# - transcript effective length
# - expected number of reads
# - TPM

# Example of usage
# cd $dir/data
# awk -v fileRef=Mus_musculus.GRCm38.cdna.all.fa.fai -v readlg=75 \
# -f make_expected_trexpr_file_from_cdna_and_art_files.awk \ 
# RandomTranscriptCoverage.for.art.tsv > RandomTranscriptCoverage.for.art.length.se75.efflg.expreads.tpm.tsv

# fileRef=Mus_musculus.GRCm38.cdna.all.fa.fai is formatted like this
# ENSMUST00000177564.1	16	246	16	17
# ENSMUST00000196221.1	9	509	9	10
# 119414 (5 fields)

# main input file = RandomTranscriptCoverage.for.art.tsv is formatted like this
# #ID	COV
# ENSMUST00000000001.4	11.8791630416756
# 119415 (2 fields)

# main output file = 
# RandomTranscriptCoverage.for.art.length.se75.efflg.expreads.tpm.tsv
# is formatted like this
# #ID	COV	tr.length	tr.efflength	expected.reads	TPM
# ENSMUST00000000001.4	11.8791630416756	3262	3188	504.944	35.6711
# 119415 (6 fields)


BEGIN{
    OFS="\t"; 
    while (getline < fileRef >0)
    {
        # since the coverage file that will be read next as main input might have a transcript id without dot the length is stored for $1 and for the tr id that is before the dot
        split($1,a,".");
        lg[$1]=$2;
        lg[a[1]]=$2;
    }
} 

NR==1{
    print $0, "tr.length", "tr.efflength", "expected.reads", "TPM";
} 

NR>=2{
    n++; 
    tid[n]=$1; 
    row[n]=$0; 
    cov[$1]=$2; 
    efflg[$1]=pos(lg[$1]-readlg+1); 
    readcount[$1]=cov[$1]*(efflg[$1]/readlg); 
    s+=((readcount[$1]*1000)/lg[$1]);
} 

END{
    for(i=1; i<=n; i++)
    {
        t=tid[i]; 
        print row[i], lg[t], efflg[t], readcount[t], (((readcount[t]*1000)/lg[t])/s)*1000000;
    }
} 

function pos(x){
    return (x>0 ? x : 0)
}