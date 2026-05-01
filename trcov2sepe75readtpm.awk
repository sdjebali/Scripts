# trcov2sepe75readtpm.awk
# this script is used to covert yasim transcript coverages (or depths) to transcript read count and transcript tpm
# in the case of se75 rnaseq and pe75 rnaseq, using the following formulas, and assuming for pe75 an avg frag length of 200
# !!! only works if the generation of as events was not done, otherwise transcript ids get some weird additions !!!
# for se75
##########
# Expected_number_of_reads(t) = coverage(t) * effective_length(t) /  read_length
# Effective_length(t) = length(t)-read_length+1
# Expected_TPM(t) = ((Expected_number_of_reads(t)*1000/length(t))/S)*1000000 
# S=sum(over transcripts ti) (Expected_number_of_reads(ti)*1000/length(ti))
# for pe75
##########
# Expected_number_of_reads(t) = coverage(t) * effective_length(t) /  avg_frag_length
# Effective_length(t) = length(t)-avg_frag_length+1
# Expected_TPM(t) = ((Expected_number_of_reads(t)*1000/length(t))/S)*1000000 
# S=sum(over transcripts ti) (Expected_number_of_reads(ti)*1000/length(ti))

# example on the genotool cluster
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/work/duplicons/rnaseq.simul/fifth.varcovacrossgn.se75.pe75
# pgm=~/fragencode/tools/multi/Scripts/trcov2sepe75readtpm.awk
# trlg=/work/project/fragencode/data/species/mus_musculus/GRCm38.ens102/Mus_musculus.GRCm38.102.chr.exons.tr.id.lg.tsv
# time awk -v fileRef=$trlg -f $pgm Mus_musculus.GRCm38.102.chr.exons.isoform_depth.tsv > Mus_musculus.GRCm38.102.chr.exons.trid_depth_readcount_tpm_se75_pe75.tsv

# fileRef=$trlg is like this
# ENSMUST00000003502 1053
# ENSMUST00000124574 482
# 117640 (2 fields)

# main input = Mus_musculus.GRCm38.102.chr.exons.isoform_depth.tsv is like this
# TRANSCRIPT_ID	DEPTH
# ENSMUST00000024880	0.8778014572639276
# 117641 (2 fields)

# main output = Mus_musculus.GRCm38.102.chr.exons.trid_depth_readcount_tpm_se75_pe75.tsv
# TRANSCRIPT_ID	DEPTH	se75.readcount	se75.tpm	pe75.readcount	pe75.tpm
# ENSMUST00000187790	0.015784933062471022	0.124175	0.0130821	0.0367	0.0122812
# 117641 (6 fields)

BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	# stores the transcript length as well as its effective length in se75 and pe75 modes (in case it is negaive or null we put to 0) which are the following resp
	# Effective_length(t) = length(t)-read_length+1
        # Effective_length(t) = length(t)-avg_frag_length+1
	lg[$1]=$2;
	seefflg[$1]=nn(lg[$1]-75+1);
	peefflg[$1]=nn(lg[$1]-200+1);
    }
}

# print the header of the output file
NR==1{
    print $1, $2, "se75.readcount", "se75.tpm", "pe75.readcount", "pe75.tpm";
}

# for each new transcript, stores its coverage/depth and then calculates its read count in se75 and in pe75 mode using the above formulas
# and also does the sesum and the pesum needed to comput the tr tpm in se75 and in pe75 modes given by the same formula below but the efflg is different
# S=sum(over transcripts ti) (Expected_number_of_reads(ti)*1000/length(ti))
# S=sum(over transcripts ti) (Expected_number_of_reads(ti)*1000/length(ti))
NR>=2{
    trcov[$1]=$2;
    sereadcount[$1]=($2*seefflg[$1])/75;
    sesum+=((sereadcount[$1]*1000)/lg[$1]);
    pereadcount[$1]=($2*peefflg[$1])/200;
    pesum+=((pereadcount[$1]*1000)/lg[$1]);
}


END{
    for(t in trcov)
    {
	print t, trcov[t], sereadcount[t], ((sereadcount[t]*1000/lg[t])/sesum)*1000000, pereadcount[t], ((pereadcount[t]*1000/lg[t])/pesum)*1000000;
    }
}

function nn(x){
    return (x>0 ? x : 0);
}
