# atacnorm.to.oksampleforcorr.awk

# example
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp5/multi/atac.rna.correlation/allpeaks
# samples=~/fragencode/workspace/sdjebali/geneswitch/wp5/multi/atac.rna.correlation/tsspeaks/common.samples.scid.myid.atacnormid.josid.rnanormid.tsv
# pgm=~/fragencode/tools/multi/Scripts/atacnorm.to.oksampleforcorr.awk
# atac=~/fragencode/workspace/sdjebali/geneswitch/wp5/atacseq/peaks/pairwise.correlation/peaks.loess.each.sample.bedlike.tsv
# time awk -v fileRef=$samples -f $pgm $atac > allpeaks.loessnorm.eachcommonsample.bedlike.tsv 
# real	0m50.698s


# fileRef=$samples
# m_f_F_70_C_E6192_913_913.4	SSC_WU_WP5_FT_70dpf_913.4:skeletalmuscle	skelmusfetuscontrol_R28	Fetus913_4_H10_S272_L004	skelmusfetuscontrol_28
# m_f_M_70_L_E6046_117_117.1	SSC_WU_WP5_FT_70dpf_117.1:skeletalmuscle	skelmusfetuslowfibre_R9	Fetus117_1_A12_S281_L004	skelmusfetuslowfibre_9
# 283 (5 fields)

# main input $atac
# #chr	gbeg	gend	peakid	skelmusfetuscontrol_R28	skelmusfetuslowfibre_R9	...	skelmusfetuslowfibre_R7	skelmuspigletlowfibre_R10
# 1	21547	24249	1:21547:24249	44.4458579633872	41.8624506001273	...	38.8445857723696	33.4147760411802
# 245472 (294 fields)

# main output allpeaks.loessnorm.eachcommonsample.bedlike.tsv
# #chr	gbeg	gend	peakid	m_f_F_70_C_E6192_913_913.4	m_f_M_70_L_E6046_117_117.1	...	m_f_M_70_L_E6211_115_115.3	m_p_M_60_L_E6021_201_5768    
# 1	21547	24249	1:21547:24249	44.4458579633872	41.8624506001273	...	38.8445857723696	33.4147760411802	
# 245472 (287 fields)


BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	i++;
	id[i]=$1;
	newid[$3]=$1;
	s=(s)($1)("\t");
    }
}

NR==1{
    for(j=5; j<=NF; j++)
    {
	idx[newid[$j]]=j;
    }
    print $1, $2, $3, $4, s;
}

NR>=2{
    s=$1"\t"$2"\t"$3"\t"$4"\t";
    for(k=1; k<=i; k++)
    {
	s=(s)($(idx[id[k]]))("\t");
    }
    print s;
}
