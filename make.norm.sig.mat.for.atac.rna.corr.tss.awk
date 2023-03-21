# make.norm.sig.mat.for.atac.rna.corr.tss.awk
# given:
########
# - a tsv file of common samples between atacseq and rnaseq w/o header with 5 columns that has the wanted ids first
#   then an intermediate id from me, then the id it has in the header of the atacseq file/matrix, then an intermediate
#   rnaseq id from Jos and tagada pipeline and then the id it has in the header of the rna matrix (fileRef1)
# - a tsv file with header and n-1 columns in the header representing the atacseq matrix (normalized signal) with the
#   peak ids and the sample ids with double quotes but the normalized values w/o double quotes (fileRef2)
# - the same for the rnaseq matrix (fileRef3)
# - an input file that is a tsv file w/o header with 2 columns, the peak ids and the gene ids for the peaks at the tss
#   of genes and with a 1 to 1 correspondance between the peaks and the genes
# outputs two files that are tsv files with headers and that are named always the same:
#######################################################################################
# - tsspeaks1to1.loessnorm.eachcommonsample.tsv
# - gnoftsspeaks1to1.tmmnorm.eachcommonsample.tsv
# those files are the subsets of fileRef2 and fileRef3 but that only have the peaks and the genes of the input file
# and with only the experiments in common between atac and rna as indicated in fileRef1 but with the ids in column1
# of fileRef1. The two files will be in the order given in the input file. This is to be able to apply the compute
# corr python script by CM with two files as input (in same order)
# Note: here we discard the peak, gene pairs for which there is no normalized values in one of the two input matrices (atac or rna)


# example
#########
# !!! finally run with sbatch to be sure of the ressources !!!
# srun --mem=32G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp5/multi/atac.rna.correlation/tsspeaks
# atac=~/geneswitch/results/atacseq/sus_scrofa/wp5/complete/nf-core.atacseq.v1.2.1.Sscrofa11.1.102.22-03-07/normalization/res/normfeatures_normc_loess.tsv
# rna=~/fragencode/workspace/sdjebali/geneswitch/wp5/rnaseq/tagada.results/normalization/res/normfeatures_normc_TMM.tsv
# pgm=~/fragencode/tools/multi/Scripts/make.norm.sig.mat.for.atac.rna.corr.tss.awk
# time awk -v fileRef1=common.samples.scid.myid.atacnormid.josid.rnanormid.tsv -v fileRef2=$atac -v fileRef3=$rna -f $pgm promtss.peak.id.geneid.1to1.tsv

# common.samples.scid.myid.atacnormid.josid.rnanormid.tsv (fileRef1)
# m_f_F_70_C_E6192_913_913.4	SSC_WU_WP5_FT_70dpf_913.4:skeletalmuscle	skelmusfetuscontrol_R28	Fetus913_4_H10_S272_L004	skelmusfetuscontrol_28
# m_f_M_70_L_E6046_117_117.1	SSC_WU_WP5_FT_70dpf_117.1:skeletalmuscle	skelmusfetuslowfibre_R9	Fetus117_1_A12_S281_L004	skelmusfetuslowfibre_9
# 283 (5 fields)

# $atac (fileRef2)
# "skelmusfetuscontrol_R28"	"skelmusfetuslowfibre_R9"	...	"skelmusfetuslowfibre_R7"	"skelmuspigletlowfibre_R10"
# "1:21547:24249"	44.4458579633872	41.8624506001273	...	38.8445857723696	33.4147760411802
# 1 (290 fields)
# 245471 (291 fields)

# $rna (fileRef3)
# "skelmusfetuscontrol_1"	"liverfetuscontrol_1"	"...	"skelmuspigletlowfibre_20"	"liverpigletlowfibre_20"
# "ENSSSCG00000000002"	5.62773874468582	18.2041290058372	...	1.48877558582336	0.714765727974258
# 1 (312 fields)
# 26872 (313 fields)

# promtss.peak.id.geneid.1to1.tsv (main input)
# 16:60331499:60333062	ENSSSCG00000017022
# 15:45189478:45191872	ENSSSCG00000015780
# 10121 (2 fields)


BEGIN{
    OFS="\t";
    # reads the common sample file and record in which order we want the samples in the final files and how they are called
    # records the correspondance between the id we want and the other two ids
    # m_f_F_70_C_E6192_913_913.4	SSC_WU_WP5_FT_70dpf_913.4:skeletalmuscle	skelmusfetuscontrol_R28	Fetus913_4_H10_S272_L004	skelmusfetuscontrol_28
    while (getline < fileRef1 >0)
    {
	i++;
	samp[i]=$1;
	corrmyid1[$3]=$1;
	corrmyid2[$5]=$1;
    }

    # record the list of samples to be in the headers of the two wanted files
    for(k=1; k<i; k++)
    {
	header=(header)(samp[k])("\t");
    }
    header=(header)(samp[k]);
    # write the headers of the two files
    print "peakid\t"(header) > "tsspeaks1to1.loessnorm.eachcommonsample.tsv";
    print "gnid\t"(header) > "gnoftsspeaks1to1.tmmnorm.eachcommonsample.tsv";

    # record the atacseq matrix = first in which order appear my ids and then the normalized value of each peak
    # in each sample named by my ids, after removing the double quotes in the sample ids and in the peak id
    # "skelmusfetuscontrol_R28"	"skelmusfetuslowfibre_R9"	...	"skelmusfetuslowfibre_R7"	"skelmuspigletlowfibre_R10"
    # "1:21547:24249"	44.4458579633872	41.8624506001273	...	38.8445857723696	33.4147760411802
    while (getline < fileRef2 >0)
    {
	n1++;
	if(n1==1)
	{
	    gsub(/\"/,"", $0);
	    for(k=1; k<=NF; k++)
	    {
		myid1[k+1]=corrmyid1[$k];
	    }
	}
	else
	{
	    gsub(/\"/,"",$1);
	    ok1[$1]=1;
	    for(k=2; k<=NF; k++)
	    {
		val1[$1,myid1[k]]=$k;
	    }
	}
    }

    # do the same with the rna matrix
    # "skelmusfetuscontrol_1"	"liverfetuscontrol_1"	"...	"skelmuspigletlowfibre_20"	"liverpigletlowfibre_20"
    # "ENSSSCG00000000002"	5.62773874468582	18.2041290058372	...	1.48877558582336	0.714765727974258
    while (getline < fileRef3 >0)
    {
	n2++;
	if(n2==1)
	{
	    gsub(/\"/,"", $0);
	    for(k=1; k<=NF; k++)
	    {
		myid2[k+1]=corrmyid2[$k];
	    }
	}
	else
	{
	    gsub(/\"/,"",$1);
	    ok2[$1]=1;
	    for(k=2; k<=NF; k++)
	    {
		val2[$1,myid2[k]]=$k;
	    }
	}
    }
}

# for each tss peak, gene pair of the input file, write in the final peak and gene matrices and in the same order
# the submatrix corresponding to the common samples and in the order of my sample seen in fileRef1
# !!! only write the peak and the gene if both of them have normalized values in the input file !!!
# 16:60331499:60333062	ENSSSCG00000017022
{
    if((ok1[$1]==1)&&(ok2[$2]==1))
    {
	s1=$1"\t";
	s2=$2"\t";
	for(k=1; k<i; k++)
	{
	    s1=(s1)(val1[$1,samp[k]])("\t");
	    s2=(s2)(val2[$2,samp[k]])("\t");
	}
	s1=(s1)(val1[$1,samp[k]]);
	s2=(s2)(val2[$2,samp[k]]);
	print s1 > "tsspeaks1to1.loessnorm.eachcommonsample.tsv";
	print s2 > "gnoftsspeaks1to1.tmmnorm.eachcommonsample.tsv";
    }
}
