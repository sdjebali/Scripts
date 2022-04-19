# add_3classes_feelnc_tr_gn_bt_to_smerge.awk
# this script aims at adding a two additional key,value pairs to a stringtie merge file after having run feelnc
# asking for 3 classes (lncRNA, noORF and mRNA) (no TUCp class here)
# It takes as input 5 files:
############################
# - the lncRNA gtf file from feelnc classification (2nd) step as fileRef1
# - the noORF gtf file from feelnc classification (2nd) step as fileRef2
# - the mRNA gtf file from feelnc classification (2nd) step as fileRef3
# - the initial reference gtf file used by feelnc to retrieve the transcripts that are protein_coding and lncrnas 
#    and not output as such by feelnc's second step as fileRef4
# - the stringtie merge file given as input to feelnc as main input file (only exon and transcript rows, gene id and
#   transcript id as two first keys in the 9th field)
# It outputs a single gtf file which is the same as the stringtie merge file except that it has
# a feelnc_transcript_biotype and a feelnc_gene_biotype 

# example of usage
# srun --mem=8G --pty bash
# cd /work/genphyse/genroc/Carine/MONOPOLY/Bovin/feelnc
# lnc=feelnc_codpot_out/intergenic/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf
# noorf=feelnc_codpot_out/intergenic/feelnc_codpot_out/candidate_lncRNA.gtf.noORF.gtf
# mrna=feelnc_codpot_out/intergenic/feelnc_codpot_out/candidate_lncRNA.gtf.mRNA.gtf
# ref=/work/genphyse/genroc/Carine/MONOPOLY/Bovin/annotation/Bos_taurus.ARS-UCD1.2.101.gtf
# pgm=~/fragencode/tools/multi/Scripts/add_3classes_feelnc_tr_gn_bt_to_smerge.awk
# outdir=~/fragencode/workspace/sdjebali/monopoly
# time awk -v fileRef1=$lnc -v fileRef2=$noorf -v fileRef3=$mrna -v fileRef4=$ref -f $pgm stringtie_merged.gtf > $outdir/smerge.with.feelnc.biotypes.gtf

# fileRef1=$lnc
# 9	ensembl	exon	100009218	100009539	.	-	.	gene_id "ENSBTAG00000050769"; transcript_id "ENSBTAT00000070195"; exon_number "1";
# 6370 (14 fields)
# 1578 (16 fields)
# 198 (18 fields)

# fileRef2=$noorf
# 4	StringTie	exon	85372298	85372625	1000	.	.	gene_id "MSTRG.12933"; transcript_id "MSTRG.12933.1"; exon_number "1";
# 8 (14 fields)

# fileRef3=$mrna
# 8	StringTie	exon	16009426	16010079	1000	-	.	gene_id "MSTRG.15482"; transcript_id "ENSBTAT00000024401"; exon_number "1"; ref_gene_id "ENSBTAG00000018340";
# 4737 (14 fields)
# 2278 (16 fields)
# 613 (18 fields)

# fileRef4=$ref 
# #!genome-build ARS-UCD1.2
# #!genome-version ARS-UCD1.2
# 1	ensembl	gene	339070	350389	.	-	.	gene_id "ENSBTAG00000006648"; gene_version "6"; gene_source "ensembl"; gene_biotype "protein_coding";
# 5 (2 fields)
# 7751 (16 fields)
# 19856 (18 fields)
# 12357 (24 fields)
# 10651 (26 fields)
# 88847 (28 fields)
# 122049 (30 fields)
# 22 (32 fields)
# 782492 (34 fields)
# 158 (36 fields)

# input file stringtie_merged.gtf
# # stringtie --merge -G /work2/genphyse/genroc/Carine/MONOPOLY/Bovin/annotation/Bos_taurus.ARS-UCD1.2.101.gtf -p 8 -o /work2/genphyse/genroc/Carine/MONOPOLY/Bovin/feelnc/stringtie_merged.gtf /work2/genphyse/genroc/Carine/MONOPOLY/Bovin/feelnc/assembly_GTF_list.txt
# # StringTie version 2.1.1
# 1	StringTie	transcript	221464	221679	1000	.	.	gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; 
# 1 (4 fields)
# 1 (10 fields)
# 30528 (12 fields)
# 407388 (14 fields)
# 131946 (16 fields)
# 317528 (18 fields)

BEGIN{
    while (getline < fileRef1 >0)
    {
	lnc1[$12]=1;
    }
    while (getline < fileRef2 >0)
    {
	noorf[$12]=1;
    }
    while (getline < fileRef3 >0)
    {
	mrna1[$12]=1;
    }
    while (getline < fileRef4 >0)
    {
	# only read transcript rows from the reference annotation file since the tr bt can be found there and less rows than exons'
	if($3=="transcript")
	{
	    found=0;
	    split($0,a,"\t");
	    split(a[9],b," ");
	    k=1;
	    while(found==0&&b[k]!="")
	    {
		if(b[k]=="transcript_biotype")
		{
		    if(b[k+1]=="protein_coding")
		    {
			mrna2[$12]=1
		    }
		    else
		    {
			if(b[k+1]=="lncRNA")
			{
			    lnc2[$12]=1;
			}
		    }
		}
		k+=2;
	    }
	}
    }
}


# the stringtie merge file is supposed to only contain exon and transcript rows and to have
# gene_id and transcript_id as two first keys in the 9th field
# here we remember all the exons and transcript rows (indexed by tr id in hashtables)
# but we will output the complete file with feelnc tr and gn bt in the end part of the script
# note that here transcripts filtered out in 1st feelnc step and that overlap with ref pcg tr
# will probably be of the other feelnc_transcript_biotype (but we will know it is likely to be mRNA)
{
    if($3=="exon")
    {
	nbex[$12]++;
	exrow[$12,nbex[$12]]=$0;
    }
    else
    {
	if($3=="transcript")
	{
	    trrow[$12]=$0;
	    trlist[$10]=(trlist[$10])($12)(",")
	    trbt[$12]=(((mrna1[$12]==1)||(mrna2[$12]==1)) ? "mRNA" : (((lnc1[$12]==1)||(lnc2[$12]==1)) ? "lncRNA" : (noorf[$12]==1 ? "noORF" : "other")));
	}
    }
}

END{
    for(t in trrow)
    {
	split(trrow[t],a,"\t");
	split(a[9],b," ");
	# look for the tr list of the gene of the current transcript t
	# and make the tr bt list out of it in order to compute the gn bt afterwards
	split(trlist[b[2]],c,",");
	k=1;
	while(c[k]!="")
	{
	    trbtlist[b[2]]=(trbtlist[b[2]])(trbt[c[k]])(",");
	    k++;
	}
	gnbt[t]=((trbtlist[b[2]]~/mRNA/) ? "mRNA" : ((trbtlist[b[2]]~/lncRNA/) ? "lncRNA" : ((trbtlist[b[2]]~/noORF/) ? "noORF" : "other")));
	print trrow[t]" feelnc_transcript_biotype \""trbt[t]"\"\; feelnc_gene_biotype \""(gnbt[t])"\"\;";
	for(i=1; i<=nbex[t]; i++)
	{
	    print exrow[t,i]" feelnc_transcript_biotype \""trbt[t]"\"\; feelnc_gene_biotype \""(gnbt[t])"\"\;";
	}
    }

}

