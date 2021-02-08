# els.gn.tsv.to.bedpe.awk
# this is to convert a tsv file without header of ccre-distal els id, gene id pairs
# as provided by bengi, into a bedpe file, provided:
# - fileRef1 has the ccre-distal els coordinates in approx bed format (els id is in column 5 instead of 4)
# - fileRef2 has the gtf or gff2 file with at least gene rows

# example
# cd ~/regenet/results/multi/homo_sapiens/hg19/BENGI/Benchmark/All-Pairs.Natural-Ratio
# cres=~/regenet/results/multi/dnase.chipseq/ccREs/homo_sapiens/hg19/V1/gm12878/gm12878_ccRE_dELSs.bed
# genes=~/fragencode/data/species/homo_sapiens/hg19.gencv19/homo_sapiens.gtf
# pgm=~/fragencode/tools/multi/Scripts/els.gn.tsv.to.bedpe.awk
# awk -v fileRef1=$cres -v fileRef2=$genes -f $pgm bengi.gm.chic.pos.tsv | sort -V -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n > bengi.gm.chic.pos.bedpe

# inputs
# chr1	927848	928157	.	EH37E0064164	distal-ELS
# chr1	974131	974485	.	EH37E0064194	distal-ELS
# 25100 (6 fields)
# ##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)
# ##provider: GENCODE
# 4 (2 fields)
# 1 (12 fields)
# 10136 (26 fields)
# 46794 (28 fields)
# 147401 (30 fields)
# 89091 (32 fields)
# 533498 (34 fields)
# 398706 (36 fields)
# 412595 (38 fields)
# 361883 (40 fields)
# 513781 (42 fields)
# 95108 (44 fields)
# 9483 (46 fields)
# 852 (48 fields)
# 116 (50 fields)
# EH37E0064164	ENSG00000131591.13
# EH37E0064164	ENSG00000187634.6
# 88245 (2 fields)

# output
# chr1	927848	928157	chr1	803451	812283	EH37E0064164:ENSG00000230368.2	0
# 88245 (8 fields)

BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	coord1[$5]=$1":"$2":"$3;
    }
    while (getline < fileRef2 >0)
    {
	if($3=="gene")
	{
	    split($10,a,"\"");
	    coord2[a[2]]=$1":"$4":"$5;
	}
    }
}

{
    c1=coord1[$1];
    c2=coord2[$2];
    split(c1,a,":");
    split(c2,b,":");
    print a[1], a[2], a[3], b[1], b[2], b[3], $1":"$2, 0;
}
