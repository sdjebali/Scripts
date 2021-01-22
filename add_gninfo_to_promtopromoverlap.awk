# add_gninfo_to_promtopromoverlap.awk
# this script takes as input
# - the result of intersectBed -wao between a prom frag file and a tss+-500bp file
# - the corresponding gene file as fileRef
# and outputs a tsv file with no header that has
# - the prom frag coord (chr:beg:end)
# - the number of genes whose tss+-500bp overlaps the prom frag
# - the comma separated list of such genes (in this form chr1:948802:949920:ENSG00000187608.5:ISG15,chr1:955502:991496:ENSG00000188157.9:AGRN,)
# on Sept 2nd 2020 adding the size of the promoter at the end of the row
# !!! be careful the information of which tss of the gene extended by 500bp on each side actually overlaps the prom frag is lost !!!
# !!! and therefore we also loose the info of which tss+-500bp overlaps the prom fragment more, and therefore also the associated gene !!!

# example
# dir=~/regenet/workspace/sdjebali/egprediction/from3D/predictions/capturehic/jung.ren.2019
# ct=HCmerge
# cd $dir/$ct/score2/pall
# prom=/work/project/fragencode/data/species/homo_sapiens/hg19.gencv19/homo_sapiens.pcg.exons.capped_sites.nr.ext500bpeachside.gff
# genes=/work/project/fragencode/data/species/homo_sapiens/hg19.gencv19/pcggn.gff
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/add_gninfo_to_promtopromoverlap.awk
# module load bioinfo/bedtools2-2.29.0
# time zcat $ct.pall.score2.prom.uniq.bed.gz | intersectBed -a stdin -b $prom -wao | awk -v fileRef=$genes -f $pgm > $ct.pall.score2.prom.uniq.over.pcgtssext500eachside.gninfo.tsv
# real	0m0.882s

# inputs
# zcat $ct.pall.score2.prom.uniq.bed.gz
# chr1	1183838	1209657

# $prom
# chr1	.	TSSext500	10002326	10003326	.	-	.	gene_id	"ENSG00000162441.7";	trlist	"ENST00000400903.2,";
# 130069 (12 fields)

# $genes
# chr1	HAVANA	gene	69091	70008	.	+	.	gene_id "ENSG00000186092.4"; transcript_id "ENSG00000186092.4"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5"; level 2; havana_gene "OTTHUMG00000001094.1";
# chr1	ENSEMBL	gene	134901	139379	.	-	.	gene_id "ENSG00000237683.5"; transcript_id "ENSG00000237683.5"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1"; level 3;
# 539 (26 fields)
# 17527 (28 fields)
# 2279 (30 fields)

# intermediate input
# chr1	1183838	1209657	chr1	.	TSSext500	1191215	1192215	.	-	.	gene_id"ENSG00000160087.16";	trlist	"ENST00000467339.1,";	1001
# for some reason the geneid is next to the gene_id string

# output
# chr1:927395:936954	1	chr1:934341:935552:ENSG00000188290.6:HES4,
# chr1:943677:957199	2	chr1:948802:949920:ENSG00000187608.5:ISG15,chr1:955502:991496:ENSG00000188157.9:AGRN,
# 25590 (3 fields) 

BEGIN{
    # read the gene gff file and keep the gene name ($18) corresponding to the gene id ($10) and the gene coord in bed format
    while (getline < fileRef >0)
    {
	split($10,a,"\"");
	split($18,b,"\"");
	name[a[2]]=b[2];
	coord[a[2]]=$1":"($4-1)":"$5":"$7;
    }
}

# in case of overlap between the prom frag and the tss extended by 500 bp on each side
# remember the prom frag and associate its gene list to it
$NF!=0{
    split($12,a,"\"");
    seen[$1":"$2":"$3]++;
    if(seen[$1":"$2":"$3]==1)
    {
	i++;
	prom[i]=$1":"$2":"$3;
    }
    ok[$1":"$2":"$3,a[2]]++;
    if(ok[$1":"$2":"$3,a[2]]==1)
    {
	nbgn[$1":"$2":"$3]++;
	gnlist[$1":"$2":"$3]=(gnlist[$1":"$2":"$3])(coord[a[2]]":"a[2]":"name[a[2]])(",");
    }
}

# in case of no overlap between the prom frag and the tss extended by 500 bp on each side, put a NA gene list
$NF==0{
    i++;
    prom[i]=$1":"$2":"$3;
    nbgn[$1":"$2":"$3]=0;
    gnlist[$1":"$2":"$3]="NA";
}

END{
    OFS="\t";
    for(j=1; j<=i; j++)
    {
	p=prom[j];
	split(p,a,":");
	print p, nbgn[p], gnlist[p], a[3]-a[2];
    }
}
