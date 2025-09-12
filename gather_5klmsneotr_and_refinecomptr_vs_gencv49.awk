# gather_5klmsneotr_and_refinecomptr_vs_gencv49.awk

# example
##########
# srun --x11 --mem=8G --pty bash
# cd ~/work/sarcomas/novel.annot/fromLucile/5kneotr.vs.gencv49
# gencv49=~/fragencode/data/species/homo_sapiens/GRCh38.gencv49/gencv49.gene.id.bt.name.tsv
# pgm=~/fragencode/tools/multi/Scripts/gather_5klmsneotr_and_refinecomptr_vs_gencv49.awk
# time awk -v fileRef1=LMS_neoTranscriptNotInNormal_u_lociFilter_unstrover_gencv49_ok_exons_comp_refinedclass_nbex_intermclass.tsv -v fileRef2=$gencv49 -f $pgm ../LMS_neoTranscriptNotInNormal_u_lociFilter_exons.gff > LMS_neoTranscriptNotInNormal_u_lociFilter_gnid_trlist_id_class1_clas2_gencv49tridlist_gencv49gnidlist_theirbt_theirname.tsv
# very fast

# inputs
########
# fileRef1=LMS_neoTranscriptNotInNormal_u_lociFilter_unstrover_gencv49_ok_exons_comp_refinedclass_nbex_intermclass.tsv# trid	comptrclass	annottrlist	refinedtrclass	nbex	interm_class	interm_gnlist
# tTCONS_00391991	Exact	ENST00000762149.1,	annot	n2	extension	ENSG00000299274.1,
# 4056 (7 fields)  *** we have a gene id when $(NF-1) is alternative, extension or known, and no gene id when it is novel

# fileRef2=$gencv49
# ENSG00000000003.17	protein_coding	TSPAN6
# 78899 (3 fields)

# main input: ../LMS_neoTranscriptNotInNormal_u_lociFilter_exons.gff
# chr2	scallop	exon	100004088	100004186	.	.	.	gene_id "XLOC_049628"; transcript_id "TCONS_00785818";
# 12597 (12 fields)

# output
########
# LMS_neoTranscriptNotInNormal_u_lociFilter_gnid_trlist_id_class1_clas2_gencv49tridlist_gencv49gnidlist_theirbt_theirname.tsv
# XLOC_024440	TCONS_00400211,	NA,	NA,	NA,	NA,	NA,	NA,
# XLOC_024441	TCONS_00384086,TCONS_00384087,TCONS_00384088,TCONS_00389100,TCONS_00384086,TCONS_00384087,TCONS_00384088,TCONS_00389100,	NA,NA,NA,NA,NA,NA,NA,NA,	NA,NA,NA,NA,NA,NA,NA,NA,	NA,NA,NA,NA,NA,NA,NA,NA,	NA,NA,NA,NA,NA,NA,NA,NA,	NA,NA,NA,NA,NA,NA,NA,NA,	NA,NA,NA,NA,NA,NA,NA,NA,
# ...
# XLOC_099509	TCONS_01540452,TCONS_01540452,	Overlap,Overlap,	novel,novel,	ENST00000812034.1:ENST00000416956.1,ENST00000812034.1:ENST00000416956.1,	.,.,	.,.,	.,.,
# 3132 (8 fields)


BEGIN{
    # reading refine comptr output of 5k neotr that over gencv49 with gencv49
    # and storing neotr class1, neotr class2, antrlist with : as separator when several
    # angnlist with : as separator when several
    while (getline < fileRef1 >0)
    {
	split($1,a,"t");
	c1[a[2]]=$2;
	c2[a[2]]=$6;
	$3=substr($3,1,(length($3)-1));
	gsub(/\,/,":",$3);
	antrlist[a[2]]=$3;
	# only if the gencv47 gene id exists do we remove the end separator
	if($7!=".")
	{
	    $7=substr($7,1,(length($7)-1));
	    gsub(/\,/,":",$7);
	    angnlist[a[2]]=$7;
	}
	angnlist[a[2]]=$7;
    }
    
    # reading gencv49 tsv file with gene info and storing biotype and name for each gene
    while (getline < fileRef2 >0)
    {
	gnbt[$1]=$2;
	gnname[$1]=$3;
    }
    # for gencv49 genes that are . like for novel class2 neotr we put . for bt and name
    gnbt["."]=".";
    gnname["."]=".";
}

# reading neotr exon file and storing trid list for each gene (comma separated)
{
    split($10,a,"\"");
    split($12,b,"\"");
    trlist[a[2]]=(trlist[a[2]])(b[2])(",");
}

# at the end, for each neotr gene, displays all the info we want
# - gene id
# - trid list (comma separated)
# - tr class1 list (comma separated) and inside : separation)
# - tr clas2 list (comma separated)
# - gencv49 tr id list (comma separated and inside : separation)
# - gencv49 gene list (comma separated and inside : separation)
# - same for their biotypes
# - same for their names  (but those last 3 do not work very well specially when we have several gencv49 tr associated to a neotr like below)
# XLOC_099509	TCONS_01540452,TCONS_01540452,	Overlap,Overlap,	novel,novel,	ENST00000812034.1:ENST00000416956.1,ENST00000812034.1:ENST00000416956.1,	NA,NA,	NA,NA,	NA,NA,
END{
    OFS="\t";
    for(g in trlist)
    {
	s1="";
	s2="";
	s3="";
	s4="";
	s5="";
	s6="";
	split(trlist[g],a,",");
	k=1;
	# a[k] is the id of a given neotr
	while(a[k]!="")
	{
	    s1=(s1)(na(c1[a[k]]))(",");
	    s2=(s2)(na(c2[a[k]]))(",");
	    s3=(s3)(na(antrlist[a[k]]))(",");
	    split(angnlist[a[k]],b,":");
	    l=1;
	    s51="";
	    s61="";
	    while(b[l]!="")
	    {
		s51=(s51)(gnbt[b[l]])(":");
		s61=(s61)(gnname[b[l]])(":");
		l++;
	    }
	    s4=(s4)(na(angnlist[a[k]]))(",");
	    s5=(s5)(na(substr(s51,1,(length(s51)-1))))(",");
	    s6=(s6)(na(substr(s61,1,(length(s61)-1))))(",");
	    k++;
	}
	print g, trlist[g], na(s1), na(s2), na(s3), na(s4), na(s5), na(s6);
    }
}

function na(x){
    return (x=="" ? "NA" : x);
}
