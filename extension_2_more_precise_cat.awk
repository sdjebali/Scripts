# extension_2_more_precise_cat.awk
# this script aims at knowing when a predicted tr is an extension of a ref tr, how many additional introns we have in 5' and in 3'
# with the final aim to know if the extensions occur more in 5' or in 3' (when extending an isoseq annot we expect the extensions
# to be more in 5' for example)
# !!! Note that some extension tr could have 2 exons and extend a tr with 1 exon only !!!
# !!! Note that the main input file has a t in front of the predicted transcript id

# this script takes as input:
#############################
# - a tsv file with header from the refine_comptr script with indication for each predicted transcript of its class (Exact, Extension, Inclusion, ...)
#   and the list of ref transcripts it is associated to (the list of ref tr it extends if the class is Extension)
# - two auxiliary files taken as fileRef1 and fileRef2 that are transcript with intron chain tsv files output by gff2trintsv.sh
#   on the reference transcripts and the predicted transcripts
# and outputs:
##############
# - a tsv files with header that has the following information
#   * predicted tr id
#   * corresponding ref tr list
#   * ref tr with more introns and in case of equality pick 1st
#   * strand of the locus
#   * intron chain of the predicted tr
#   * intron chain of the reference tr (monoex coord if it is a monoex)
#   * nb of introns of the predicted tr
#   * nb of introns of the ref tr
#   * list of additional 5' introns of the predicted tr tr wrt ref tr
#   * listr of additional 3' introns of the tagada tr tr wrt ref tr
#   * nb of additional 5' introns of the predicted tr wrt ref tr
#   * nb of additional 3' introns of the predicted tr wrt ref tr

# example
# cd /work/project/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/TAGADA.v2.1.1.GRCg6a.isoseq.23-01-13/control/annotation/string
# isoseq=/work/project/fragencode/workspace/geneswitch/results/isoseq/gallus_gallus/chicken_ensembl_gs.alltr.tr.gn.id.chr.str.intrnb.intrlistorexcoord.tsv
# tagiso=/work/project/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/TAGADA.v2.1.1.GRCg6a.isoseq.23-01-13/annotation/novel.alltr.tr.gn.id.chr.str.intrnb.intrlistorexcoord.tsv
# pgm=~/fragencode/tools/multi/Scripts/extension_2_more_precise_cat.awk
# time awk -v fileRef1=$isoseq -v fileRef2=$tagiso -f $pgm novel_complete_comp_refinedclass_nbex_intermclass.tsv > novel.extensions.predtrid.reftrlist.withmorex.strand.predtr.reftr.intr.list.nb.5p.3p.intrlist.nbintr.tsv
# *** real	0m4.162s  *** seems to be working, OK nb rows, ok nb columns, just when the reftr is monoex less good
 
# fileRef1=/work/project/fragencode/workspace/geneswitch/results/isoseq/gallus_gallus/chicken_ensembl_gs.alltr.tr.gn.id.chr.str.intrnb.intrlistorexcoord.tsv
# 221677 (7 fields)   *** those are tagada transcripts and column 3 provides a list of isoseq transcripts
# G30820.2	G30820	7	-	7	30341293:30342197,30342340:30346596,30346709:30359082,30359180:30375331,30375424:30376520,30376565:30416979,30417100:30417236,
# G20797.6	G20797	3	+	5	41221770:41223970,41224078:41224243,41224375:41224862,41224947:41226309,41226444:41226557,
# 139659 (6 fields)

# fileRef2=/work/project/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/TAGADA.v2.1.1.GRCg6a.isoseq.23-01-13/annotation/novel.alltr.tr.gn.id.chr.str.intrnb.intrlistorexcoord.tsv
# TM_000001878500	LOC_000000010115	5	-	3	40291851:40300308,40300448:40301200,40301344:40317547,
# G17938.1	LOC_000000060139	25	-	1	2326073:2327136,
# 221676 (6 fields)

# input=novel_complete_comp_refinedclass_nbex_intermclass.tsv 
# trid	comptrclass	annottrlist	refinedtrclass	nbex	interm_class	interm_gnlist
# tTM_000000000060	Monoexonic	.	intergenic	n1	novel	.
# tTM_000000000852	Extension	G13.1,G13.2,	extension	n37	alternative	G13,

# output=novel.extensions.predtrid.reftrlist.withmorex.strand.predtr.reftr.intr.list.nb.5p.3p.intrlist.nbintr.tsv
# predtrid	reftrlist	reftr	strand	predtr.intr.list	reftr.intr.list	predtr.intr.nb	reftr.intr.nb	5p.intr.list	3p.intr.list	5p.intr.nb	3p.intr.nb
# TM_000000000061	G2.1,	G2.1	+	21996:22761,	22474:24155	1	0	21996:22761,	NA	1	0   *** OK 
# TM_000000000209	G4.40,	G4.40	+	30503:31389,31485:32229,32332:53983,54091:54189,	54091:54189,	4	1	30503:31389,31485:32229,32332:53983,	NA	3	0  *** OK
# TM_000000000215	G4.29,	G4.29	+	30503:31389,31485:32229,32332:32539,32722:33080,33188:33286,	32332:32539,32722:33080,33188:33286,	5	3	30503:31389,31485:32229,	NA	2	0  *** OK
# TM_000000001822	G17.6,	G17.6	-	110846:112021,112118:112689,	110996:112757	2	0	NA	110846:112021,	0	1   *** OK but in case the ref tr is monoexonic it is less reliable, see below
# TM_000000001892	G5.24,G5.31,G5.29,G5.36,G5.35,	G5.24	-	194912:195000,195108:195301,195484:195691,195785:199824,	194912:195000,	4	1	195108:195301,195484:195691,195785:199824,	NA  *** OK	3	0   *** OK
# TM_000000001891	G5.24,G5.29,G5.32,	G5.24	-	194912:195000,195108:195301,195484:195691,195785:197149,197343:200013,	194912:195000,	5	1	195108:195301,195484:195691,195785:197149,197343:200013,	NA	4	0   *** OK
# 41747 (12 fields)




BEGIN{
    OFS="\t";
    # remember the intron list and the number of introns of the reference transcripts
    # !!! here monoex transcripts have 0 intron but the intron list is in fact the monoex coordinates !!!
    while (getline < fileRef1 >0)
    {
	intlist[$1]=$6;
	nbint[$1]=$5;
    }

    # remember the strand, the intron list and the number of introns of the predicted transcripts
    while (getline < fileRef2 >0)
    {
	str[$1]=$4;
	intlist[$1]=$6;
	nbint[$1]=$5;
    }
    # print the header of the file we want
    print "predtrid", "reftrlist", "reftr", "strand", "predtr.intr.list", "reftr.intr.list", "predtr.intr.nb", "reftr.intr.nb", "5p.intr.list", "3p.intr.list", "5p.intr.nb", "3p.intr.nb";
}

# only retain the predicted transcripts that are extensions of the ref transcripts
$2=="Extension"{
    # get access to the real predicted transcript id (a[2])
    split($1,a,"t");
    # get access to the reference transcript corresponding to it and with the max number of introns
    split($3,b,",");
    nbmax=-1;
    trmax="NA";
    k=1;
    while((nbmax==-1)&&(b[k]!=""))
    {
	if(nbint[b[k]]>nbmax)
	{
	    nbmax=nbint[b[k]];
	    trmax=b[k];
	}
	k++;
    }

    # get the list of 5' introns of the predicted transcripts that are not present in the ref tr
    # and the same for the 3' introns
    list5p="";
    list3p="";

    # c will contain the introns of the predicted transcript
    split(intlist[a[2]],c,",");

    # d will contain the introns of the reference transcript
    n=split(intlist[trmax],d,",");

    # go over the introns of the predicted transcript and for a given intron of the predicted tr
    # if it is absent from the list of introns of the ref tr then retain it either in the 5' or the 3' list
    # according to its position wrt the beg of the end of the ref tr and the strand of the locus
    l=1;
    while(c[l]!="")
    {
	found=0;
	m=1;
	while(found==0&&d[m]!="")
	{
	    if(c[l]==d[m])
	    {
		found=1;
	    }
	    m++;
	}
	# get in c1 the coordinates of the current intron of the predicted transcript
	split(c[l],c1,":");
	# get in d1 the coordinates of the 1st intron of the ref tr 
	split(d[1],d1,":");
	# get in d2 the coordinates of the last intron of the ref tr
	split(d[n],d2,":");

	# in case the current intron of the predicted transcript is not found
	if(found==0)
	{
	    if(((c1[1]<d1[1])&&(str[a[2]]=="+"))||(c1[2]>d1[2])&&(str[a[2]]=="-"))
	    {
		list5p=(list5p)(c[l])(",");
	    }
	    else
	    {
		if(((c1[1]<d1[1])&&(str[a[2]]=="-"))||(c1[2]>d1[2])&&(str[a[2]]=="+"))
		{
		    list3p=(list3p)(c[l])(",");
		}	
	    }
	}
	l++;
    }

    # in case we found no intron in 5' and/or in 3' we put correct strings and numbers in the variables to print
    if(list5p=="")
    {
	list5p="NA";
	nb5p=0;
    }
    else
    {
	nb5p=split(list5p,e,",")-1;
    }
    if(list3p=="")
    {
	list3p="NA";
	nb3p=0;
    }
    else
    {
	nb3p=split(list3p,f,",")-1;
    }
    
    print a[2], $3, trmax, str[a[2]], intlist[a[2]], intlist[trmax], nbint[a[2]], nbint[trmax], list5p, list3p, nb5p, nb3p;
}
