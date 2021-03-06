# gnprom2uniqgnprom.awk
# note: it is possible to use this script with elt2 as well, which have an indication of type (elt2A and elt2B) at the end of $4
# of the input file, and for this reason we test for the existence of this last piece of info to be put in promlist

# example
# ct=HCmerge
# dir=~/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019
# pgm=/work/project/fragencode/tools/multi/Scripts/gnprom2uniqgnprom.awk
# cd $dir/$ct/score2/pall
# awk -f $pgm $ct.pall.score2.gene.of.prom.uniq.bed | sort -V -k1,1 -k2,2n -k3,3n > $ct.pall.score2.uniq.gene.of.prom.bed

# input
# chr1	1189288	1209265	ENSG00000160087.16:UBE2J2:chr1:1183838:1209657	.	-
# 14671 (6 fields)  *** with one row per gene,prom pair (prom is given at the end of $4), $1, $2 and $3 are the gene coordinates

# output
# chr1	1189288	1209265	ENSG00000160087.16:UBE2J2:chr1:1183838:1209657,	.	-
# 13923 (6 fields)  *** with one row per gene (with prom list at the end of $4), $1, $2 and $3 are the gene coordinates



{
    split($4,a,":");
    coord[a[1]":"a[2]]=$1"\t"$2"\t"$3;
    str[a[1]":"a[2]]=$6;
    addstr=(a[6]=="" ? "" : ":"a[6]);
    promlist[a[1]":"a[2]]=(promlist[a[1]":"a[2]])((a[3]":"a[4]":"a[5])(addstr))(",");
}

END{
    OFS="\t";
    for(g in promlist)
    {
	print coord[g], g":"promlist[g], ".", str[g];
    }
}
