# pchic.to.sc.dist.type.category.awk
# given a bedpe file of pchic relations with indication of:
###########################################################
# - score in 7th column
# - p-p or p-o type in 11th column
# - one category out of 16 that gather ones that concern overlap of prom frag, of other end frag or of both
#   from columns 19 to 33
# provides a tsv file with header for ggplot2 that has:
#######################################################
# - score
# - distance between the two fragments (from middle to middle)
# - type (pp or po, how it is written in col 11)
# - category from 16 with a number before to have them in the correct order in the plot

# example of usage
##################
# cd ~/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018
# pgm=~/fragencode/tools/multi/Scripts/pchic.to.sc.dist.type.category.awk
# time awk -f $pgm pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.bedpe > pchic.ipsccm.montefiori2018.score.dist.class.16catwrtplsels.tsv

# input file = pchic.ipsccm.montefiori2018.hg19.overgencv49tss.part1.part2.reltype.gxorder.sorted.overpls.els.dels.first.second.part.16catwrtplsels.bedpe
# chr1	834807	835706	chr1	893041	896871	6.6	NOC2L*NM	1	0	p-o	rev.order	0	0	0	1	1	0	1	0	1	0	0	0	00	0	0	0	0	0	1	0	0
# 401098 (34 fields)

# output file = pchic.ipsccm.montefiori2018.score.dist.class.16catwrtplsels.tsv
# score	dist	class	category
# 6.6	59699	p-o	01_pp
# 2 024 919 (4 fields) 

BEGIN{
    OFS="\t";
    cat[1]="01_pp";
    cat[2]="02_op";
    cat[3]="03_pe";
    cat[4]="04_oe";
    cat[5]="05_pde";
    cat[6]="06_ode";
    cat[7]="07_ppoe";
    cat[8]="08_ppode";
    cat[9]="09_peop";
    cat[10]="10_pdeop";
    cat[11]="11_ppop";
    cat[12]="12_peoe";
    cat[13]="13_pdeode";
    cat[14]="14_pona";
    cat[15]="15_pnao";
    cat[16]="16_pnaona";
    print "score", "dist", "class", "category";
}

{
    if($1==$4)
    {
	d=dist($2,$3,$5,$6);
    }
    else
    {
	d="NA";
    }
    for(i=1; i<=16; i++)
    {
	if($(18+i)==1)
	{
	    print $7, d, $11, cat[i];
	}
    }
}

# computes distance between two fragments of (x,y) and (z,t) coordinates on the same chromosome
# as middle to middle
function dist(x,y,z,t){
    m1=(x+y)/2;
    m2=(z+t)/2;
    return ((m2-m1) > 0 ? int(m2-m1) : int(m1-m2));
}
