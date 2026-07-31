# compute_jaccard_common_enh_gnpair.awk

# From a tsv filw with header that has for each gene pair and in the last two columns (8th and 9th)
# the comma separated lists of enhancers contacted by any of the two genes and by exactly the two genes
# computes and provides as output the cardinal of the union and of the intersection as well as the
# Jaccard percentage

# example
# cd ~/work/duplicons/common_enhancers
# pgm=~/fragencode/tools/multi/Scripts/compute_jaccard_common_enh_gnpair.awk
# awk -f $pgm supplementary_datasets7.closer4Mb.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.union.inter.tsv > supplementary_datasets7.closer4Mb.gn1.gn2.tss1.tss2.dist.contactedenhancers.union.inter.nb.jaccardpcent.tsv

# input file = supplementary_datasets7.closer4Mb.gn1.gn2.tss1.tss2.dist.nrenhokdist.gn1.gn2.union.inter.tsv 
# gnid1	gnid2	tsspos1	tsspos2	tssdist	nrenhdistok1	nrenhdistok2	enhunionlist	enhinterlist
# ENSMUSG00000001225	ENSMUSG00000020651	12:31390871	12:31559969	169098	31700960:31701110,31703360:31703510,31701360:31701510,31121860:31122010,31721560:31721710,31719260:31719410,31720560:31720710,31720960:31721410,31719560:31719710,31927460:31927610,	31261260:31261410,31268560:31268710,	31700960:31701110,31703360:31703510,31701360:31701510,31121860:31122010,31721560:31721710,31719260:31719410,31720560:31720710,31720960:31721410,31719560:31719710,31927460:31927610,31261260:31261410,31268560:31268710,	NA
# 330 (9 fields)

# output file = supplementary_datasets7.closer4Mb.gn1.gn2.tss1.tss2.dist.contactedenhancers.union.inter.nb.jaccardpcent.tsv
# gnid1	gnid2	tsspos1	tsspos2	tssdist	nunion	ninter	jaccard.pcent
# ENSMUSG00000001225	ENSMUSG00000020651	12:31390871	12:31559969	169098	12	0	0
# 330 (8 fields)

NR==1{
    OFS="\t";
    print $1, $2, $3, $4, $5, "nunion", "ninter", "jaccard.pcent";
}


NR>=2{
    nu=split($8,a,",");
    ni=split($9,b,",");
    if($9=="NA")
    {
	inter=0;
	union=($8=="NA" ? 0 : nu-1);
	jac=($8=="NA" ? "NA" : 0);
    }
    else
    {
	# here the intersection is not empty so in theory the union should not be empty
	if($8=="NA")
	{
	    union=0;
	    jac="pb";
	}
	else
	{
	    union=nu-1;
	    inter=ni-1;
	    jac=inter/union*100;
	}
    }
    print $1, $2, $3, $4, $5, union, inter, jac;
}
