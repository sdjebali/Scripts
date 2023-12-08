# chinputtech2biorep.awk

# takes as input the concatenation of chinput files without headers representing technical replicate of a sample
# with indication of number of occurence in column 3 and other coordinates and distance in columns 1, 2, 4 and 5
# and produces a chinput file with a single time a given bait-rf pair but with a number of occurence
# taking into account the one seen in each tech rep (so adding them all for a given pair)
# it also puts the header back

{
    nbreads[$1":"$2":"$4":"$5]+=$3;
}

END{
    OFS="\t";
    print "baitID", "otherEndID", "N", "otherEndLen", "distSign";
    for(id in nbreads)
    {
	split(id,a,":");
	print a[1], a[2], nbreads[id], a[3], a[4];
    }
}
