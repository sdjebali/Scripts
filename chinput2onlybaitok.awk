# chinput2onlybaitok.awk

# takes as fileRef a baitmap file with the identifiers of the baits in column 4
# and a chinput file with the bait id in column 1 and only selects the rows with a bait id in fileRef


BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	ok[$4]=1;
    }
}

$1!="baitID"&&ok[$1]==1{
    print
}
