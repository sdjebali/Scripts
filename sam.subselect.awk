# sam.subselect.awk
# this takes as input a sam file (with header) and as fileRef a bed file of features
# and makes a sam file that only has the chr present in fileRef

# example
# time

BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	ok[$1]=1
    }
}

$1~/@/{
    if($1=="@SQ")
    {
	split($2,a,":");
	if(ok[a[2]]==1)
	{
	    print
	}
    }
    else
    {
	print
    }
}

$1!~/@/{
    if(ok[$3]==1)
    {
	print
    }
}
