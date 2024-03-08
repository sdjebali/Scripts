# phicconnoveratac.to.simplertsv.awk

# starts from a file of phic connections with header that has the conection id as coord separated by :
# and then cis category (1.close, 2.far or 3.trans) and then nbtissstage also separated by :
# with indication of overlap of atac-seq peaks on either end and both ends
# and provides a simpler tsv file for a ggplot2 stacked barplot that gives a combined class for both the number
# of tissue, stage and the cis class (from lower to higher number of tissue, stage and with cis class in order)
# and then the atac.either class (as 1.ok or 2.ko) and the number of connections of this class

# srun --x11 --mem=32G --pty bash
# basedir=~/fragencode/workspace/geneswitch/results/hic
# sp=gallus_gallus
# outdir=$basedir/$sp/chicago
# cd $outdir/multi
# pgm=~/fragencode/tools/multi/Scripts/phicconnoveratac.to.simplertsv.awk
# time awk -v sp=$sp -f $pgm alltissstage.conn.overatac.eitherend.bothends.tsv > alltissstage.nbtsissciscat.ataceitherclass.nbconn.tsv

# alltissstage.conn.overatac.eitherend.bothends.tsv 
# connid	atac.eitherend	atac.bothends
# 1:7193:8745:1:8911:9717:1.close:1	0	1
# 2012580 (3 fields)

# alltissstage.nbtsissciscat.ataceitherclass.nbconn.tsv
# class	atac.either	nbconn
# 1_1_1.close	1.atac.either.ok	84944
# 37 (3 fields)

NR==1{
    OFS="\t";
    print "class", "atac.either", "nbconn";
    imax=(sp=="gallus_gallus" ? 6 : 7);
    for(i=1; i<=imax; i++)
    {
	no[i"_1.close"]=1+(i-1)*3;
	no[i"_2.far"]=2+(i-1)*3;
	no[i"_3.trans"]=3+(i-1)*3;
    }
}


NR>=2{
    split($1,a,":");
    nbclass=no[a[8]"_"a[7]];
    tot[nbclass"_"a[8]"_"a[7]]++;
    either[nbclass"_"a[8]"_"a[7]]+=$2;
}


END{
    for(i=1; i<=imax; i++)
    {
	print twodigit(no[i"_1.close"])"_"i"_1.close", "1.atac.either.ok", either[no[i"_1.close"]"_"i"_1.close"];
	print twodigit(no[i"_1.close"])"_"i"_1.close", "2.atac.either.ko", tot[no[i"_1.close"]"_"i"_1.close"]-either[no[i"_1.close"]"_"i"_1.close"];
	print twodigit(no[i"_2.far"])"_"i"_2.far", "1.atac.either.ok", either[no[i"_2.far"]"_"i"_2.far"];
	print twodigit(no[i"_2.far"])"_"i"_2.far", "2.atac.either.ko", tot[no[i"_2.far"]"_"i"_2.far"]-either[no[i"_2.far"]"_"i"_2.far"];
	print twodigit(no[i"_3.trans"])"_"i"_3.trans", "1.atac.either.ok", either[no[i"_3.trans"]"_"i"_3.trans"];
	print twodigit(no[i"_3.trans"])"_"i"_3.trans", "2.atac.either.ko", tot[no[i"_3.trans"]"_"i"_3.trans"]-either[no[i"_3.trans"]"_"i"_3.trans"];
    }
} 

function twodigit(x){
    if(x<10)
    	return "0"x;
    else
	return x;
}
