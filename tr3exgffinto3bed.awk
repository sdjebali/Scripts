# tr3exgffinto3bed.awk

# example
# sp=gallus_gallus
# cd ~/fragencode/workspace/geneswitch/results/rnaseq/$sp/tagiso/annotation
# pgm=~/fragencode/tools/multi/Scripts/tr3exgffinto3bed.awk
# time zcat novel.clean.full.tr3ex.exons.sortbytr.gff.gz | awk -f $pgm
# real	0m24.029s
# for c in init intern term
# do
# c novel.clean.full.tr3ex.lncmrna.exons.$c.bed
# done

# 3 outputs like this (could have redundance and are not ordered according to chr, start, end)
# 4	73527045	73527251	TM_000001742337.1	.	+
# 7	30417236	30417544	G30820.2.8	.	-
# 152655 (6 fields)
# 4	73528081	73528236	TM_000001742337.2	.	+
# 4	73529793	73529912	TM_000001742337.3	.	+
# 1363301 (6 fields)
# 4	73547924	73550833	TM_000001742337.11	.	+
# 7	30340171	30341292	G30820.2.1	.	-
# 152655 (6 fields) *** same nb in init and term

BEGIN{
    OFS="\t";
}

{
    split($12,a,"\"");
    split($NF,b,"\"");
    if(b[2]=="mRNA"||b[2]=="lncRNA")
    {
	bt[a[2]]=b[2];
	nbex[a[2]]++;
	str[a[2]]=$7;
	row[a[2],nbex[a[2]]]=$0;
    }
}

END{
    for(t in bt)
    {
	for(i=1; i<=nbex[t]; i++)
	{
	    split(row[t,i],a,"\t");
	    if((str[t]=="+"&&i==1)||(str[t]=="-"&&i==nbex[t]))
	    {
		c="init";
	    }
	    else
	    {
		if((str[t]=="+"&&i==nbex[t])||(str[t]=="-"&&i==1))
		{
		    c="term";
		}
		else
		{
		    c="intern";
		}
	    }
	    print a[1], a[4]-1, a[5], t"."i, ".", str[t] > "novel.clean.full.tr3ex.lncmrna.exons."c".bed";
	}
    }
}
