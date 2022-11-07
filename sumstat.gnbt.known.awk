
# sumstat.gnbt.known.awk

# example
#########
# cd /work/project/fragencode/results/rnaseq/bos_taurus/TAGADA.2.1.0.ARS-UCD1.2.102.22-04-20/annotation
# pgm=~/fragencode/tools/multi/Scripts/sumstat.gnbt.known.awk
# time awk -v fst=bos_taurus -v snd="tagada" -f $pgm novel.with.gnbt.gtf
# real	0m3.819s

# input novel.with.gnbt.gtf
# 3	tagada	gene	94278673	94314092	0	+	.	gene_id "LOC_000000013702"; ref_gene_id "."; feelnc_gene_biotype "mRNA";
# 3	tagada	transcript	94278673	94314092	0	+	.	gene_id "LOC_000000013702"; transcript_id "TM_000000238815"; contains "testis.filtered:STRG.8824.1,muscle.filtered:STRG.3536.2,cd4.filtered:STRG.6368.2,cd8.filtered:STRG.5188.2"; contains_count "4"; 3p_dists_to_3p "0,0,1567,1567"; 5p_dists_to_5p "0,390,211,211"; flrpm "1.334196"; longest "testis.filtered:STRG.8824.1"; longest_FL_supporters "cd8.filtered:STRG.5188.2,testis.filtered:STRG.8824.1,cd4.filtered:STRG.6368.2,muscle.filtered:STRG.3536.2"; longest_FL_supporters_count "4"; mature_RNA_length "3374"; meta_3p_dists_to_5p "1,1,0.53556609365738,0.53556609365738"; meta_5p_dists_to_5p "0,0.115589804386485,0.0625370480142264,0.0625370480142264"; rpm "1.334196"; spliced "1"; ref_gene_id "."; feelnc_gene_biotype "mRNA";
# 34409 (14 fields)
# 349460 (16 fields)
# 35446 (18 fields)
# 396330 (20 fields)
# 7381 (22 fields)
# 32963 (42 fields)
# 12595 (44 fields)
# 38727 (46 fields)
# 2329 (48 fields)

# output
# bos_taurus  tagada  allgnbt  34409  13561  39.4112  9214  26.7779  11634  33.8109
# bos_taurus  tagada  mRNA     20771  9267   44.6151  349   1.68023  11155  53.7047
# bos_taurus  tagada  lncRNA   7735   798    10.3167  6479  83.7621  458    5.92114
# bos_taurus  tagada  other    5903   3496   59.2241  2386  40.4201  21     0.355751


$3=="gene"{
    totgn++;
    split($10,a,"\"");
    split($NF,b,"\"");
    gnbt[a[2]]=b[2];
}

$3=="transcript"{
    split($10,a,"\"");
    split($12,b,"\"");
    trlist[a[2]]=(trlist[a[2]])(b[2])(",");
    split($0,c,"\t");
    split(c[9],d," ");
    k=1;
    while(d[k]!="")
    {
	if(d[k]=="transcript_biotype")
	{
	    known[b[2]]=1;
	}
	k+=2;
    }
}


END{
    OFS="\t";
    for(g in gnbt)
    {
	tot[gnbt[g]]++;
	nbtrkn=0;
	n=split(trlist[g],a,",");
	k=1;
	while(a[k]!="")
	{
	    if(known[a[k]])
	    {
		nbtrkn++;
	    }
	    k++;
	}

	if(nbtrkn==(n-1))
	{
	    alltrkn++;
	    atk[gnbt[g]]++;
	}
	else
	{
	    if(nbtrkn>0)
	    {
		mix++;
		m[gnbt[g]]++;
	    }
	    else
	    {
		alltrunkn++;
		atuk[gnbt[g]]++;
	    }
	}
    }

    print fst, snd, "allgnbt", totgn, alltrkn, alltrkn/totgn*100, alltrunkn, alltrunkn/totgn*100, mix, mix/totgn*100;
    c[1]="mRNA";
    c[2]="lncRNA";
    c[3]="other";
    for(i=1; i<=3; i++)
    {
	print fst, snd, c[i], tot[c[i]], atk[c[i]], divna(atk[c[i]],tot[c[i]])*100, atuk[c[i]], divna(atuk[c[i]],tot[c[i]])*100, m[c[i]], divna(m[c[i]],tot[c[i]])*100;
    }
}

function divna(x,y){
    return (y==0 ? "NA" : x/y);
}
