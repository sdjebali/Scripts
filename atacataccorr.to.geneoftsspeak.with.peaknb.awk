# atacataccorr.to.geneoftsspeak.with.peaknb.awk
# todo: take the threshold out of the script (pos 9 here)


# example
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp5/atacseq/peaks/pairwise.correlation
# pgm=~/fragencode/tools/multi/Scripts/atacataccorr.to.geneoftsspeak.with.peaknb.awk
# time awk -v dir=pos -v t=0.7 -f $pgm peaks.loess.each.sample.corr.1Mb.cat.gene.peak1.2.tsv > atac.atac.corr.more0.7.tsspeak.gnid.totpeak.percat.percatandgene.tsv
# real	0m52.217s


# input: peaks.loess.each.sample.corr.1Mb.cat.gene.peak1.2.tsv 
# #chr	start	end	ID	chr	start	end	ID	pearson_r	spearman_r	peak1.cat	peak1.gntss	peak2.cat	peak2.gntss
# 1	21547	24249	1:21547:24249	1	37863	38387	1:37863:38387	0.4977	0.5124	promoter-TSS	ENSSSCG00000037372	intron	ENSSSCG00000027257
# 32891977 (14 fields) 

# output: atac.atac.corr.more0.7.tsspeak.gnid.totpeak.percat.percatandgene.tsv
# ENSSSCG00000015780	2	1	0	1	0	0	1	0	0	0	0	0	0	1	0	0
# ENSSSCG00000015781	24	15	3	4	2	0	1	0	0	0	0	14	3	4	2	0
# 7239 (17 fields)


BEGIN{
    OFS="\t";
    if(dir=="")
    {
	dir="pos"
    }
    if(t=="")
    {
	t=0.7;
    }
    cat[1]="intron";
    cat[2]="Intergenic";
    cat[3]="promoter-TSS";
    cat[4]="exon";
    cat[5]="TTS";
}

(NR>=2)&&((dir=="pos"&&$9>=t)||(dir=="neg"&&$9<=t)){
    if($11=="promoter-TSS")
    {
	tot[$12]++;
	tot1[$12,$13]++;
	tot2[$12,$13,($14==$12 ? "same" : "diff")]++;
    }
    if($13=="promoter-TSS")
    {
	tot[$14]++;
	tot1[$14,$11]++;
	tot2[$14,$11,($12==$14 ? "same" : "diff")]++;
    }
}

END{
    for(g in tot)
    {
	s1="";
	s21="";
	s22="";
	for(i=1; i<=5; i++)
	{
	    s1=(s1)(nn(tot1[g,cat[i]]))("\t")
	}
	for(i=1; i<=4; i++)
	{
	    s21=(s21)(nn(tot2[g,cat[i],"same"]))("\t");
	    s22=(s22)(nn(tot2[g,cat[i],"diff"]))("\t");
	}
	print g, tot[g], (s1)(s21)(nn(tot2[g,cat[5],"same"]))("\t")(s22)(nn(tot2[g,cat[5],"diff"]));
    }
}

function nn(x){
    return (x!="" ? x : 0);
}
