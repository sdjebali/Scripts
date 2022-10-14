# rnaataccorr.to.gene.with.peaknb.awk
# todo: take the threshold out of the script (pos 11 here)

# example
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp5/multi/atac.rna.correlation/tsspeaks
# pgm=~/fragencode/tools/multi/Scripts/rnaataccorr.to.gene.with.peaknb.awk
# time awk -v dir=pos -v t=0.7 -f $pgm tsspeaks1to1.gn.loess.tmm.norm.log10.twocorr.gnfirst.cat.gene.peak.tsv > rna.atac.corr.more0.7.gnid.totpeak.percat.percatandgene.tsv
# real	0m0.251s

# input: tsspeaks1to1.gn.loess.tmm.norm.log10.twocorr.gnfirst.cat.gene.peak.tsv 
# #chr	start	end	ID	group	chr	start	end	ID	group	pearson_r	spearman_r	peak.cat	peak.gntss
# 1	22151	22152	ENSSSCG00000037372	B	1	21547	24249	1:21547:24249	A	0.5950	0.6167	promoter-TSS	ENSSSCG00000037372
# 153956 (14 fields)


# output: rna.atac.corr.more0.7.gnid.totpeak.percat.percatandgene.tsv
# ENSSSCG00000017024	1	0	0	1	0	0	0	0	1	0	0	0	0	0	0	0
# ENSSSCG00000015781	1	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0
# 3456 (17 fields) 


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

(NR>=2)&&((dir=="pos"&&$11>=t)||(dir=="neg"&&$11<=t)){
    tot[$4]++;
    tot1[$4,$13]++;
    tot2[$4,$13,($14==$4 ? "same" : "diff")]++;
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
