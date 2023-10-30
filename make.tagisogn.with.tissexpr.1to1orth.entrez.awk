# make.tagisogn.with.tissexpr.1to1orth.entrez.awk

# usage
# srun --x11 --mem=8G --pty bash
# dir=~/fragencode/workspace/geneswitch/results/rnaseq
# sp=gallus_gallus
# tagisodir=nf-feelnc.tagiso.23-03-03
# pgm=~/fragencode/tools/multi/Scripts/make.tagisogn.with.tissexpr.1to1orth.entrez.awk
# cd $dir/$sp/$tagisodir/results/annotation
# quantif=$dir/$sp/$tagisodir/quantification/novel_genes_TPM.tsv 
# orth=~/fragencode/data/species/multi/homo_sapiens-$sp/ens102_orth_genes_homo_sapiens-$sp.tsv
# corr=~/fragencode/data/species/homo_sapiens/multi/gene.id.ens.ncbi.name.tsv
# time awk -v fileRef1=$quantif -v fileRef2=$orth -v fileRef3=$corr -f $pgm novel.gnid.coding.novelty.class.ens108refgnid.tsv

# fileRef1=$quantif 
# gene	cerebellum_stage1_rep1	cerebellum_stage1_rep2	...	skin_stage3_rep3	skin_stage3_rep4
# LOC_000000000000	12.04361	19.133354	...	1.641696	1.77886
# 33933 (85 fields)
# fileRef2=$orth
# Gene stable ID	Chicken gene stable ID	Chicken gene name	Chicken chromosome/scaffold name	Chicken chromosome/scaffold start (bp)	Chicken chromosome/scaffold end (bp)	Last common ancestor with Chicken	Chicken homology type	%id. target Chicken gene identical to query gene	%id. query gene identical to target Chicken gene
# ENSG00000198888	ENSGALG00010000007	ND1	MT	2824	3798	Amniota	ortholog_one2one	69.1824	67.9012
# 4240 (9 fields)
# 13413 (10 fields)
# 1 (45 fields)
# fileRef3=$corr
# Gene stable ID	NCBI gene (formerly Entrezgene) ID	Gene name
# ENSG00000210049		MT-TF
# 21039 (1 fields)
# 20189 (2 fields)
# 33890 (3 fields)
# 1 (10 fields)
# main input novel.gnid.coding.novelty.class.ens108refgnid.tsv
# LOC_000000004151	mRNA	known	ENSGALG00000005767
# 33932 (4 fields)   *** all the known genes have an ensembl gene id in last column (fine)

# output is 7 files such as this one for lung
# novel.gnid.ensid.expr012samp.tisslung.1to1hsid.entrezid.tsv
# LOC_000000004151	ENSGALG00000005767	1	ENSG00000131653	84231
# LOC_000000046247	NA	0	NA	NA
# 33932 (5 fields)

BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	n++;
	if(n==1)
	{
	    i=2;
	    while($i!="")
	    {
		split($i,a,"_");
		nb[a[1]]++;
		tiss[i]=a[1];
		i++;
	    }
	}
	else
	{
	    for(t in nb)
	    {
		ok[t]=0;
	    }
	    i=2;
	    while($i!="")
	    {
		if($i>=0.1)
		{
		    ok[tiss[i]]++;
		}
		i++;
	    }
	    for(t in nb)
	    {
		if(ok[t]>=2)
		{
		    expr[$1,t]=1;
		}
	    }
	}
    }

    n=0;
    while (getline < fileRef2 >0)
    {
	split($0,a,"\t");
	n++;
	if(n==1)
	{
	    found=0;
	    i=1;
	    while(found==0&&a[i]!="")
	    {
		if(a[i]~/homology\ type/)
		{
		    found=1;
		}
		i++;
	    }
	    if(found==1)
	    {
		col=(i-1);
	    }
	}
	else
	{
	    if(a[col]=="ortholog_one2one")
	    {
		hsid[$2]=$1;
	    }
	}
    }
    
    while (getline < fileRef3 >0)
    {
	split($0,a,"\t");
	if(a[2]!="")
	{
	    entrez[a[1]]=a[2];
	}
    }
}

{
    hs=(hsid[$4]!="" ? hsid[$4] : "NA");
    for(t in nb)
    {
	print $1, $4, (expr[$1,t]==1 ? 1 : 0), hs, (entrez[hs]!="" ? entrez[hs] : "NA") > "novel.gnid.ensid.expr012samp.tiss"t".1to1hsid.entrezid.tsv";
    }
}
