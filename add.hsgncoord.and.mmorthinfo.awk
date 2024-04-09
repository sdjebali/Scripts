# add.hsgncoord.and.mmorthinfo.awk
# to a file of genes activated by same enhancers than 13 genes important in iron from ABC model and phic predictions
# from 5 tissues (liver and intestine related), add information about the chr of gene 1 and gene 2 in human GRCh38
# as well as the same for the mouse orth on the mouse genome when it exists as 1 to 1 (and from mm39) and then the
# distance between the two genes in human and in mouse when they exist and are on the same chr


# example
# cd ~/regenet/workspace/sdjebali/iron/eg.graph.to.find.new.genes
# hs=~/fragencode/data/species/homo_sapiens/GRCh38.108/homo_sapiens.gtf
# orth=~/fragencode/data/species/multi/homo_sapiens-mus_musculus/ens11_orth_genes_homo_sapiens-mus_musculus.tsv
# pgm=~/fragencode/tools/multi/Scripts/add.hsgncoord.and.mmorthinfo.awk
# time awk -v fileRef1=$hs -v fileRef2=$orth -f $pgm new.genes.from.13genes.abc5ct.chicliver.tsv > new.genes.from.13genes.abc5ct.chicliver.chr.gn1.gn2.hs.mm.dist.hs.mm.tsv
# *** real	0m3.159s

# fileRef1=$hs
# #!genome-build GRCh38.p13
# #!genome-version GRCh38
# #!genome-date 2013-12
# #!genome-build-accession GCA_000001405.28
# #!genebuild-last-updated 2022-07
# 1	ensembl_havana	gene	1471765	1497848	.	+	.	gene_id "ENSG00000160072"; gene_version "20"; gene_name "ATAD3B"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
# 1	ensembl_havana	transcript	1471765	1497848	.	+	.	gene_id "ENSG00000160072"; gene_version "20"; transcript_id "ENST00000673477"; transcript_version "1"; gene_name "ATAD3B"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "ATAD3B-206"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select";

# fileRef2=$mm 
# Gene stable ID	Gene start (bp)	Gene end (bp)	Strand	Gene name	Chromosome/scaffold name	Mouse gene stable ID	Mouse gene name	Mouse chromosome/scaffold name	Mouse chromosome/scaffold start (bp)	Mouse chromosome/scaffold end (bp)	Mouse homology type	%id. target Mouse gene identical to query gene	%id. query gene identical to target Mouse gene	Mouse orthology confidence [0 low, 1 high]
# ENSG00000198888	3307	4262	1	MT-ND1	MT	ENSMUSG00000064341	mt-Nd1	MT	2751	3707	ortholog_one2one	77.044	77.044	1
# 5 (13 fields)
# 901 (14 fields)
# 25454 (15 fields)
# 1 (58 fields)

# main input file: new.genes.from.13genes.abc5ct.chicliver.tsv
# gene	gene.id	score.adj	gene.source	score.tot	abc.enhancer.intersect.label	ABC.product.label	chic.enhancer.intersect.label	contact.product.label	count	datatype
# COQ9	ENSG00000088682	15	CIAPIN1	31	2	3	2	3	21	ABC,CHiC
# 454 (11 fields)  

# main output new.genes.from.13genes.abc5ct.chicliver.chr.gn1.gn2.hs.mm.dist.hs.mm.tsv
# gene	gene.id	score.adj	gene.source	score.tot	abc.enhancer.intersect.label	ABC.product.label	chic.enhancer.intersect.label	contact.product.label	count	datatype	hs.gn1.chr	hs.gn2.chr	mm.gn1.chr	mm.gn2.chr	hs.gn12.dist	mm.gn12.dist
# COQ9	ENSG00000088682	15	CIAPIN1	31	2	3	2	3	21	ABC,CHiC	16	16	8	8	16544	17527
# 454 (17 fields)    *** with 154 rows with NA so 1/3 of the genes

BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	if($3=="gene")
	{
	    split($0,a,"\t");
	    split(a[9],b," ");
	    k=1;
	    found=0;
	    while(found==0&&b[k]!="")
	    {
		if(b[k]=="gene_name")
		{
		    split(b[k+1],c,"\"");
		    found=c[2];
		}
		k+=2;
	    }
	    
	    if(found!=0)
	    {
		split($10,d,"\"");
		hchr[d[2]]=$1;
		hbeg[d[2]]=$4;
		hend[d[2]]=$5;
		hchr[found]=$1;
		hbeg[found]=$4;
		hend[found]=$5;
	    }
	}
    }

    while (getline < fileRef2 >0)
    {
	split($0,a,"\t");
	if(a[12]=="ortholog_one2one")
	{
	    orth[a[5]]=a[8];
	    mchr[a[5]]=a[9];
	    mbeg[a[5]]=a[10];
	    mend[a[5]]=a[11];
	}
    }
    hchr["HFE2"]=1;
    hbeg["HFE2"]=146017470;
    hend["HFE2"]=146021735;
}

NR==1{
    print $0, "hs.gn1.chr", "hs.gn2.chr", "mm.gn1.chr", "mm.gn2.chr", "hs.gn12.dist", "mm.gn12.dist";
}

NR>=2{
    hdist=((hchr[$1]!=""&&hchr[$4]!=""&&hchr[$1]==hchr[$4]) ? (abs((hbeg[$1]+hend[$1])/2-(hbeg[$4]+hend[$4])/2)) : ((hchr[$2]!=""&&hchr[$4]!=""&&hchr[$2]==hchr[$4]) ? (abs((hbeg[$2]+hend[$2])/2-(hbeg[$4]+hend[$4])/2)) : "NA"));
    mdist=((mchr[$1]!=""&&mchr[$4]!=""&&mchr[$1]==mchr[$4]) ? (abs((mbeg[$1]+mend[$1])/2-(mbeg[$4]+mend[$4])/2)) : ((mchr[$2]!=""&&mchr[$4]!=""&&mchr[$2]==mchr[$4]) ? (abs((mbeg[$2]+mend[$2])/2-(mbeg[$4]+mend[$4])/2)) : "NA"))
    print $0, (hchr[$1]!="" ? hchr[$1] : (hchr[$2]!="" ? hchr[$2] : "NA")), nn(hchr[$4]), nn(mchr[$1]), nn(mchr[$4]), hdist, mdist;
}


function nn(x){
	return (x!="" ? x : "NA");
    }

function abs(x){
    return (x>=0 ? x : (-1)*x);
}
