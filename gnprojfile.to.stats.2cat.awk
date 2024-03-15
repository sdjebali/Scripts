# gnprojfile.to.stats.2cat.awk
# from a tagiso gene file with indication of gene id, feelnc biotype, overlapwithens108
# and also 5 booleans for projectability to human genome, to human ens108 genes, in a 1 to 1 way
# to a single hs gene, on a hs gene associated to a single tagiso gene and then the hs gene in question
# make stats regarding these 5 projectability numbers for all genes, for genes of the 4 gene biotypes,
# for genes of the 2 classes of overlap with livestock ens108 genes and then for the 8 combined classes

# example
# cd ~/fragencode/workspace/geneswitch/results/rnaseq/multi/projection.to.human/geneswitch/sus_scrofa$
# indir=~/fragencode/workspace/geneswitch/results/rnaseq/multi
# pgm=~/fragencode/tools/multi/Scripts/gnprojfile.to.stats.2cat.awk
# sp=gallus_gallus
# p=tagiso
# gtf=/work/project/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/TAGADA.v2.1.1.GRCg6a.isoseq.23-01-13/annotation/novel.gtf
# cd $indir/projection.to.human/geneswitch/$sp
# base=`basename ${gtf%.gtf}`
# awk -v p=$p -v sp=$sp -f $pgm $p.$sp.$base.gn.id.bt.nov.hsproj.hsgn.1to1.1hsgn.1taggn.hsid1to1.tsv 

# $p.$sp.$base.gn.id.bt.nov.hsproj.hsgn.1to1.1hsgn.1taggn.hsid1to1.tsv 
# LOC_000000004151	mRNA	OverEns108	1	1	1	1	1	ENSG00000131653
# LOC_000000046247	lncRNA	NotOverEns108	0	0	0	0	0	NA
# 33932 (9 fields)

# output
# tagiso  gallus_gallus  allgn                 33932  21359  62.9465  14530  42.8209  11005  75.7398  12563  86.4625  12790  88.0248
# tagiso  gallus_gallus  mRNA                  15495  13493  87.0797  12858  82.9816  9963   77.4848  10972  85.3321  11682  90.8539
# tagiso  gallus_gallus  lncRNA                16948  7509   44.3061  1585   9.35214  997    62.9022  1504   94.8896  1063   67.0662
# tagiso  gallus_gallus  TUCp                  1423   347    24.3851  84     5.90302  44     52.381   84     100      44     52.381
# tagiso  gallus_gallus  noORF                 66     10     15.1515  3      4.54545  1      33.3333  3      100      1      33.3333
# tagiso  gallus_gallus  OverEns108            19168  15472  80.7179  13335  69.5691  10321  77.3978  11403  85.5118  12077  90.5662
# tagiso  gallus_gallus  NotOverEns108         14764  5887   39.874   1195   8.09401  684    57.2385  1160   97.0711  713    59.6653
# tagiso  gallus_gallus  mRNA.OverEns108       13933  12947  92.9233  12548  90.0596  9774   77.8929  10676  85.0813  11479  91.4807
# tagiso  gallus_gallus  mRNA.NotOverEns108    1562   546    34.9552  310    19.8464  189    60.9677  296    95.4839  203    65.4839
# tagiso  gallus_gallus  lncRNA.OverEns108     4996   2455   49.1393  749    14.992   523    69.8264  689    91.9893  574    76.6355
# tagiso  gallus_gallus  lncRNA.NotOverEns108  11952  5054   42.2858  836    6.99465  474    56.6986  815    97.488   489    58.4928
# tagiso  gallus_gallus  TUCp.OverEns108       238    69     28.9916  37     15.5462  24     64.8649  37     100      24     64.8649
# tagiso  gallus_gallus  TUCp.NotOverEns108    1185   278    23.4599  47     3.96624  20     42.5532  47     100      20     42.5532
# tagiso  gallus_gallus  noORF.OverEns108      1      1      100      1      100      0      0        1      100      0      0
# tagiso  gallus_gallus  noORF.NotOverEns108   65     9      13.8462  2      3.07692  1      50       2      100      1      50

BEGIN{
    c1[1]="mRNA";
    c1[2]="lncRNA";
    c1[3]="TUCp";
    c1[4]="noORF";
    c2[1]="OverEns108";
    c2[2]="NotOverEns108";
    for(i=1; i<=4; i++)
    {
	for(j=1; j<=2; j++)
	{
	    l++;
	    c3[l]=c1[i]"."c2[j];
	    nb30[c3[l]]=0;
	    nb31[c3[l]]=0;
	    nb32[c3[l]]=0;
	    nb33[c3[l]]=0;
	    nb34[c3[l]]=0;
	}
    }
}

{
    n++;
    nb1[$2]++;
    nb2[$3]++;
    nb3[$2"."$3]++;
    n0+=$4;
    n1+=$5;
    n2+=$6;
    n3+=$7;
    n4+=$8;
    nb10[$2]+=$4;
    nb11[$2]+=$5;
    nb12[$2]+=$6;
    nb13[$2]+=$7;
    nb14[$2]+=$8;
    nb20[$3]+=$4;
    nb21[$3]+=$5;
    nb22[$3]+=$6;
    nb23[$3]+=$7;
    nb24[$3]+=$8;
    nb30[$2"."$3]+=$4;
    nb31[$2"."$3]+=$5;
    nb32[$2"."$3]+=$6;
    nb33[$2"."$3]+=$7;
    nb34[$2"."$3]+=$8;
}

END{
    OFS="\t";
    print p, sp, "allgn", n, n0, div(n0,n), n1, div(n1,n), n2, div(n2,n1), n3, div(n3,n1), n4, div(n4,n1);
    for(i=1; i<=4; i++)
    {
	c=c1[i];
	print p, sp, c, nb1[c], nb10[c], div(nb10[c],nb1[c]), nb11[c], div(nb11[c],nb1[c]), nb12[c], div(nb12[c],nb11[c]), nb13[c], div(nb13[c],nb11[c]), nb14[c], div(nb14[c],nb11[c]);
    }
    for(j=1; j<=2; j++)
    {
	c=c2[j];
	print p, sp, c, nb2[c], nb20[c], div(nb20[c],nb2[c]), nb21[c], div(nb21[c],nb2[c]), nb22[c], div(nb22[c],nb21[c]), nb23[c], div(nb23[c],nb21[c]), nb24[c], div(nb24[c],nb21[c]);
    }
    for(i=1; i<=8; i++)
    {
	c=c3[i];
	print p, sp, c, nb3[c], nb30[c], div(nb30[c],nb3[c]), nb31[c], div(nb31[c],nb3[c]), nb32[c], div(nb32[c],nb31[c]), nb33[c], div(nb33[c],nb31[c]), nb34[c], div(nb34[c],nb31[c]);
    }
}

function div(x,y){return ((y!=0&&y!="") ? x/y*100 : "NA")}
