# transdecoder.gff3.to.gtf.awk
# this is to make a gtf with only CDS and exon rows from the final gff3 file from transdecoder
# the same way as gffread from cufflinks did it but when the files have homology information
# since in this case gffread does not work even when increasing the ram
# uses only mrna, cds and exon rows
# and takes information of gene id, score and type of orf in mrna rows
# it also provides transcript_id and then gene_id and in gene_name the type of orf and score
# be careful of the mrna rows with homology since they have a much more complete form with spaces
# in this case just take the orf type and score otherwise too big 

# srun --x11 --mem=16G --pty bash
# cd /work/project/fragencode/workspace/geneswitch/results/rnaseq/gallus_gallus/nf-feelnc.tagiso.23-03-03/results/annotation/transdecoder
# pgm=/work/project/fragencode/tools/multi/Scripts/transdecoder.gff3.to.gtf.awk
# time awk -f $pgm tagisotr.fasta.transdecoder.genome.gff3 | sort -V -k1,1 -k10,10 -k3,3 > tagisotr.fasta.transdecoder.genome.gtf
# real	5m37.416s

# tagisotr.fasta.transdecoder.genome.gff3
# !!! be careful in mrna rows only 3 subfields in 9th field but then if splitting on = dangerous since we have another = inside the last string !!!
# 2	transdecoder	mRNA	58581620	58867438	.	-	.	ID=TM_000000938490.p1;Parent=LOC_000000000000^2^-;Name="ORF type:5prime_partial (+),score=28.34"
# 2	transdecoder	CDS	58867299	58867438	.	-	0	ID=cds.TM_000000938490.p1;Parent=TM_000000938490.p1
# 2	transdecoder	exon	58867299	58867438	.	-	.	ID=TM_000000938490.p1.exon1;Parent=TM_000000938490.p1
# 24	transdecoder	mRNA	5405480	5407228	.	-	.	ID=G17684.3.p1;Parent=LOC_000000018251^24^-;Name="ORF type:complete (+),score=33.17,UniRef90_UPI00129D2F60|96.3|4.46e-28,UniRef90_A0A7K9YKV7|90.7|5.07e-26,UniRef90_UPI001C66E67B|90.7|3.30e-24,UniRef90_A0A7L3J8F9|76.6|3.89e-24,UniRef90_A0A7K9UNK3|85.2|1.92e-23,UniRef90_A0A7L4I703|64.4|4.46e-23,UniRef90_A0A7K7NSY9|85.2|5.45e-23,UniRef90_A0A7K5H0A4|85.2|7.73e-23,UniRef90_A0A7L1XYQ2|85.2|1.27e-22,UniRef90_A0A7K5URQ9|73.3|1.27e-22,UniRef90_A0A7K7UMB5|65.8|2.55e-22,UniRef90_A0A7K9SH62|84.9|3.12e-22,UniRef90_A0A7L4N8X1|85.2|4.42e-22,UniRef90_A0A7K6NV84|83.3|6.26e-22,UniRef90_A0A7L4MCQ6|69.0|7.24e-22,UniRef90_A0A7L1L412|81.5|1.03e-21,UniRef90_A0A851ZCP9|76.2|1.04e-21,UniRef90_A0A7L0RS88|81.5|1.11e-21,UniRef90_A0A7L2Z472|81.5|1.45e-21,UniRef90_UPI0022875E70|87.0|2.19e-21,UniRef90_A0A7L1MBZ0|75.9|3.58e-21,UniRef90_A0A7L2V2H8|79.6|3.58e-21,UniRef90_A0A7L0DBC3|79.6|4.13e-21,UniRef90_A0A851CNC2|77.8|5.07e-21,UniRef90_A0A7L4D1N4|79.6|7.18e-21,ATP1G1_PLM_MAT8|PF02038.21|5.7e-11"
# 21532 (0 fields)  
# 5709947 (9 fields)  
# 266500 (11 fields)

# tagisotr.fasta.transdecoder.genome.gtf
# 1	transdecoder	CDS	27072	27346	.	-	2	transcript_id "G3.1.p2"; gene_id "LOC_000000090338^1^-"; gene_name "ORF type:5prime_partial (+),score=17.47";
# 1	transdecoder	CDS	27608	27767	.	-	0	transcript_id "G3.1.p2"; gene_id "LOC_000000090338^1^-"; gene_name "ORF type:5prime_partial (+),score=17.47";
# 4300073 (16 fields)

$3=="mRNA"{
    trid="";
    split($0,a,"\t");
    split(a[9],b,";");
    k=1;
    while(b[k]!="")
    {
	split(b[k],c,"=");
	if(c[1]=="ID")
	{
	    trid=c[2]; 
	}
	else
	{
	    if(c[1]=="Parent"&&trid!="")
	    {
		    gnid[trid]=c[2];
	    }
	    else
	    {
		# for gene_name just take the orf type and score, not the homology otherwise file too big
		if(c[1]=="Name"&&trid!="")
		{
		    split(c[3],d,",");
		    gnname[trid]=c[2]("="d[1]);
		}
	    }
	}
	k++;
    }
}

# 2	transdecoder	CDS	58867299	58867438	.	-	0	ID=cds.TM_000000938490.p1;Parent=TM_000000938490.p1
# 2	transdecoder	exon	58867299	58867438	.	-	.	ID=TM_000000938490.p1.exon1;Parent=TM_000000938490.p1
$3=="CDS"||$3=="exon"{
    trid="";
    split($0,a,"\t");
    split(a[9],b,";");
    k=1;
    while(b[k]!="")
    {
	split(b[k],c,"=");
	if(c[1]=="ID")
	{
	    if($3=="CDS")
	    {
		split(c[2],d,"cds.");
		trid=d[2];
	    }
	    else
	    {
		if($3=="exon")
		{
		    split(c[2],d,".exon");
		    trid=d[1];
		}
	    }
	}
	k++;
    }
    if(trid!="")
    {
	if($3=="CDS")
	{
	    nbcds[trid]++;
	    cdsrow[trid,nbcds[trid]]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;
	}
	else
	{
	    if($3=="exon")
	    {
		nbex[trid]++;
		exrow[trid,nbex[trid]]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;
	    }
	}
    }
}

END{
    for(t in gnid)
    {
	for(i=1; i<=nbex[t]; i++)
	{
	    print exrow[t,i]"\ttranscript_id \""(t)"\"; gene_id \""(gnid[t])"\"; gene_name "(gnname[t])"\";";
	}
	for(i=1; i<=nbcds[t]; i++)
	{
	    print cdsrow[t,i]"\ttranscript_id \""(t)"\"; gene_id \""(gnid[t])"\"; gene_name "(gnname[t])"\";";
	}
    }
}
