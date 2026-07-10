# add_pos_neg_contrib_to_GO_file.awk

# Given a file of enriched GO terms for BP, MF or CC obtained by Cervin's wrapper of topGO 
# and some operations on the main output file such as 
# - a different header
# - erasing last column
# - with underscores instead of spaces
# ie with the following information:
# GO.ID Term    Annotated       Significant     Expected        pval
# that was obtained on a list of genes that were the X best of spls comp no i,
# adds two columns with the gene names of the GO term  that were tested split 
# according to whether they contributed positively or negatively to the spls 
# component no i
# For this the script uses 3 different files used before reading the main input file
# - 3 col tsv file wo header with gene id, gene biotype and gene name for the universe
#   that here both serves as a correspondance file and to know which genes are expressed 
# - 2nd tsv file obtained by Cervin's wrapper of topGO that has for all go terms, 
#   not only the signif, the term id and the number and comma separated gene id 
#   list of the terms tested there
# - 2 column tsv file of all genes of the universe with gene name and signed 
#   contrib to comp1, with 0 when it is not among the top X genes of comp no i

# usage
# cd $basedir/results/go.analysis/long/males
# i=1
# go=BP
# time awk -v fileRef1=../../../../data/block.spls.quant.allgn.id.bt.name.tsv 
# -v fileRef2=hs_go_enrichment_male_sig_genes_per_term.tsv 
# -v fileRef3=../../../mixomics/22males.blockspls.long.allgn.comp1.signed.tsv 
# -f ../../../../scripts/go.analysis/long/add_pos_neg_contrib_to_GO_file.awk 
# hs_go_enrichment_male_$go.tsv > hs_go_enrichment_male_$go\_pos_neg_contribtocomp_gnnames.tsv
# real    0m1.809s

# fileRef1=../../../../data/block.spls.quant.allgn.id.bt.name.tsv is like this
# ENSG00000000003	protein_coding	TSPAN6
# 13682 (3 fields)

# fileRef2=hs_go_enrichment_male_sig_genes_per_term.tsv is like this
# GO_ID nb_sig_genes    sig_genes
# GO:0000002    2       ENSG00000143799,ENSG00000177302
# 9143 (2 fields)
# 11506 (3 fields)

# fileRef3=../../../mixomics/22males.blockspls.long.allgn.comp1.signed.tsv  is like this
# TSPAN6        0
# TNMD  0
# 13682 (2 fields)

# main input file = hs_go_enrichment_male_$go.tsv  is like this
# GO.ID Term    Annotated       Significant     Expected        pval
# GO:0006338    chromatin_remodeling    784     187     119.36  6.4e-07
# 56 (6 fields)

# main output file = hs_go_enrichment_male_$go\_pos_neg_contribtocomp_gnnames.tsv
# GO.ID Term    Annotated       Significant     Expected        pval    testedgn.withposcontrib.oncomp  testedgn.withnegcontrib.oncomp
# GO:0006338    chromatin_remodeling    784     187     119.36  6.4e-07 CREBBP,KMT2E,BAZ1B,MRE11,UBR2,ARID4A,JADE2,TPR,RSF1,ARID1B,ARID4B,USP36,KMT2C,WNK1,SPEN,SLK,ROCK1,RPS6KA2,PTPN18,MARK2,SMARCE1,KDM5A,NUAK1,BAZ2A,SMARCA2,RSBN1,PTPRC,CDC14B,PHLPP1,KAT6A,YTHDC1,REST,NCOA1,ATRX,MECOM,KDM2B,CBX5,SRPK1,MKNK2,MTMR3,EP300,YY1,CHD8,STK4,RIOK3,STK3,UBR5,TASOR2,SUPT6H,BCL7A,CHD4,SUDS3,HCFC2,CILK1,PRP4K,ARB2A,RAD50,CSNK1A1,NEK11,PPM1G,INO80B,ASH1L,PRDM2,RBBP5,ARID1A,KMT2A,KDM3B,RPAP2,BAZ2B,NCOA3,CHD6,DEK,ZBTB1,KDM4B,PTPN12,EIF2AK4,TTBK2,INO80,MINDY2,PARP2,RAF1,PRMT7,MBD2,ROCK2,CLOCK,SRPK2,USP15,KAT7,CDK9,ANP32B,NEK1,PPP3CA,SMARCC2,PML,PDPK1,BRD4,MTF2,SHPRH,TLK2,NSD3,ATM,TAOK2,SMARCA5,CFDP1,TRIP12,CHD1,TTN,RNF20,USP16,ATAD2,SMG1,GATAD1,BRAF,KALRN,PRMT2,TAOK1,LMNA,TAF6L,MYSM1,IWS1,PBRM1,HMGB2,NIPBL,UBLCP1,STK17A,PSIP1,NSD1,BRD7,IGF2,ING2,MECP2,CSNK1G1,BRD3,PWWP2A,NFRKB,CHD7,BPTF,JMJD1C,KDM2A,SMARCC1,CHD2,MTHFR,SUZ12,ZBTB7A,PTPN11,SETD2,KMT5A,TSPYL2,H2AX,LRRK2,ARID2,MPHOSPH8,SUPT5H,PTPN1,H4C3,MYO1C,MIER1,RPS6KL1,TLK1,BAZ1A,CTR9,TOP1,ZDBF2,H2AJ,PRKDC,KCNQ1OT1,KMT2B,     HPF1,HDAC4,SIRT6,PRMT5,PTP4A1,DUSP4,CSNK1G2,PSKH1,YEATS2,PHKG1,NUDT5,PHB1,TP53RK,NAP1L5,METTL23,MTOR,SUPT4H1,PPM1N,H2BC15,DUSP14,
# 56 (8 fields)

BEGIN{
    OFS="\t"; 
    # ENSG00000000003	protein_coding	TSPAN6
    while (getline < fileRef1 >0)
    {
        name[$1]=$3;
    } 
    
    # GO:0000002    2       ENSG00000143799,ENSG00000177302
    while (getline < fileRef2 >0)
    {
        split($0,a,"\t"); 
        if(a[2]>0)  # added after test so need to retest
        {
            split(a[3],b,","); 
            k=1; 
            while(b[k]!="")
            {
                seen[a[1],name[b[k]]]++; 
                if(seen[a[1],name[b[k]]]==1)
                {
                    nbgn[a[1]]++; 
                    gn[a[1],nbgn[a[1]]]=name[b[k]];
                } 
                k++;
            }
        }   
    } 

    # TSPAN6        0
    while (getline < fileRef3 >0)
    {
        if($2!=0)
        {
            ok[$1]=1;
            pos[$1]=($2>0 ? 1 : 0);
        }
    }
} 

# GO.ID Term    Annotated       Significant     Expected        pval
# GO:0006338    chromatin_remodeling    784     187     119.36  6.4e-07
{
    split($0,a,"\t"); 
    if(NR==1)
    {
        print $0, "testedgn.withposcontrib.oncomp", "testedgn.withnegcontrib.oncomp";
    } 
    else
    {
        st=$0"\t"; 
        spos=""; 
        sneg=""; 
        for(k=1; k<=nbgn[a[1]]; k++)
        {
            if(ok[gn[a[1],k]]==1)
            {
                if(pos[gn[a[1],k]]==1)
                {
                    spos=(spos)(gn[a[1],k])(",");
                }
                else
                {
                    sneg=(sneg)(gn[a[1],k])(",");
                }
            }
        } 
        print (st)(nn(spos))"\t"(nn(sneg));
    }
} 
    
    function nn(x){
        return (x=="" ? "NA" : x);
    }