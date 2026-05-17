# make_duplgnfile_with_deltaexpr_expected_predicted.awk
# this script takes as input two tsv files with our paralogs
# - one with header that is read beforehand and that has at least the genes of our paralogs, one gene on each row, with information about ageclass and expected and predicted TPM
# - one with header that is the main input file and that has each pair of paralog on a different row and where the two pieces of information used here are the two gene ids in col 1 and 2
# and it makes a tsv file with header that has each gene pair of our paralogs on each row with indication of gene ids, age class, expected and predicted expression of each gene
# expected and predicted average expressions of the two genes and exected and predicted delta expressions of the two genes

# example of usage
# cd $dir/$m/kallisto
# paralogs=filtered_paralogs_dataset.tsv
# awk -v fileRef=Mus_musculus.GRCm38.cdna.realparalogs.gn.id.bt.inourparalogs.agecl.tr.list.nb.tpm.exp.kall.tsv -f make_duplgnfile_with_deltaexpr_expected_predicted.awk $paralogs > Mus_musculus.GRCm38.cdna.realparalogs.gnid.expkallexpr.1.2.expkallavgexpr.expkallnormdeltaexpr.tsv

# fileRef=Mus_musculus.GRCm38.cdna.realparalogs.gn.id.bt.inourparalogs.agecl.tr.list.nb.tpm.exp.kall.tsv 
# gnid	gnbt	inparalogs	ageclass	trlist	nbtr	expected.TPM	kallisto.TPM
# ENSMUSG00000030510.11	protein_coding	1	2.old	ENSMUST00000066475.10,ENSMUST00000208521.1,	2	11.3762	10.6057
# 2841 (8 fields)

# main input file = $paralogs = filtered_paralogs_dataset.tsv
# GeneID	ParalogueGeneID	LastCommonAncestor	ParalogueType	Id_ID2onID1	Id_ID1onID2	Term_DupliNode	root_id	ohnologs	class_age	protein_coding_1	protein_coding_2	both_protein_coding	targeted_1	targeted_2	Bait_targeting_ID1	Bait_targeting_ID2	Nbr_Bait_targeting_ID1	Nbr_Bait_targeting_ID2	nbr_common_bait	Contact_gene_enhancers	MaxExonsGene1MaxExonsGene2	class_exons_number	exons_difference_ratio	Gene1_chr	Gene2_chr	relationship	Gene1_strand	Gene2_strand	Gene1_TSS	Gene2_TSS	Gene1_Start	Gene2_Start	Gene1_End	Gene2_End	DistanceTSS	distance_between_genes	Not_a_retrotranscription_duplication	MeanTPM_GeneID	MeanTPM_ParalogGeneID	MeanTPM_pair	class_perfivetil	ID_dist	ID_expr	ID_dist_expr
# ENSMUSG00000032079	ENSMUSG00000032080	Gnathostomata	within_species_paralog	18.2065	16.962	yes	608290	non	very_old	yes	yes	yes	yes	yes	9:46267916-46274007	9:46240462-46245005	1	1	0	yes	4	3	multi-diff	0.25	9	9	cis	1	1	46268633	46240696	46268633	46240696	46271919	46243459	28460	25174	yes	0.0647482798611111	0.0306306402777778	0.0476894600694444	2	chr9_chr9_0	chr9_chr9_2	chr9_chr9_0_2
# 1405 (46 fields)
# 16 (47 fields)

# output file = Mus_musculus.GRCm38.cdna.realparalogs.gnid.expkallexpr.1.2.expkallavgexpr.expkallnormdeltaexpr.tsv
# gene1	gene2	ageclass	expexpr1	kallexpr1	expexpr2	kallexpr2	expavgexpr	kallavgexpr	deltaexpexpr	deltakallexpr
# ENSMUSG00000032079	ENSMUSG00000032080	3.veryold	27.4888	18.8386	33.748	15.0234	30.6184	16.931	0.204419	0.225325
# 1421 (11 fields)


BEGIN{
    OFS="\t"; 
    while (getline < fileRef >0)
    {
        split($1,a,"."); 
        expect[a[1]]=$7; 
        kall[a[1]]=$8; 
        agecl[a[1]]=$4;
    }
}  


NR==1{
    print "gene1", "gene2", "ageclass", "expexpr1", "kallexpr1", "expexpr2", "kallexpr2", "expavgexpr", "kallavgexpr", "deltaexpexpr", "deltakallexpr";
} 

NR>=2{
    split($0,a,"\t"); 
    s1=expect[a[1]]; 
    s2=expect[a[2]]; 
    savg=(s1+s2)/2; 
    k1=kall[a[1]]; 
    k2=kall[a[2]]; 
    kavg=(k1+k2)/2; 
    print a[1], a[2], agecl[a[1]], s1, k1, s2, k2, savg, kavg, (abs(s1-s2))/(savg+0.001), (abs(k1-k2))/(kavg+0.001);
} 

function abs(x){
    return (x>=0 ? x : ((-1)*x));
}