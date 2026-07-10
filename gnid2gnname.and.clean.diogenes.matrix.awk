# gnid2gnname.and.clean.diogenes.matrix.awk
# Given
# - a fileRef1 tsv file without header that has the correspondance between long gene id (col1) and long gene name (col3)
# - a fileRef2 tsv file without header that has the correspondance between initial small gene probe id (col1) and smaller small gene probe id
# - a main tsv file with header that has the information we want
# Produces
# - a tsv file with header that has the same information as the input file except that
#  * long gene ids are replaced by long gene names (and gene id is kept when the gene name does not exist, 122 cases here)
#  * small probe ids are replaced by their shorter id 
#  * 1 and 2 in the sex columns are replaced by Men and Woman for later plotting purposes
#  * Delta are added for the clinical variables representing deltas

# Example
# cd $basedir/results
# pgm=$basedir/scripts/gathering/gnid2gnname.and.clean.diogenes.matrix.awk
# ouruniv=$basedir/data/block.spls.quant.allgn.id.bt.name.tsv
# small=$basedir/data/smallrna.id.simpleid.tsv
# time awk -v fileRef1=$ouruniv -v fileRef2=$small -f $pgm 
# gathering/66patients.baseclin.quantseq.mirna.fm.wstcir.matsuda.crp.adipo.normdiffcid23.sc2.cid2.cid3.delta.tsv 
# > 66patients.baseclin.quantseqgnname.mirnasimpleid.fm.wstcir.matsuda.crp.adipo.normdiffcid23.sc2.cid2.cid3.delta.tsv

# input file fileRef1=$ouruniv 
# ENSG00000000003	protein_coding	TSPAN6
# 13682 (3 fields)

# input file fileRef2=$small 
# hsa-let-7a-5p.MIMAT0000062_st	hsa-let-7a-5p
# hsa-let-7a-3p.MIMAT0004481_st	hsa-let-7a-3p
# 6609 (2 fields)

# main input file = gathering/66patients.baseclin.quantseq.mirna.fm.wstcir.matsuda.crp.adipo.normdiffcid23.sc2.cid2.cid3.delta.tsv 
# SUBJID center sex diet age ENSG00000000003 ... ENSG00000283633 hsa-let-7a-5p.MIMAT0000062_st ... gi:555853.gi555853_copy9_st norm.fm.diff.23 norm.wstcir.diff.23 norm.matsuda.diff.23 norm.crp.diff.23 norm.adipo.diff.23 cid2.sc2 cid3.sc2 norm.sc2.diff.23
# 2_090_1 2 1 1 35 -0.278904 ... -0.20481 0.0853377 ... -0.0132179 -0.105433879748406 -0.131299842187894 0.249915190783208 0.229811134231047 0.189048763669598 -0.00209233573493544 0.265692935698436 0.718210527020228
# 67 (20304 fields) *** 5 anthrop + 13,682 long genes + 6609 small genes + 8 clinical = 20304

# main output file = 66patients.baseclin.quantseqgnname.mirnasimpleid.fm.wstcir.matsuda.crp.adipo.normdiffcid23.sc2.cid2.cid3.delta.tsv
# SUBJID center sex diet age TSPAN6 ... ENSG00000283633 hsa-let-7a-5p ... gi:555853.10 Delta.FM Delta.WSTCIR Delta.MATSUDA Delta.CRP Delta.ADIPO CID2.SC2 CID3.SC2 Delta.SC2
# 2_090_1 2 Man 1 35 -0.278904 ... -0.20481 0.0853377 ... -0.0132179 -0.105433879748406 -0.131299842187894 0.249915190783208 0.229811134231047 0.189048763669598 -0.00209233573493544 0.265692935698436 0.718210527020228
# 67 (20304 fields)


BEGIN{
    OFS="\t";  
    while (getline < fileRef1 >0)
    {
        gnname[$1]=$3;
    } 
    while (getline < fileRef2 >0)
    {
        simple[$1]=$2;
    }
} 

# the header of the file
NR==1{
    # first anthropomorphic variables
    s=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t";
    # then the long genes 
    for(i=6; i<=13687; i++)
    {
        s=(s)(gnname[$i]!="NA" ? gnname[$i] : $i)("\t");
    } 
    # then the small genes
    for(i=13688; i<=20296; i++)
    {
        s=(s)(simple[$i])("\t");
    } 
    # finally the clinical variables
    for(i=20297; i<NF; i++)
    {
        if($i~/diff/)
        {
            split($i,a,"."); 
            $i="Delta."toupper(a[2]);
        }
        else
        {
            $i=toupper($i);
        } 
        s=(s)($i)("\t");
    } 
    if($i~/diff/)
    {
        split($i,a,"."); 
        $i="Delta."toupper(a[2]);
    }
    else
    {
        $i=toupper($i);
    } 
    print (s)($i);
} 

# the body of the file
NR>=2{
    $3=($3==1 ? "Men" : "Women"); 
    print;
} 