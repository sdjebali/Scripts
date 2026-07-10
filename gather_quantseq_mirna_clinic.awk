# gather_quantseq_mirna_clinic.awk
# Given:
# - a fileRef0 tsv file with header of unscaled anthropomorphic and clinical data about patients whose ids are in column 1 (- in id)
# - a fileRef1 tsv file with header of patients (col1, _ in id) with normalised long gene expression difference between cid2 and cid3 ((cid3-cid2)/(cid2+1e-3))
# - a fileRef2 tsv file with header of patients (col1, _ in id) with normalised small gene expression difference between cid2 and cid3 ((cid3-cid2)/cid2)
# - a main input file with the scaled clinical data for patients whose id is with - and in column 1
# Produces:
# -  a tsv file with header that gathers unscaled anthropomorphic information, norm delta expr of long, norm delta expr of small and scaled clinical data
#    but just for patients with no NA in those (also homogeneises the patients' ids since with - in clinical and with _ in expression data)

# example
# basedir=~/bridge/workspace/sdjebali/papers/diogenes_paper_2026
# cd $basedir/results/gathering
# pgm=$basedir/scripts/gathering/gather_quantseq_mirna_clinic.awk
# time awk -v fileRef0=$basedir/results/clinical/subjid.center.sex.diet.age.sc2.fm.wstcir.matsuda.crp.adipo.cid2.3.normdiff23.tsv 
# -v fileRef1=$basedir/results/long/patients_with_quantseqexpr_cid23_normgnexprdiff23.tsv 
# -v fileRef2=$basedir/results/small/patients_with_mirnaexpr_cid23_normgnexprdiff23.tsv -f $pgm 
# $basedir/results/clinical/215.subjects.cid2.3.23.phenotypes.centeredscaled.tsv 
# > 66patients.baseclin.quantseq.mirna.fm.wstcir.matsuda.crp.adipo.normdiffcid23.sc2.cid2.cid3.delta.tsv

# input file fileRef0 = $basedir/results/clinical/subjid.center.sex.diet.age.sc2.fm.wstcir.matsuda.crp.adipo.cid2.3.normdiff23.tsv 
# subjid	center	sex	diet	age	cid2.sc2	cid3.sc2	norm.sc2.diff.23	cid2.fm	cid3.fm	norm.fm.diff.23	cid2.wstcir	cid3.wstcir	norm.wstcir.diff.23	cid2.matsuda	cid3.matsuda	norm.matsuda.diff.23	cid2.crp	cid3.crp	norm.crp.diff.23	cid2.adipo	cid3.adipo	norm.adipo.diff.23
# 5-018-2	5	2	1	44	0.9046712451030947	1.1306941983805552	0.24984	34.5	33.40000153	-0.031884	85.85	86.25	0.00465929	11.35015392	12.32283783	0.0856979	0.949999988	0.470999986	-0.504211	9.098378181	12.54714584	0.379053
# 216 (23 fields)

# input file fileRef1 = $basedir/results/long/patients_with_quantseqexpr_cid23_normgnexprdiff23.tsv 
# patient	ENSG00000000003	ENSG00000000005	ENSG00000000419	...	ENSG00000283078	ENSG00000283632	ENSG00000283633
# 2_001_2	0.300355	-0.411881	0.309421	...	NA	2.12952	-0.799337
# 144 (13683 fields)

# input file fileRef2 = $basedir/results/small/patients_with_mirnaexpr_cid23_normgnexprdiff23.tsv 
# patient	hsa-let-7a-5p.MIMAT0000062_st	hsa-let-7a-3p.MIMAT0004481_st	hsa-let-7a-2-3p.MIMAT0010195_st	...	gi:555853.gi555853_copy7_st	gi:555853.gi555853_copy8_st	gi:555853.gi555853_copy9_st
# 2_001_2	0.119966	0.00173605	-0.0171799	...	-0.0097605	-0.0524913	-0.0904286
# 164 (6610 fields)

# main input file = $basedir/results/clinical/215.subjects.cid2.3.23.phenotypes.centeredscaled.tsv 
# SUBJID	center	sex	diet	age	cid2.sc2	cid3.sc2	norm.sc2.diff.23	cid2.fm	cid3.fm	norm.fm.diff.23	cid2.wstcir	cid3.wstcir	norm.wstcir.diff.23	cid2.matsuda	cid3.matsuda	norm.matsuda.diff.23	cid2.crp	cid3.crp	norm.crp.diff.23	cid2.adipo	cid3.adipo	norm.adipo.diff.23
# 5-018-2	0.262457143888782	0.737712347812099	-1.49870274429222	0.317715343088411	-0.531108172816617	-0.415469901153674	0.693037740366694	-0.20621546022137	-0.221290223789479	-0.094325032741949	-1.09405198982156	-1.03494160476201	0.0642783206057395	1.1031069033284	1.4750301076408	0.233850980752308	-0.716646859325442	-0.871354249622079	-0.454630672558315	0.0434144458381082	0.486907426321029	0.167906168851558
# 216 (23 fields) 

# main output file = 66patients.baseclin.quantseq.mirna.fm.wstcir.matsuda.crp.adipo.normdiffcid23.sc2.cid2.cid3.delta.tsv
# SUBJID center sex diet age ENSG00000000003 ... ENSG00000283633 hsa-let-7a-5p.MIMAT0000062_st ... gi:555853.gi555853_copy9_st norm.fm.diff.23 norm.wstcir.diff.23 norm.matsuda.diff.23 norm.crp.diff.23 norm.adipo.diff.23 cid2.sc2 cid3.sc2 norm.sc2.diff.23
# 2_090_1 2 1 1 35 -0.278904 ... -0.20481 0.0853377 ... -0.0132179 -0.105433879748406 -0.131299842187894 0.249915190783208 0.229811134231047 0.189048763669598 -0.00209233573493544 0.265692935698436 0.718210527020228
# 67 (20304 fields) *** 5 anthrop + 13,682 long genes + 6609 small genes + 8 clinical = 20304


BEGIN{
    OFS="\t";  
    while (getline < fileRef0 >0)
    {
        gsub(/-/,"_",$1); 
        center[$1]=$2; 
        sex[$1]=$3; 
        diet[$1]=$4; 
        age[$1]=$5;
    } 
    n=0; 
    while (getline < fileRef1 >0)
    {
        n++; 
        if(n==1)
        {
            for(i=2; i<=NF; i++)
            {
                sh1=(sh1)($i)("\t");
            }
        }
        else
        {
            for(i=2; i<=NF; i++)
            {
                s1[$1]=(s1[$1])($i)("\t");
            }
        }
    } 
    n=0; 
    while (getline < fileRef2 >0)
    {
        n++; 
        if(n==1)
        {
            for(i=2; i<NF; i++)
            {
                sh2=(sh2)($i)("\t");
            } 
            sh2=(sh2)($i);
        }
        else
        {
            for(i=2; i<NF; i++)
            {
                s2[$1]=(s2[$1])($i)("\t");
            } 
            s2[$1]=(s2[$1])($i);
        }
    }
} 

NR==1{
    s=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"; 
    s=(s)(sh1)(sh2)("\t"); 
    print (s)($11)("\t")($14)("\t")($17)("\t")($20)("\t")($23)("\t")($6)("\t")($7)("\t")($8);
} 

NR>=2{
    gsub(/-/,"_",$1); 
    if(s1[$1]!=""&&s2[$1]!="")
    {
        s=($1)("\t")(center[$1])("\t")(sex[$1])("\t")(diet[$1])("\t")(age[$1])("\t"); 
        s=(s)(s1[$1])(s2[$1])("\t"); 
        print (s)($11"\t"$14"\t"$17"\t"$20"\t"$23"\t"$6"\t"$7"\t"$8);
    }
} 

function idem(x,y){
    return (x==y ? 1 : 0);
}