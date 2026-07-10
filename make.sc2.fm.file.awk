# make.sc2.fm.file.awk

# Given
# - the complete file of clinical data for all patients
# - the file of basic info and score2 for all (patient,cid) pairs
# Produces 
# - a file with all patients with score2 and for which we also have 
#   fat mass at cid2 and cid3 and their normalised differences

# Example
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# basedir=~/bridge/workspace/sdjebali/papers/diogenes_paper_2026
# cd $basedir/clinical
# pgm=$basedir/scripts/clinical/make.sc2.fm.file.awk
# awk -v fileRef=$basedir/data/622.subjects.phenotypes.okheader.tsv 
# -f $pgm $basedir/data/subjid.cidno.center.sex.diet.age.score2.indok.tsv 
# > subjid.center.sex.diet.age.sc2.fm.cid123.normdiff23.tsv

# fileRef=$basedir/data/622.subjects.phenotypes.okheader.tsv 
# SUBJID	Age	cid1.bmi	cid1.mean.wstcir	cid1.sbp	cid1.dbp	wk8.total.loss.lcd	cid2.bmi	cid2.mean.wstcir	cid2.sbp	cid2.dbp	diet	cid3.bmi	cid3.mean.wstcir	cid3.sbp	cid3.dbp	cid1.fast.chol	cid1.fast.tg	cid1.fast.hdlc	cid1.fast.fructos	cid1.fast.adipo	cid1.fast.crp	cid2.fast.chol	cid2.fast.tg	cid2.fast.hdlc	cid2.fast.fructos	cid2.fast.adipo	cid2.fast.crp	cid3.fast.chol	cid3.fast.tg	cid3.fast.hdlc	cid3.fast.fructos	cid3.fast.adipo	cid3.fast.crp	cid1.fast.ldlc	cid2.fast.ldlc	cid3.fast.ldlc	scr.total.eintake	wk4.total.eintake	cid3.total.eintake	cid1.matsuda	cid2.matsuda	cid3.matsuda	cid1.ffm	cid1.fm	cid2.ffm	cid2.fm	cid3.ffm	cid3.fm	cid1.fast.insulin	cid1.fast.glucose	cid2.fast.insulin	cid2.fast.glucose	cid3.fast.insulin	cid3.fast.glucose	cid1.insulin.auc	cid2.insulin.auc	cid3.insulin.auc	cid1.glucose.auc	cid2.glucose.auc	cid3.glucose.auc	cid1.homa.ir	cid2.homa.ir	cid3.homa.ir	cid1.mets	cid2.mets	cid3.mets
# 1-001-2	42	30.6794419649824	110.75	120.5	81	10.69999981	27.3399706625886	105.45	103	75	3	26.0291501513686	104.3	122.5	71.5	2.884000063	0.860000014	0.629999995	132	9.880000114	0.418000013	4.494999886	0.910000026	1.24000001	182	14.31461239	3.904999971	4.801000118	1.159999967	1.149999976	235	10.29192162	2.164999962	1.870000005	2.849999905	3.130000114	10385.09961	8351.099609	7903.600098	5.255000114	NA	8.37774086	58.89999771	40.1	55.19999695	37	54.29999924	34.90000153	11.85000038	4.400000095	6.420000076	NA	7.279999733	4.199999809	5799.75	NA	4159.200195	711	NA	666	2.32	NA	1.36	00	0
# 623 (67 fields)

# $basedir/data/subjid.cidno.center.sex.diet.age.score2.indok.tsv 
# subjid	cidno	center	sex	diet	age	score2
# 1-001-2	1	1	2	3	42	1.095968450285778
# 1810 (7 fields)

# subjid.center.sex.diet.age.sc2.fm.cid123.normdiff23.tsv
# subjid	center	sex	diet	age	cid1.sc2	cid2.sc2	cid3.sc2	cid1.fm	cid2.fm	cid3.fm	norm.sc2.diff.23	norm.fm.diff.23
# 5-018-2	5	2	1	44	1.151870311886427	0.9046712451030947	1.1306941983805552	45.1	41.70000076	41.90000153	0.24984	-0.031884
# 604 (13 fields)  


BEGIN{
    OFS="\t"; 
    while (getline < fileRef >0)
    {
        fm[$1,1]=$45; 
        fm[$1,2]=$47; 
        fm[$1,3]=$49;
    }
} 

NR==1{
    print $1, $3, $4, $5, $6, "cid1.sc2", "cid2.sc2", "cid3.sc2", "cid1.fm", "cid2.fm", "cid3.fm", "norm.sc2.diff.23", "norm.fm.diff.23";
} 

NR>=2{
    patient[$1]++; 
    center[$1]=$3; 
    sex[$1]=$4; 
    diet[$1]=$5; 
    age[$1]=$6; 
    sc2[$1,$2]=$7;
} 

END{
    for(p in patient)
    {
        print p, center[p], sex[p], diet[p], age[p], sc2[p,1], sc2[p,2], sc2[p,3], fm[p,1], fm[p,2], fm[p,3], nn(sc2[p,2],sc2[p,3]), nn(fm[p,2],fm[p,3]);
    }
} 

function nn(x,y){
    return ((x!="NA"&&y!="NA") ? (y-x)/x : "NA");
}