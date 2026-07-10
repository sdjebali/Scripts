# join_two_clinical_files.awk

# Given:
# - a fileRef tsv file with header with all initial clinical data for all patients
#   and where we know that cid2.mean.wstcir, cid3.mean.wstcir, cid2.fast.adipo, 
#   cid2.fast.crp, cid3.fast.adipo, cid3.fast.crp, cid2.matsuda, cid3.matsuda are 
#   in columns 9, 14, 27, 28, 33, 34, 42 and 43 respectively (and never take the null value)
# - a main input tsv file with header with score2 information for all patients
# Provides:
# - a tsv file with header that has the anthropometric information, the non score2 info 
#   we are interested in (wstcir, matsuda, crp, adipo at cid2 and 3) 
#   and the score2 at cid2, cid3 and the normalized difference between cid2 and cid3

# example
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# basedir=~/bridge/workspace/sdjebali/papers/diogenes_paper_2026
# cd $basedir/results/clinical
# pgm=$basedir/scripts/clinical/join_two_clinical_files.awk
# time awk -v fileRef=$basedir/data/622.subjects.phenotypes.okheader.tsv 
# -f $pgm subjid.center.sex.diet.age.sc2.fm.cid123.normdiff23.tsv 
# > subjid.center.sex.diet.age.sc2.fm.wstcir.matsuda.crp.adipo.cid2.3.normdiff23.tsv

# fileRef=$basedir/data/622.subjects.phenotypes.okheader.tsv 
# SUBJID	Age	cid1.bmi	cid1.mean.wstcir	cid1.sbp	cid1.dbp	wk8.total.loss.lcd	cid2.bmi	cid2.mean.wstcir	cid2.sbp	cid2.dbp	diet	cid3.bmi	cid3.mean.wstcir	cid3.sbp	cid3.dbp	cid1.fast.chol	cid1.fast.tg	cid1.fast.hdlc	cid1.fast.fructos	cid1.fast.adipo	cid1.fast.crp	cid2.fast.chol	cid2.fast.tg	cid2.fast.hdlc	cid2.fast.fructos	cid2.fast.adipo	cid2.fast.crp	cid3.fast.chol	cid3.fast.tg	cid3.fast.hdlc	cid3.fast.fructos	cid3.fast.adipo	cid3.fast.crp	cid1.fast.ldlc	cid2.fast.ldlccid3.fast.ldlc	scr.total.eintake	wk4.total.eintake	cid3.total.eintake	cid1.matsuda	cid2.matsuda	cid3.matsuda	cid1.ffm	cid1.fm	cid2.ffm	cid2.fm	cid3.ffm	cid3.fm	cid1.fast.insulin	cid1.fast.glucose	cid2.fast.insulin	cid2.fast.glucose	cid3.fast.insulin	cid3.fast.glucose	cid1.insulin.auc	cid2.insulin.auc	cid3.insulin.auc	cid1.glucose.auc	cid2.glucose.auc	cid3.glucose.auc	cid1.homa.ir	cid2.homa.ir	cid3.homa.ir	cid1.mets	cid2.mets	cid3.mets
# 1-001-2	42	30.6794419649824	110.75	120.5	81	10.69999981	27.3399706625886	105.45	103	75	3	26.0291501513686	104.3	122.5	71.5	2.884000063	0.860000014	0.629999995	132	9.880000114	0.418000013	4.494999886	0.910000026	1.24000001	182	14.31461239	3.904999971	4.801000118	1.159999967	1.149999976	235	10.29192162	2.164999962	1.870000005	2.849999905	3.130000114	10385.09961	8351.099609	7903.600098	5.255000114	NA	8.37774086	58.89999771	40.1	55.19999695	37	54.29999924	34.90000153	11.85000038	4.400000095	6.420000076	NA	7.279999733	4.199999809	5799.75	NA	4159.200195	711	NA	666	2.32	NA	1.36	0	0	0
# 623 (67 fields)

# subjid.center.sex.diet.age.sc2.fm.cid123.normdiff23.tsv
# subjid	center	sex	diet	age	cid1.sc2	cid2.sc2	cid3.sc2	cid1.fm	cid2.fm	cid3.fm	norm.sc2.diff.23	norm.fm.diff.23
# 5-018-2	5	2	1	44	1.151870311886427	0.9046712451030947	1.1306941983805552	38.4	34.5	33.40000153	0.24984	-0.031884
# 604 (13 fields)

# subjid.center.sex.diet.age.sc2.fm.wstcir.matsuda.crp.adipo.cid2.3.normdiff23.tsv
# subjid	center	sex	diet	age	cid2.sc2	cid3.sc2	norm.sc2.diff.23	cid2.fm	cid3.fm	norm.fm.diff.23	cid2.wstcir	cid3.wstcir	norm.wstcir.diff.23	cid2.matsuda	cid3.matsuda	norm.matsuda.diff.23	cid2.crp	cid3.crp	norm.crp.diff.23	cid2.adipo	cid3.adipo	norm.adipo.diff.23
# 5-018-2	5	2	1	44	0.9046712451030947	1.1306941983805552	0.24984	34.5	33.40000153	-0.031884	85.85	86.25	0.00465929	11.35015392	12.32283783	0.0856979	0.949999988	0.470999986	-0.504211	9.098378181	12.54714584	0.379053
# 216 (23 fields) 


BEGIN{
    OFS="\t"; 
    while (getline < fileRef >0)
    {
        wstcir2[$1]=$9; 
        wstcir3[$1]=$14; 
        matsuda2[$1]=$42; 
        matsuda3[$1]=$43; 
        crp2[$1]=$28; 
        crp3[$1]=$34; 
        adipo2[$1]=$27; 
        adipo3[$1]=$33;
    }
} 

NR==1{
    print $1, $2, $3, $4, $5, $7, $8, $12, $10, $11, $13, "cid2.wstcir", "cid3.wstcir", "norm.wstcir.diff.23", "cid2.matsuda", "cid3.matsuda", "norm.matsuda.diff.23", "cid2.crp", "cid3.crp", "norm.crp.diff.23", "cid2.adipo", "cid3.adipo", "norm.adipo.diff.23";
} 

NR>=2&&$7!="NA"&&$8!="NA"&&$10!="NA"&&$11!="NA"&&wstcir2[$1]!="NA"&&wstcir3[$1]!="NA"&&matsuda2[$1]!="NA"&&matsuda3[$1]!="NA"&&crp2[$1]!="NA"&&crp3[$1]!="NA"&&adipo2[$1]!="NA"&&adipo3[$1]!="NA"{
    print $1, $2, $3, $4, $5, $7, $8, $12, $10, $11, $13, wstcir2[$1], wstcir3[$1], (wstcir3[$1]-wstcir2[$1])/wstcir2[$1], matsuda2[$1], matsuda3[$1], (matsuda3[$1]-matsuda2[$1])/matsuda2[$1], crp2[$1], crp3[$1], (crp3[$1]-crp2[$1])/crp2[$1], adipo2[$1], adipo3[$1], (adipo3[$1]-adipo2[$1])/adipo2[$1];
}