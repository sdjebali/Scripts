# clinvar66to24.awk

# In the context of the diogenes project, and to use mixomics pls with multilevel option
# on rnaseq data and clinical data from the same 908 patient/cids as rows
# I need to make a clinical variable file from the 623 row * 67 column file (available for the
# 622 patients and the 67 clinvar including patient id) into a 909 row * 25 column file
# (for the 908 patient/cids common to rnaseq and clinical and the labExpId + age + diet + loss
# (the 3 I will repeat at each cid) and 21 normal var which are available at the 3 cids 

# srun --x11 --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/mixomics/Viguerie.Moro.Obesity/rnaseq.clinic
# pgm=~/fragencode/tools/multi/Scripts/clinvar66to24.awk
# time awk -v fileRef1=../clinical.phenotypes/clinvar.name.cid1.2.3.tsv -v fileRef2=../clinical.phenotypes/622.subjects.phenotypes.okheader.tsv -f $pgm rnaseq.clinic.samples.tsv > rnaseq.clinic.samples.24phenotypes.tsv

# fileRef1=../clinical.phenotypes/clinvar.name.cid1.2.3.tsv
# clinvar	name.at.cid1	name.at.cid2	name.at.cid3
# total.eintake	scr.total.eintake	wk4.total.eintake	cid3.total.eintake
# 22 (4 fields)

# fileRef2=../clinical.phenotypes/622.subjects.phenotypes.okheader.tsv
# SUBJID	Age	cid1.bmi	cid1.mean.wstcir	cid1.sbp	cid1.dbp	wk8.total.loss.lcd	cid2.bmi	cid2.mean.wstcir	cid2.sbp	cid2.dbp	diet	cid3.bmi	cid3.mean.wstcir	cid3.sbp	cid3.dbp	cid1.fast.chol	cid1.fast.tg	cid1.fast.hdlc	cid1.fast.fructos	cid1.fast.adipo	cid1.fast.crp	cid2.fast.chol	cid2.fast.tg	cid2.fast.hdlc	cid2.fast.fructos	cid2.fast.adipo	cid2.fast.crp	cid3.fast.chol	cid3.fast.tg	cid3.fast.hdlc	cid3.fast.fructos	cid3.fast.adipo	cid3.fast.crp	cid1.fast.ldlc	cid2.fast.ldlc	cid3.fast.ldlc	scr.total.eintake	wk4.total.eintake	cid3.total.eintake	cid1.matsuda	cid2.matsuda	cid3.matsuda	cid1.ffm	cid1.fm	cid2.ffm	cid2.fm	cid3.ffm	cid3.fm	cid1.fast.insulin	cid1.fast.glucose	cid2.fast.insulin	cid2.fast.glucose	cid3.fast.insulin	cid3.fast.glucose	cid1.insulin.auc	cid2.insulin.auc	cid3.insulin.auc	cid1.glucose.auc	cid2.glucose.auc	cid3.glucose.auc	cid1.homa.ir	cid2.homa.ir	cid3.homa.ir	cid1.mets	cid2.mets	cid3.mets
# 1-001-2	42	30.6794419649824	110.75	120.5	81	10.69999981	27.3399706625886	105.45	103	75	3	26.0291501513686	104.3	122.5	71.5	2.884000063	0.860000014	0.629999995	132	9.880000114	0.418000013	4.494999886	0.910000026	1.24000001	182	14.31461239	3.904999971	4.801000118	1.159999967	1.149999976	235	10.29192162	2.164999962	1.870000005	2.849999905	3.130000114	10385.09961	8351.099609	7903.600098	5.255000114	NA	8.37774086	58.89999771	40.1	55.19999695	37	54.29999924	34.90000153	11.85000038	4.400000095	6.420000076	NA	7.279999733	4.199999809	5799.75	NA	4159.200195	711	NA	666	2.32	NA	1.36	00	0
# 623 (67 fields)

# main input file rnaseq.clinic.samples.tsv
# 1_001_2_1
# 1_001_2_3
# 908 (1 fields)

# main output file rnaseq.clinic.samples.24phenotypes.tsv
# labExpId	age	diet	loss.lcd	total.eintake	fast.ldlc	fast.crp	fast.hdlc	glucose.auc	fast.insulin	matsuda	bmi	fm	mean.wstcir	fast.adipo	sbp	ffm	fast.tg	fast.chol	mets	fast.fructos	fast.glucose	dbp	homa.ir	insulin.auc
# 1-001-2-1	42	3	10.69999981	10385.09961	1.870000005	0.418000013	0.629999995	711	11.85000038	5.255000114	30.6794419649824	40.1	110.75	9.880000114	120.5	58.89999771	0.860000014	2.884000063	0	132	4.400000095	81	2.32	5799.75
# 909 (25 fields)

BEGIN{
    OFS="\t";
    # start writing the header of the output file with labExpId (= will be patient/cid), the age, the diet and the loss
    # which are the clin var that belong to the patient and that I will repeat at the 3 cids
    s="labExpId\tage\tdiet\tloss.lcd\t";
    # then for the other 21 variables of fileRef1, put them in the header one after the other
    # but remember in which order so that when I read the big clinical file, I know in which order I need to put them
    while (getline < fileRef1 >0)
    {
	n++;
	# in case we are after the header and before the last clinvar (just because we do not want a tab at the end)
	if(n>=2&&n<=21)
	{
	    s=(s)($1)("\t");
	    varname[n-1]=$1;
	    vname[n-1,1]=$2;
	    vname[n-1,2]=$3;
	    vname[n-1,3]=$4;
	}
	# in case we are at the last clinvar, we print the end of the header and remember the name of the clinvar we want at that point of the file
        # clinvar	name.at.cid1	name.at.cid2	name.at.cid3
        # total.eintake	scr.total.eintake	wk4.total.eintake	cid3.total.eintake
        # 22 (4 fields)
	if(n==22)
	{
	    print (s)($1);
	    varname[n-1]=$1;
	    vname[n-1,1]=$2;
	    vname[n-1,2]=$3;
	    vname[n-1,3]=$4;
	}
    }
    # when we read the big clinical file, we want to remember the value of each variable for each patient
    while (getline < fileRef2 >0)
    {
	m++;
	# if we are in the header, we just remember the variable at each column and the index of each variable
	if(m==1)
	{
	    for(i=3; i<=NF; i++)
	    {
		var[i]=$i;
	    }
	}
	# if we are in the body of the patient file, we remember the age, diet and loss of the patient first
	# and then the 3 clinvar of each patient that can be in any order
	else
	{
	    # here $1 is a patient id
	    age[$1]=$2;
	    diet[$1]=$12;
	    loss[$1]=$7;
	    for(i=3; i<=NF; i++)
	    {
		val[$1,var[i]]=$i;
	    }
	}
    }
}

{
    # when we read the main input file, each row is a patient/cid, so need to remove the last bit
    # to have the patient id (a[1]"-"a[2]"-"a[3])
    gsub(/\_/,"-",$1);
    split($1,a,"-");
    # first we put the patient/cid, then the age of the patient, its diet and its loss (will be present
    # at most 3 times, if the patient is present at the 3 cids
    ind=a[1]"-"a[2]"-"a[3];
    cidno=a[4];
    s=$1"\t"age[ind]"\t"diet[ind]"\t"loss[ind]"\t";
    # for the other 21 variables present in fileRef1, we will write in s their value for the patient and cid in question
    for(i=1; i<=(n-2); i++)
    {
	s=(s)(val[ind,vname[i,cidno]])("\t");
    }
    # at the end we print s with the last clinvar
    print (s)(val[ind,vname[i,cidno]]);
}

