# techrepcount2biorepcount.awk
# takes as input a matrix raw count for tech rep with header and a fileRef file with header that has biorep id
# and comma separated list of two tech ids as in matrix file in two first columns, and outputs a matrix raw count
# for bio rep by adding counts of the two tech reps

# example of usage
# cd ~/fragencode/workspace/sdjebali/geneswitch/wp5/atacseq/qc.measures/firstrun.techrep.as.biorep
# pgm=~/fragencode/tools/multi/Scripts/techrepcount2biorepcount.awk
# time awk -v fileRef=peakid.rawcount.eachtechrep.tsv -f $pgm biorepid.techreps.tiss.age.mid.mdiet.animal.id.name.gender.mbreed.fid.read1.2.nb.bp.trim.nb.pcent.gc.duplev.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.frip.tsv > peakid.rawcount.eachbiorep.tsv 
# real	9m55.478s

# peakid.rawcount.eachtechrep.tsv
# peakid	liver_fetus_lowfibre_R34	liver_fetus_highfibre_R45	...	liver_fetus_lowfibre_R24	liver_fetus_lowfibre_R13
# 1:21547:24249:+	276	484	...	1123	321
# 245453 (633 fields)


# biorepid.techreps.tiss.age.mid.mdiet.animal.id.name.gender.mbreed.fid.read1.2.nb.bp.trim.nb.pcent.gc.duplev.read12.clnreadnb.map.nb.pcent.mt.nb.pcent.fracdupl.frip.tsv
# biorep	techreps	tissue	age	motherid	motherdiet	animalid	animalname	gender	motherbreed	fatherid	read1.initnb	read2.initnb	read1.initbp.nb	read1.trimbp.nb	read1.trimbp.pcent	read2.initbp.nb	read2.trimbp.nb	read2.trimbp.pcent	read1.gcsummary	read2.gcsummary	read1.duplevsummary	read2.duplevsummary	clnread.tot	clnreadmap.nb	clnreadmap.pcent	clnreadmt.nb	clnreadmt.pcent	frac.dupl.picard	frip
# SSC_WU_WP5_Piglet_6447:skeletalmuscle	skelmus_piglet_highfibre_R38,skelmus_piglet_highfibre_R37,	skeletalmuscle	piglet	194	highfibre	6447	SSC_WU_WP5_Piglet_6447	female	TN60	E6021	72091288	72091288	10813693200	7031819374	65.027	10813693200	6992384383	64.6623	PASS	PASS	PASS	PASS	89828778	89318716	99.4322	2469796	2.74945	0.206196	0.226743
# 317 (30 fields)

# output peakid.rawcount.eachbiorep.tsv 
# peakid	SSC_WU_WP5_Piglet_6447:skeletalmuscle	SSC_WU_WP5_FT_70dpf_196.1:skeletalmuscle	...	SSC_WU_WP5_FT_70dpf_194.4:skeletalmuscle	SSC_WU_WP5_FT_70dpf_856.2:liver	SSC_WU_WP5_FT_70dpf_196.1:liver
# 1:21547:24249:+	495	482	...	314	408
# 245453 (317 fields)

BEGIN{
    OFS="\t";
    # read the matrix read count by tech rep
    while (getline < fileRef >0)
    {
	n++;
	# when in the header, remember the tech rep id corresponding to each column
	if(n==1)
	{
	    for(i=2; i<=NF; i++)
	    {
		tid[i]=$i;
	    }
	}
	# when in the body, remember the order of each peak as well as the raw read count of the peak for each tech rep id
	else
	{
	    for(i=2; i<=NF; i++)
	    {
		peakid[n-1]=$1;
		peakreadnb[$1,tid[i]]=$i;
	    }
	}
    }
}

# read the body of the metadata file that has biorep id followed by two tech rep id as a comma separated list
NR>=2{
    split($2,a,",");
    # remember the order of the biorep ids
    m++;
    biorepid[m]=$1;
    # for each biopep number, remember its two tech rep ids
    techrepid[m,1]=a[1];
    techrepid[m,2]=a[2];
}

END{
    # at the end, first write all the biorep ids in the order they were in the metadata file as a header
    s="peakid\t";
    for(k=1; k<=(m-1); k++)
    {
	s=(s)(biorepid[k])("\t");
    }
    print (s)(biorepid[k]);

    # and then write each peak in the order it came in the peak file, together with the raw read count of each biorep in order
    # and obtained as the sum of the raw read count of its two techrep
    for(l=1; l<=(n-1); l++)
    {
	s=(peakid[l])("\t");
	for(k=1; k<=(m-1); k++)
	{
	    s=(s)(peakreadnb[peakid[l],techrepid[k,1]]+peakreadnb[peakid[l],techrepid[k,2]])("\t");
	}
	print (s)(peakreadnb[peakid[l],techrepid[k,1]]+peakreadnb[peakid[l],techrepid[k,2]]);
    }
}
