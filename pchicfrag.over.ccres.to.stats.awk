# pchicfrag.over.ccres.to.stats.awk

# this script takes as input:
#############################
# - ftype with -v which is the type of pchic fragment
# - the 4 column tsv file which is the output of an intersectBed between the pchic fragments in bed format
#   and ccres of the V3 catalog in bed format (V3) done with -wao and gathering results within the same bed file
#   as provided as input to intersectBed but with a 4th column with the list of ccres overlapping the fragment
#   (with indication of their type PLS, ELS, dELS..., see below)
# and outputs:
##############
# - a tsv file of stats on one row with indication of the type of fragment and then
#   * their number in the file
#   * the number of them that overlap a pls and the corresponding %
#   * the number of them that overlap an els and the corresponding %
#   * the number of them that overlap both a pls and an els and the corresponding %
#   * the number of them that overlap a distal els (dels) and the corresponding %
#   * the number of them that overlap both a pls and a dels and the corresponding %

# example of usage
##################
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/bridge/results/pchic/homo_sapiens/hg19/montefiori_2018
# pgm=~/fragencode/tools/multi/Scripts/pchicfrag.over.ccres.to.stats.awk
# ftype=promfrag
# fbeg=pchic.ipsccm.montefiori2018.hg19
# fend=extto5kb.sorted.over.plsels.bed
# time awk -v ftype=$ftype -f $pgm $fbeg.$ftype.$fend
# real	0m0.162s

# input file = $fbeg.$ftype.$fend = pchic.ipsccm.montefiori2018.hg19.promfrag.extto5kb.sorted.over.plsels.bed
# chr1	117110968	117115968	117110793:117111015:dELS,117111033:117111230:dELS_CTCF-bound,117111324:117111639:dELS,117111781:117111969:pELS,117111989:117112190:pELS,117112432:117112781:pELS_CTCF-bound,117112973:117113322:pELS_CTCF-bound,117113341:117113505:pELS_CTCF-bound,117113522:117113871:PLS_CTCF-bound,117113886:117114118:pELS_CTCF-bound,117114164:117114508:pELS,117114526:117114750:pELS,117114774:117115107:pELS,117115440:117115716:pELS_CTCF-bound,117115781:117115953:dELS_CTCF-bound,
# chr3	14579823	14584822	14579606:14579952:dELS,14580266:14580469:pELS_CTCF-bound,14580961:14581213:pELS_CTCF-bound,14581240:14581466:pELS_CTCF-bound,14581511:14581686:pELS_CTCF-bound,14582331:14582678:pELS,14582773:14582929:pELS,14583082:14583297:pELS,14581748:14582086:PLS_CTCF-bound,
# 16909 (4 fields)

# output
# promfrag	16909	15127	89.4612	16112	95.2865	14918	88.2252	7525	44.5029	6909	40.8599
# which means
# - 89% overlap a pls
# - 95% overlap an els
# - 88% overlap both
# - 44% overlap a dels
# - 40% overlap both a pls and a dels

{
    n++;
    if($4~/PLS/)
    {
	p++;
	if($4~/ELS/)
	{
	    e++;
	    pe++;
	    if($4~/dELS/)
	    {
		de++;
		pde++;
	    }
	}
    }
    else
    {
	if($4~/ELS/)
	{
	    e++;
	    if($4~/dELS/)
	    {
		de++;
	    }
	}
    }
}

END{
    OFS="\t";
    print ftype, n, p, p/n*100, e, e/n*100, pe, pe/n*100, de, de/n*100, pde, pde/n*100;
}
