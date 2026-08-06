# make_2_hichip_partitions.awk

# Given a hichip bedpe file with 24 fields of which the last 10 correspond to categories of overlap of the two fragments
# with ccre pls, els and dels, provides the same file but with 2 additional columns for 2 partitions of the hichip relations:
# - one based on one frag overlap with pls, els and dels
# - one based on the two frag overlap with pls, els and dels

# In the input file we have the last 10 columns, from col 15 to 24 that represent the following information:
############################################################################################################
# - col 15 (a) = p (one frag over pls) / col 16 (b) = e (one frag over els) / col 17 (c) = de (one frag over dels)
# - col 18 (d) = pe (one frag over pls and one frag over els) / col 19 (e) = pde (...) / col 20 (f) = pp
# - col 21 (g) = ee / col 22 (h) = dede / col 23 (i) = na (one frag overlap nothing and the other pls or els or dels) / col 24 (j) = nana

# The 2 paritions will be made using these 2 priority orders:
#############################################################
# - Hichip relation partitioning based on single fragment overlap:
#   p (a) > de (c) > e (b) > na  (4 classes)
# - Hichip relation partitioning based on both fragments overlap:
#   pde (e) > pp (f) > pe (d) > dede (h) > ee (g) > na (i) > nana (j)  (7 classes)

# Example of usage:
###################
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/bridge/results/hichip/homo_sapiens/hg19/wilson_2020
# pgm=~/fragencode/tools/multi/Scripts/make_2_hichip_partitions.awk
# time awk -f $pgm hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.10catwrtplsels.bedpe > hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.10catwrtplsels.2partitions.bedpe2
# real	0m0.752s

# input file = hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.10catwrtplsels.bedpe 
# chr1	1250000	1255000	chr1	1310000	1315000	chr1:1250000-1255000?chr1:1310000-1315000	7.139784871	0	1	1	1	1	0	1	1	1	1	1	0	10	0	0
# 61136 (24 fields)

# output file = hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.10catwrtplsels.2partitions.bedpe
# chr1	1250000	1255000	chr1	1310000	1315000	chr1:1250000-1255000?chr1:1310000-1315000	7.139784871	0	1	1	1	1	01	1	1	1	1	0	1	0	0	0	p	pde
# 61136 (26 fields)

BEGIN{
    OFS="\t";
}

{
    print $0, onefragpart($15,$16,$17), twofragpart($18,$19,$20,$21,$22,$23,$24);
}

function onefragpart(a,b,c)
{
    return (a==1 ? "p" : (c==1 ? "de" : (b==1 ? "e" : "na")));
}

function twofragpart(d,e,f,g,h,i,j)
{
    return (e==1 ? "pde" : (f==1 ? "pp" : (d==1 ? "pe" : (h==1 ? "dede" : (g==1 ? "ee" : (i==1 ? "na" : (j==1 ? "nana" : "NA")))))));
}
