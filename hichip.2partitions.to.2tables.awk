# hichip.2partitions.to.2tables.awk
# given a bedpe file of hichip relations and the number of a column in this file representing a partition
# (based on one or the two frag), outputs a tsv file of 4 fields without header
# that has for each class of the partition
# - the total nb of relations
# - the class
# - the number of relations in this class
# - the corresponding %

# Example of usage
##################
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/bridge/results/hichip/homo_sapiens/hg19/wilson_2020
# pgm=~/fragencode/tools/multi/Scripts/hichip.2partitions.to.2tables.awk
# col=26
# awk -v col=$col -f $pgm hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.10catwrtplsels.2partitions.bedpe | sort -k2,2 

# input file = hPSC-CM_HiChIP_H3K27ac_combined_FitHiCq7_5kb_WashU.overpls.els.dels.first.second.part.10catwrtplsels.2partitions.bedp
# chr1	1250000	1255000	chr1	1310000	1315000	chr1:1250000-1255000?chr1:1310000-1315000	7.139784871	0	1	1	1	1	01	1	1	1	1	0	1	0	0	0	p	pde
# 61136 (26 fields)

# output 
# 61136	1_pde	30100	49.2345
# 61136	2_pp	1901	3.10946
# 61136	3_pe	979	1.60135
# 61136	4_dede	19151	31.3252
# 61136	5_ee	1027	1.67986
# 61136	6_na	7147	11.6903
# 61136	7_nana	831	1.35926
# 11 (4 fields)

BEGIN{
    cl["p"]="1_p";
    cl["de"]="2_de";
    cl["e"]="3_e";
    cl["na"]=(col==25 ? "4_na" : "6_na");
    cl["pde"]="1_pde";
    cl["pp"]="2_pp";
    cl["pe"]="3_pe";
    cl["dede"]="4_dede";
    cl["ee"]="5_ee";
    cl["nana"]="7_nana";
}

{
    tot++;
    tot1[cl[$col]]++;
}

END{
    OFS="\t";
    for(c in tot1)
    {
	print tot, c, tot1[c], tot1[c]/tot*100;
    }
}
