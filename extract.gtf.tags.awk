# usage
# awk -v fld=transcript_id -f extract.gtf.tags.awk file.gtf > trid.txt
# note: will not work for a key that is the last one on the row unless there is a space at the end
# note: inspired from script bash with the same name but here can only extract one key

# example
# cd ~/fragencode/workspace/geneswitch/results/rnaseq/sus_scrofa/allbig/2020-08-20/temp/work/40/0bb79e5c41b1efbb1b81fe505ef2ee
# pgm=~/fragencode/workspace/sdjebali/geneswitch/pipelines/rnaseq/proj-gs-rna-seq/bin/extract.gtf.tags.awk
# awk -v fld=transcript_id -f $pgm sus_scrofa.gtf > trid.txt 

# input
# #!genome-build Sscrofa11.1
# #!genome-version Sscrofa11.1
# 5 (2 fields)
# 7872 (16 fields)
# 18008 (18 fields)
# 24510 (24 fields)
# 17465 (26 fields)
# 122614 (28 fields)
# 201259 (30 fields)
# 32 (32 fields)
# 900563 (34 fields)
# 190 (36 fields)
 

# output
# ENSSSCT00000065539
# ENSSSCT00000065539
# 1266633 (1 fields)


$1!~/#/{
    split($0,a,"\t");
    n=split(a[9],b,"; "); 
    for(i=1;i<=n;i++) 
    {
        split(b[i],c," "); 
        if(c[1]==fld) 
        {
            gsub(/\"/,"",c[2])
            print c[2];
        }       
    }
}
