# split.pefastq.into.twofiles.awk
# takes as input a PE fastq file and splits it into two files, one for R1 and one for R2
# it also takes as input the name of the file to split, which is passed as a variable "filename" to the awk command

# example of usage
# awk -v filename=RandomTranscriptCoverage.for.art.pe75.fastq -f split.pefastq.into.twofiles.awk RandomTranscriptCoverage.for.art.pe75.fastq


BEGIN{
        split(filename,f,".fastq");
}

NR%4==1{
    split($0,a,"/"); 
    if(a[2]==1)
    {
        ok1=1; 
        ok2=0; 
        print $0 > f[1]".R1.fastq";
    }
    else
    {
        if(a[2]==2)
        {
            ok2=1; 
            ok1=0; 
            print $0 > f[1]".R2.fastq";
        }
    }
} 

NR%4!=1{
    if(ok1==1)
    {
        print $0 > f[1]".R1.fastq";
    } 
    else
    {
        if(ok2==1)
        {
            print $0 > f[1]".R2.fastq";
        }
    }
} 