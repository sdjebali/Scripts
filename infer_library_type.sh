#!/bin/bash
set -Eexo pipefail

# Description
##############
# Takes as input a bam file with a set of mapped reads, a reference gene annotation and infers the sequencing library protocol (Unstranded, Mate2_sense & Mate1_sense) used to generate the rna-seq data. It does it by comparing the mapping strand in 1% of the aligments with the strand of the gene the read maps. Finally it produces three numbers: 

# 1) Fraction of reads explained by "1++,1--,2+-,2-+"
# 2) Fraction of reads explained by "1+-,1-+,2++,2--"
# 3) Fraction of reads explained by other combinations

# They give information regarding the library. They contain several strings of three characters, i.e. 1+-, where:
#   Character 1. 1 and 2 are mate1 and mate2 respectively.
#   Character 2. + and - is the strand where the read maps.
#   Character 3. + and - is the strand where the gene in which the read overlaps is annotated.

# You can apply the following rules to infer the used library from this information:

#    NONE. Not strand-specific protocol (unstranded data). Fraction of reads explained by “1++,1–,2+-,2-+” and “1+-,1-+,2++,2–” close to 0.5000 in both cases.

# Strand-specific protocols (stranded data):
#    MATE1_SENSE. Fraction of reads explained by “1++,1–,2+-,2-+” close to 1.0000.
#    MATE2_SENSE. Fraction of reads explained by “1+-,1-+,2++,2–” close to 1.0000.


# usage
#######
# infer_library_type.sh alignments.bam annotation.gff

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts
# - uses bedtools
# - comes from chimpipe but changed a bit on July 3rd 2020 to call awk script at the same level as the bash script
# - could be made faster by parallelizing the samtools part

# In case the user does not provide any input file, an error message is raised
##############################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
echo "" >&2
echo "infer_library_type.sh"
echo "" >&2
echo "Takes as input a bam file with a set of mapped reads, a reference gene annotation and infers the sequencing library protocol (Unstranded, Mate2_sense & Mate1_sense) used to generate the rna-seq data. It does it by comparing the mapping strand in 1% of the aligments with the strand of the gene the read maps."
echo "" >&2
echo Usage:  infer_library_type.sh alignments.bam annotation.gff >&2
echo "" >&2
echo "" >&2
exit 1
fi

# GETTING INPUT ARGUMENTS
#########################
bamfile=$1
annot=$2

# Directories 
#############
## Set root directory
path="`dirname \"$0\"`"              # relative path
rootdir="`( cd \"$path\" && pwd )`"  # absolute path

if [ -z "$rootdir" ] ; 
then
  # error; for some reason, the path is not accessible
  # to the script
  log "Path not accessible to the script\n" "ERROR" 
  exit 1  # fail
fi

# PROGRAMS
##########
cutgff=$rootdir/cutgff.awk

# START
########
## Comment: samtools view -F 260 ** filter out unmapped reads (4) + secondary alignments (254) = 260
## Note about random sampling with samtools view -s:
# Current versions of samtools (ab)use the integer part of the -s value to set the seed:
# -s FLOAT    Integer part is used to seed the random number generator. Part after the decimal point sets the fraction of templates/pairs to subsample
# So two runs with e.g. view -s 123.1 / view -s 456.1 should select two distinct randomly-selected 10% subsets.
# I will take a random number between 0 and 10 as a seed
randomSeed=`shuf -i1-10 -n1`

# read ids are either like this
# xxxx/1 or xxxx/2
# or
# xxxx 1:N:0:ATTACTCG+AGGCTATA or xxxx 2:N:0:ATTACTCG+AGGCTATA
# so in the 1st case the kind of read is after the / and in the 2nd case before the first : of the second string
# however here we use the bam file and not the fastq file and in the bam file we do not have the second part of the id
# and therefore we have to look at the bitwise flag
# so for the moment this script only works with /1 /2 read id nomenclature
# unless we know how to get 1 or 2 from the bitwise flag
# in fact from this bam file
# LH00392:213:22LLCVLT4:6:1101:35958:1056	419	chr1	8894222	110	150M	=	26475639	17581568	*	*	RG:Z:0	NH:i:3	NM:i:2	XT:A:R	md:Z:1G121A26
# LH00392:213:22LLCVLT4:6:1101:35958:1056	419	chr1	26475591	119	150M	=	26475639	199	*	*	RG:Z:0	NH:i:3	NM:i:1	XT:A:R	md:Z:1G148
# LH00392:213:22LLCVLT4:6:1101:35958:1056	339	chr1	26475639	119	118M1D32M	=	26475591	-199	*	*	RG:Z:0	NH:i:3	NM:i:2	XT:A:R	md:Z:1G30>1+118
# LH00392:213:22LLCVLT4:6:1101:35958:1056	339	chr1	26475639	110	118M1D32M	=	8894222	-17581568	*	*	RG:Z:0	NH:i:3	NM:i:2	XT:A:R	md:Z:1G30>1+118
# to know if read 1 or 2 we need to look at the bitwise flag in column 2
# in fact
# - for 2 first rows we have read 2
#   since bitwise flag = 419 = 00110100011
#   and therefore because in 8th pos from the right (2^7=128) we have 1 it means it is read2 (or because in 7th pos we have 0)
# - for 2 next rows we have read 1
#   since bitwise flag = 339 = 101010011
#   and therefore because in 8th pos from the right (2^7=128) we have 0 it means it is read 1 (or because in 7th pos we have 1)
# so we could use that in the script but it would mean converting a number in binary first


awk -v elt='exon' '$3==elt' "$annot" | awk -v to=8 -f $cutgff | sort -k1,1 -k4,4 -k5,5 | uniq | bedtools intersect -abam <(samtools view -b -s ${randomSeed}.01 -F 260 "$bamfile") -b stdin -split -bed -wo | awk '{print $4, $6, $19}' | uniq | awk '{readCount["total"]++; split($1,a,"/"); readCount[a[2]":"$2":"$3]++} END{fraction1=(readCount["1:+:+"]+readCount["1:-:-"]+readCount["2:+:-"]+readCount["2:-:+"]); fraction2=(readCount["1:+:-"]+readCount["1:-:+"]+readCount["2:+:+"]+readCount["2:-:-"]); other=(readCount["total"]-(fraction1+fraction2)); print (fraction1/readCount["total"]*100), (fraction2/readCount["total"]*100), other;}' 



