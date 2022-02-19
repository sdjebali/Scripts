#!/bin/bash
set -Eexo pipefail


# find.polya.in.genome.sh
# script that looks for occurences of polyA sites in a genome 
# - inputs:
#   * genome file (fasta format)
# - output:
#   * file of polyA site positions in genome (bed format) in the same directory as the genome file
# - dependences:
#   * homer

# example
#########
# srun --mem=8G --pty bash
# cd ~/fragencode/workspace/sdjebali/geneswitch/analysis/annotation.quality/tts
# module load bioinfo/HOMER-v4.10.4
# genome=~/fragencode/data/species/gallus_gallus/GRCg6a/gallus_gallus.fa
# pgm=~/fragencode/tools/multi/Scripts/find.polya.in.genome.sh
# time $pgm $genome 2> find.polya.in.genome.err
# MT	3	9	.	.	-
# MT	401	407	.	.	-
# 3252081 (6 fields) *** real	7m4.661s


# Check inputs are provided
###########################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: find.polya.in.genome.sh genome.fa >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a genome file in fasta format" >&2
    echo "" >&2
    echo "produces as output:" >&2
    echo "- the positions of the polyA sites in the genome (both strands) in bed format using homer" >&2
    echo "  and in the same directory as the genome file was" >&2
    echo "" >&2
    echo "Note: Needs homer in your path (tested with version 4.10.4 on genotoul slurm cluster)" >&2
    exit 1
fi

# Variable assignment
#####################
genome=$1
genbasetmp=${genome%.fa}
genbase=${genbasetmp%.fasta}

# Start of the script
#####################
# 1. Makes a profile of the searched motifs for polyA sites
###########################################################
echo "I am making a profile of the ATTAAA motif" >&2
seq2profile.pl ATTAAA 0 polyAhex > polyAhex.motif.profile.txt
echo "done" >&2
# >ATTAAA	polyAhex	8.28973911259755
# 0.997	0.001	0.001	0.001
# 0.001	0.001	0.001	0.997
# 0.001	0.001	0.001	0.997
# 0.997	0.001	0.001	0.001
# 0.997	0.001	0.001	0.001
# 0.997	0.001	0.001	0.001

# 2. Makes a profile of the reverse complement of the motif
###########################################################
echo "I am making a profile of AATAAA motif" >&2
seq2profile.pl AATAAA 0 polyAhex >> polyAhex.motif.profile.txt
echo "done" >&2

# 3. Scans the genome for the motif and its reverse complement keeping all occurences
#####################################################################################
echo "I am scanning the genome for polyA sites using the two provided motif (ATTAAA and AATAAA) profiles" >&2
scanMotifGenomeWide.pl polyAhex.motif.profile.txt $genome -keepAll -bed | awk '{print $1,$(NF-4)-1,$(NF-3),".",".",$NF}' OFS='\t'  > $genbase.polyAsites.bed
echo "done" >&2

# 4. Cleans
###########
echo "I am cleaning" >&2
# rm polyAhex.motif.profile.txt
echo "done" >&2
