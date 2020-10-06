#!/bin/bash
set -Eeuxo pipefail

# geom_histogram_by_factor.sh (following model of geom_density_by_factor.sh)
# Make histogram plot for a numeric column of a tsv file but splitting the data by a factor provided as another column in the tsv file
# - takes as input
#   * the input tsv file with header, from which the numeric column and the factor will be taken
#   * the id of the column we want to plot in the header of the input file
#   * the id in the header of the input file of the column in which the factor is
#   * a one word label for the x axis which is the meaning of the column in question and that will give the name to the output file
#   * a string with minimum and maximum values for the column in question (min-max)
#   * a string to put as a title of the plot
#   * an optional number of bins (50 by default)
# - produces as output and in the same directory as the input:
#   * a png file (named after the one word label of the column and to the factor name), with the histogram plot of this column
# Note: needs R with ggplot2 library to be available

# example
#########
# dir=~/regenet/workspace/sdjebali/egprediction/from3D/predictions/capturehic/jung.ren.2019
# pgm=~/fragencode/tools/multi/Scripts/geom_histogram_by_factor.sh
# module load system/R-4.0.2_gcc-7.2.0
# ct=HCmerge
# e=elt2A
# cd $dir/$ct/score2/pall
# time $pgm $ct.pall.score2.elt2.$e.uniq.lg.nbgn.elt2type.tsv elt2.lg nbgn $e.frag.length 500-130000 $e.frag.length 
# real	0m6.727s

# input is like this
####################
# elt2.lg	nbgn
# 22752	2_morethan2gn
# 5087 (2 fields)

# Check all the obligatory inputs are indeed provided (should also check the input file exists and is not empty and that the ids are fine, for later)
#####################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ] || [ ! -n "$6" ]
then
   echo "" >&2
    echo "Usage: geom_histogram_by_factor.sh input.tsv colid colfactor xaxislabel min-max title <nbbins>" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- the input tsv file with header, from which the numeric column will be taken" >&2
    echo "- the id of the column we want to plot in the header of the input file" >&2
    echo "- the id in the header of the input file of the column in which the factor is" >&2
    echo "- a one word label for the x axis which is the meaning of the column in question" >&2
    echo "- a string with <min>-<max> values for the column to be plotted" >&2
    echo "- a string that will be the title for the plot" >&2
    echo "- an optional number of bins for the plot (default: 50)" >&2
    echo "" >&2
    echo "produces as output and in the same directory as the input:" >&2
    echo "- a png file (named after the one word label of the column and the factor name), with the histogram plot of this column" >&2
    echo "" >&2
    echo "Note: needs R and R ggplot2 library to be available" >&2
    exit 1
fi

# in case a bin number is provided and is indeed an integer
re='^[0-9]+$'
if [ ! -n "${7-default}" ] || ! [[ "${7-default}" =~ $re ]] 
then
    nbbins=50
else
    nbbins=$7
fi

# Set the output directory as the directory of the input tsv file and go there
##############################################################################
outdir=`dirname $1`
cd $outdir

# Set the min and max for the plot in case proper integers are given
#####################################################################
mM=$5
arg51=`echo $mM | awk '{split($1,a,"-"); print a[1]}'`
arg52=`echo $mM | awk '{split($1,a,"-"); print a[2]}'`
if [[ "$arg51" =~ ^[0-9]+$ ]] && [[ "$arg52" =~ ^[0-9]+$ ]] 
then
    m=$arg51
    M=$arg52
else
    echo "min and max values for the column to be plotted should be proper integers and be provided separated by a dash on the command line" >&2
    exit 1
fi

# Script content
#################
echo '
library(ggplot2)
theme_set(theme_bw(base_size = 16))
data = read.delim("'$1'", sep="\t", h=TRUE)
summary(data$'$2')
gp = ggplot(data=data, aes(x=data$'$2')) + geom_histogram(breaks=seq('$m', '$M', by='$nbbins'), col="red", fill="green", alpha = .2)
gp = gp + labs(title="'$6'", x="'$4'", y="Count")
gp = gp +  xlim(c('$m','$M')) + theme(plot.title = element_text(hjust = 0.5, size=24)) + theme(axis.text.x = element_text(angle = 90)) 
gp = gp + facet_grid(. ~'$3', scales = "free")
w=8
h=5
ggsave(filename="'$4'_by_'$3'.histogram.png", h=h, w=w)
' | R --vanilla

