#!/bin/bash
set -Eeuxo pipefail

# geom_density_by_factor.sh
# Make density plot for a column in a tsv file but splitting the data by a factor provided as another column in the tsv file
# - takes as input
#   * the input tsv file with header, from which the column to be plotted and the factor will be taken, the 1st one has to be numeric, the other one not
#   * the id in the header of the input file of the column we want to plot 
#   * the id in the header of the input file of the column in which the factor is
#   * a one word label for the x axis which is the meaning of the column in question
#   * a one word label for the title of the plot
#   * the min and max to zoom in for the x axis
# - produces as output and in the same directory as the input file:
#   * a png file (named after the one word label of the column and to the factor name), with the density plot of this column split according to the factor
# Note: needs reshape2 and ggplot2 libraries to be installed

# example
#########
# dir=~/regenet/workspace/sdjebali/egprediction/from3D/predictions/capturehic/jung.ren.2019
# pgm=~/fragencode/tools/multi/Scripts/geom_density_by_factor.sh
# module load system/R-4.0.2_gcc-7.2.0
# ct=HCmerge
# e=elt2A
# cd $dir/$ct/score2/pall
# time $pgm $ct.pall.score2.elt2.$e.uniq.lg.nbgn.elt2type.tsv elt2.lg nbgn $e.frag.length "$ct."$e".frag.length" 500-130000
# real	0m2.151s

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
    echo "Usage: geom_density_by_factor.sh input.tsv colid colfactor xaxislabel title xmin-xmax" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- the input tsv file with header, from which the column to be plotted and the factor will be taken, the first one with a numeric content, the second one not" >&2
    echo "- the id in the header of the input file of the numeric column we want to plot" >&2
    echo "- the id in the header of the input file of the column in which the factor is" >&2
    echo "- a one word label for the x axis which is the meaning of the column to be plotted" >&2
    echo "- a one word label for the title of the plot" >&2
    echo "- the x axis min and max to zoom in" >&2
    echo "" >&2
    echo "produces as output and in the same directory as the input:" >&2
    echo "- a png file (named after the one word label of the column and the factor name), with the density plot of this column split according to the factor" >&2
    echo "" >&2
    echo "Note: needs reshape2 and ggplot2 libraries to be installed" >&2
    exit 1
fi

# Set the output directory as the directory of the input tsv file and go there
##############################################################################
outdir=`dirname $1`
cd $outdir

# Set the min and max for the plot in case proper integers are given
#####################################################################
mM=$6
arg61=`echo $mM | awk '{split($1,a,"-"); print a[1]}'`
arg62=`echo $mM | awk '{split($1,a,"-"); print a[2]}'`
if [[ "$arg61" =~ ^[0-9]+$ ]] && [[ "$arg62" =~ ^[0-9]+$ ]] 
then
    m=$arg61
    M=$arg62
else
    echo "min and max values for the column to be plotted should be proper integers and be provided separated by a dash on the command line" >&2
    exit 1
fi

# Script content
#################
echo '
library(reshape2)
library(ggplot2)
theme_set(theme_bw(base_size = 16))
data = read.delim("'$1'", sep="\t", h=TRUE)
gp = ggplot(data, aes(x='$2', colour='$3')) + geom_density(size=1) + coord_cartesian(xlim = c('$m', '$M'))
gp = gp + labs(title = "'$5'", x="'$4'", y="Density") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=24))
w=8
h=5
ggsave(filename="'$4'_by_'$3'.density.png", h=h, w=w)
' | R --vanilla

