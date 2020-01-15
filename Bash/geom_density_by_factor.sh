#!/bin/bash

# geom_density_by_factor.sh
# Make density plot for a column in a tsv file but splitting the data by a factor provided as another column in the tsv file
# - takes as input
#   * the input tsv file with header, from which the column to be plotted and the factor will be taken, the 1st one has to be numeric, the other one not
#   * the id in the header of the input file of the column we want to plot 
#   * the id in the header of the input file of the column in which the factor is
#   * a one word label for the x axis which is the meaning of the column in question
#   * a one word label for the title of the plot
# - produces as output and in the same directory as the input file:
#   * a png file (named after the one word label of the column and to the factor name), with the density plot of this column split according to the factor
# Note: needs reshape2 and ggplot2 libraries to be installed

# example
#########
# cd ~/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019/sumstats
# module load system/R-3.6.1
# pgm=~/fragencode/tools/multi/Scripts/Bash/geom_density_by_factor.sh
# time $pgm dist.score.quantiles.tsv distance score.quantile Distance
# real	0m39.706s

# input is like this
####################
# distance	score	quantile_score	quant_index_score	score.quantile
# 121316	1.65928294258351	(0.557,79.7]	4	4_(0.557,79.7]
# 6074050 (5 fields)


# Check all the obligatory inputs are indeed provided (should also check the input file exists and is not empty and that the ids are fine, for later)
#####################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ]
then
   echo "" >&2
    echo "Usage: geom_density_by_factor.sh input.tsv colid colfactor xaxislabel title" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- the input tsv file with header, from which the column to be plotted and the factor will be taken, the first one with a numeric content, the second one not" >&2
    echo "- the id in the header of the input file of the numeric column we want to plot" >&2
    echo "- the id in the header of the input file of the column in which the factor is" >&2
    echo "- a one word label for the x axis which is the meaning of the column to be plotted" >&2
    echo "- a one word label for the title of the plot" >&2
    echo "" >&2
    echo "produces as output and in the same directory as the input:" >&2
    echo "- a png file (named after the one word label of the column and to the factor name), with the density plot of this column split according to the factor" >&2
    echo "" >&2
    echo "Note: needs reshape2 and ggplot2 libraries to be installed" >&2
    exit 1
fi

# Set the output directory as the directory of the input tsv file and go there
##############################################################################
outdir=`dirname $1`
cd $outdir

# Script content
#################
echo '
library(reshape2)
library(ggplot2)
theme_set(theme_bw(base_size = 16))
data = read.delim("'$1'", sep="\t", h=TRUE)
gp = ggplot(data, aes(x='$2', colour='$3')) + geom_density(size=1) 
gp = gp + labs(title = "'$5'", x="'$4'", y="Density") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=24))
w=8
h=5
ggsave(filename="'$4'_by_'$3'.png", h=h, w=w)
' | R --vanilla

