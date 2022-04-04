#!/bin/bash
set -Eexo pipefail

# boxplots.sh
# make boxplots for values (y) belonging to several categories (x) from an input tsv file with ggplot2
# !!! this bash script is only valid for plotting 4 boxplots corresponding to 4 transcript sets allref, exprref, allnew, exprnew !!!
# - takes as input
#   * absolute path to input tsv file
#   * header key for x values in the tsv file
#   * header key for y values in the tsv file
#   * label for y axis
#   * min y value to be plotted
#   * max y value to be plotted
#   * absolute path to output file (with extension), could be pdf, png, eps.
# Note1: needs ggplot2 libraries to be available
# Note2: in order to have the boxplots in a given order one can number the different sets for which we want boxplots

# example
# cd ~/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Stringtie/WithoutAnnot/Test
# input=~/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Plots/ExonPerTranscript/annotation_nbexintr_forggplot.tsv
# output=annotation_nbexintr_forggplot.png
# time boxplots.sh $input annotation nb_ex_in_transcript "Number of exons per transcript" 20 $output
# real    0m12.362s

# input is like this
####################
# annotation    nb_ex_in_transcript
# 1_Gencv19_all   2
# 1220427 (2 fields)

if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ] || [ ! -n "$6" ] || [ ! -n "$7" ]
then
   echo "" >&2
    echo Usage: boxplots.sh input.tsv headkey_xval headkey_yval yaxis_label min_yval max_yval output_file >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- absolute path to input tsv file" >&2
    echo "- header key for x values in the tsv file" >&2
    echo "- header key for y values in the tsv file" >&2
    echo "- label for y axis" >&2
    echo "- min y value to be plotted" >&2
    echo "- max y value to be plotted" >&2
    echo "- absolute path to output file (with extension), could be pdf, png, eps." >&2
    echo "" >&2
    echo "produces as output:" >&2
    echo "- the file specified as the last argument with as many boxplots as categories for x and with the values indicated in y" >&2
    echo "" >&2
    echo "WARNING: this bash script is only valid for plotting 4 boxplots corresponding to 4 transcript sets allref, exprref, allnew, exprnew" >&2
    echo "WARNING: this bash script will plot the boxplots in the alphabetical order of the trpred ids" >&2
    echo "Note1: needs ggplot2 libraries to be installed" >&2
    echo "Note2: in order to have the boxplots in a given order one can number the different sets for which we want boxplots" >&2
    exit 1
fi

echo '
library(ggplot2)
sessionInfo()
theme_set(theme_bw(base_size = 16))
data = read.delim("'$1'", sep="\t", h=TRUE)
gp = ggplot(data) + geom_boxplot(aes(y='$3',x=factor('$2'),fill=factor('$2')), varwidth = TRUE, notch=T)
gp = gp + scale_x_discrete(labels=c("all_refannot", "expr_refannot", "all_newannot", "expr_newannot"))
gp = gp + theme(axis.text.x=element_text(angle=35, hjust=1, vjust=1))
gp = gp + labs(x= "Transcript set", y='\'$4\'') + ylim(c('$5','$6'))   
gp = gp + scale_fill_manual(name = "Transcript set",
			    labels = c("00_ref" = "all_refannot", "01_ref_expr" = "expr_refannot", "02_string" ="all_newannot", "03_string_expr" ="expr_newannot"),
			    values = c("00_ref" = "#E41A1C", "01_ref_expr" = "#377EB8", "02_string" ="#4DAF4A", "03_string_expr" ="#984EA3"))
ggsave(filename="'$7'")
' | R --vanilla

