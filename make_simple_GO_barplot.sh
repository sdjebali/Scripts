#!/bin/bash
set -Eexo pipefail

# make_simple_GO_barplot.sh
# - accepts a single tsv file with GO enrichment results (see format below, at least GO term, FDR and fold enrichment)
#   as well as a title
# and
# - produces a file with the same name but ending with png instead of tsv that has horizontal barplots
#   with the complete GO term name followed by the FDR on the x axis (or y since turned) and then
#   the bars representing the fold enrichment the bottom one being the smaller and the top one the bigger

# this is to work from Matteo Serino's GO enrichment program output file
# that has 9 columns and is formatted like this (here GO_analysis_for_Sarah_tableno2_miR-877-3p.tsv)
# GO.Term	Mus_musculus_-_REFLIST_(21983)	upload_1_(4439)	upload_1_(expected)	upload_1_(over/under)	upload_1_(fold_Enrichment)	upload_1_(raw_P-value)	upload_1_(FDR)
# lipid_metabolic_process	1227	327	247.77	+	1.32	7.36E-06	5.49E-04
# 10 (8 fields)

# Has been tested on genobioinfo with R 4.3.1 and fragencode renv config

# Example of usage on the above file 
# srun --x11 --mem=8G --time=10:00:00 --pty bash
# cd ~/fragencode/tools/multi/Renv
# module load statistics/R/4.3.1
# pgm=/work/project/fragencode/tools/multi/Scripts/make_simple_GO_barplot.sh
# input=~/bridge/workspace/sdjebali/origamir/papers/1st.2026/GO_analysis_for_Sarah_tableno2_miR-877-3p.tsv
# time $pgm $input miR-877-3p
# real	0m2.726s
# and has produced
# /home/sdjebali/bridge/workspace/sdjebali/origamir/papers/1st.2026/GO_analysis_for_Sarah_tableno2_miR-877-3p.png miR-877-3p


if [ ! -n "$1" ]
then
    echo "" >&2
    echo "Usage: make_simple_GO_barplot.sh file.tsv <title>" >&2
    echo "" >&2
    exit 1
else
    input=$1
    output=${input%.tsv}.png
    pal="#FEB24C"
fi

if [ ! -n "$2" ]
then
    title=""
else
    title=$2
fi

echo '
     library(reshape2)
     library(dplyr)
     library(ggplot2)
     theme_set(theme_bw(base_size = 16))
     # a. read the input file into a data frame and then make a tibble from it
     #########################################################################
     df = read.delim("'$input'", sep="\t", h=TRUE)
     dft=as_tibble(df)
     # b. prepare the label for the GO terms (add pvalues) and define an order (sort by category, then enrichment)
     #############################################################################################################
     dft = dft %>% arrange(upload_1_.fold_Enrichment.) %>%
     mutate(Term.label=paste0(GO.Term, "   ", signif(upload_1_.FDR., 2)),
     Term.label=factor(Term.label, levels=unique(Term.label)))
     # c. make the plot using Jean s code with all terms
     ###################################################
     gp = ggplot(dft, aes(x=Term.label, y=upload_1_.fold_Enrichment.)) +
     geom_col(color="'$pal'", fill="'$pal'") + theme_bw() + coord_flip() +
     ylab("Fold Enrichment") + xlab("GO term (FDR < 0.05)") + ggtitle("'$title'")
     gp = gp + theme(axis.text.x=element_text(hjust=1, vjust=1), axis.text=element_text(size=100), axis.title=element_text(size=120), plot.title=element_text(size=120))
     ggsave(filename="'$output'", h=45, w=45) 
' | R --vanilla
