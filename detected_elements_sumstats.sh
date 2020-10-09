#!/bin/bash
set -Eexo pipefail

# detected_elements_sumstats.sh
###############################
# takes as input:
#################
# - a reference gene annotation file in gtf (such as the ones provided by ensembl) or gff2 format 
# - a predicted gene annotation file in gtf or gff2 format
# - a matrix of transcript TPM values in a set of samples (tsv file with header whose columns are transcript_id, gene_id and TPM values in a set of samples)
#   this file should contain transcript ids from the reference and from the predicted gene annotation
# - a matrix of gene TPM values in a set of samples (tsv file with header whose columns are gene_id and TPM values in a set of samples)
# outputs 
##########
# - a tsv file with header that has:
####################################
#   * 3 rows: one for the header, one for transcripts and one for genes
#   * 9 columns: one for the element type, one for the number of elements in the reference, the number of elements output by stringtie and its % of the total,
#     then the number and % of total of the ref elements that have tpm 0.1 in at least 2 samples, the number of novel (not in ref) elements output by stringtie,
#     and the number and % of total novel elements that have tpm 0.1 in at least 2 samples
# - the reference and predicted transcript gff file filtered for TPM above 0.1 in at least 2 samples"
# - a Plots and a Tables directory with plots and tables about ratios, sizes and positions of transcripts from 4 different sets"
#   (ref, ref_expr, string, string_expr) with respect to the reference annotation"
# - some .err files about the different processing steps"
# - a tsv file with header called make_summary_stat_from_4annots.tsv that has:
###############################################################################
#   * 5 rows, one for the header, one for the reference gene annotation, one for the ref gene annotation where only
#     transcripts with tpm above 0.1 in at least 2 samples were retained, one for the predicted gene annotation and one for the pred gene annotation+
#     where only transcripts with tpm above 0.1 in at least 2 samples were retained
#   * 14 columns with number and ratios of elements in each of the 4 gene prediction sets




# check all the inputs are there
################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ]
then
    echo "" >&2
    echo Usage: detected_elements_sumstats.sh ref_annot.gtf pred_annot.gff tr_TPM.tsv gn_TPM.tsv >&2
    echo "" >&2
    echo where >&2
    echo "- ref_annot.gtf is a gtf (or gff2) file of the reference gene annotation with at least exons rows and with gene_id and transcript_id in the 9th field" >&2
    echo "- pred_annot.gff is a gff2 file with at least exons rows and with gene_id and transcript_id in the 9th field" >&2
    echo "- tr_TPM.tsv is a tsv file with header that has transcript id, gene id and then TPM expression in a set of samples" >&2
    echo "- gn_TPM.tsv is a tsv file with header that has gene id and then TPM expression in a set of samples" >&2
    echo "and will write in current working directory" >&2
    echo "- a tsv file with header with number and ratio of elements in 4 gene annotations (ref annot, ref annot with transcripts whose TPM is at least 0.1" >&2
    echo "  in at least 2 samples, stringtie annot and stringtie annot with transcript whose TPM is at least 0.1 in at least 2 samples)" >&2
    echo "- a tsv file with header that has the number of initial reference transcripts and genes" >&2
    echo "  as well as those in the expression matrices and those with expression above 0.1 in at least 2 samples" >&2
    echo "- the reference and predicted transcript gff file filtered for TPM above 0.1 in at least 2 samples" >&2
    echo "- a Plots and a Tables directory with plots and tables about ratios, sizes and positions of transcripts from 4 different sets" >&2
    echo "  (ref, ref_expr, string, string_expr) with respect to the reference annotation" >&2
    echo "- some .err files about the different processing steps" >&2
    echo "- a tsv file with header with number and ratio of elements in 4 gene annotations (ref annot, ref annot with transcripts whose TPM is at least 0.1" >&2
    echo "  in at least 2 samples, stringtie annot and stringtie annot with transcript whose TPM is at least 0.1 in at least 2 samples)" >&2
    echo "Needs bedtools and R to be installed (tested with versionns 2.29.0 and 3.6.2 resp) and to have ggplot2 library available in R" >&2
    exit 1
fi


# Set general variables
#######################
execpath=`readlink -f $0`
rootdir=`dirname $execpath`
# path="`dirname \"$0\"`" # relative path
# rootdir="`( cd \"$path\" && pwd )`" # absolute path
ref=$1
str=$2
trexpr=$3
gnexpr=$4

# Set variables of programs
###########################
TRGNNB=$rootdir/compute_detected_tr_gn.sh
TRPOS=$rootdir/ref_pred_all_expr_2_tr_ratio_size_pos.sh 
ELTNB=$rootdir/ref_pred_annot2stats.sh


# 1. Call compute_detected_tr_gn.sh and redirect main output tsv file in a file called detected_transcripts_genes_numbers.tsv
##############################################################################################################################
echo "I am calling compute_detected_tr_gn.sh" >&2
$TRGNNB $ref $trexpr $gnexpr > detected_transcripts_genes_numbers.tsv
echo "done" >&2

# 2. Call ref_pred_all_expr_2_tr_ratio_size_pos.sh
##################################################
echo "I am calling ref_pred_all_expr_2_tr_ratio_size_pos.sh" >&2
$TRPOS $ref $str $trexpr $PWD
echo "done" >&2

# 3. Call ref_pred_annot2stats.sh and redirect main output tsv file in a file called detected_elements_numbers.tsv
###################################################################################################################
echo "I am calling ref_pred_annot2stats.sh" >&2
$ELTNB $ref $str $trexpr > detected_elements_numbers.tsv
echo "done" >&2
