
##################
# OPTION PARSING #
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-k", "--key"), default="ENTREZID", help="type of gene identifiers used in other arguments (one of \"ENTREZID\", \"ENSEMBL\", \"REFSEQ\", \"GENENAME\", \"SYMBOL\") [default=%default]"),
  make_option(c("-u", "--universe"), default="None", help="list of background gene identifiers (of same type as --key), WITHOUT header.
              Use \"None\" for default background (given by the clusterProfiler package) [default=%default]"),
  make_option(c("-G", "--genes"), default="stdin",
              help="list of foreground gene identifiers (of same type as --key), WITHOUT header [default=%default]"),
  make_option(c("-c", "--ontology"), help="choose the GO category < BP | MF | CC > [default=%default]", default="BP"),
  make_option(c("-p", "--pval"), help="p-value threshold. All GO terms found to be enriched up to this threshold will be saved as output. [default=%default]", default=0.1),
  make_option(c("-q", "--qval"), help="q-value threshold. All GO terms found to be enriched up to this threshold will be saved as output. [default=%default]", default=0.2),
  make_option(c("-a", "--padjMethod"), help="p-value adjustement method threshold. One of \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\", \"none\" [default=%default]", default="fdr"),
  make_option(c("-f", "--fdr"), help="FDR (or adjusted p-value in general) threshold. All GO terms found to be enriched up to this threshold will be saved as output in a separate table. [default=%default]", default=0.1),
  make_option(c("-v", "--verbose"), default=TRUE, help="Verbosity (TRUE/FALSE) [default=%default]"),
  make_option(c("-o", "--output"), default="out", help="additional tags for otuput [default=%default]"),
  make_option(c("-d", "--output_dir"), default="./GO_output/", help="directory for the output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
verbose = opt$verbose

#################################
# Load libraries & Import data #
#################################

if(verbose){
  library("clusterProfiler")
  library("org.Hs.eg.db")
  library(enrichplot)
  library(ggplot2)
} else{
  suppressPackageStartupMessages(library("clusterProfiler"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))
  suppressPackageStartupMessages(library(enrichplot))
  suppressPackageStartupMessages(library(ggplot2))
}

if (verbose) print("Loading input data...")
key = opt$key
if(opt$universe=="None"){
  universe = NA
  print("Warning: using defaut universe automatically provided by the clusterProfiler package")
} else{
  U = read.table(opt$universe, h=F, col.names='hs')
  universe = as.character(unique(U$hs))
}

G = read.table(opt$genes, h=F, col.names='hs')
genes = unique(G$hs)
ontology = opt$ontology
pval = opt$pval
qval = opt$qval
pAdjMethod = opt$padjMethod
fdr_threshold = opt$fdr

if (verbose) print("Done.")

######################
# Compute enrichment #
######################

to_readable = TRUE
if (key=="SYMBOL") to_readable = FALSE

if (verbose) print("Computing GO enrichment...")

res = enrichGO(gene = genes,
               OrgDb = org.Hs.eg.db,
               keyType = key,
               ont = ontology,
               universe = universe,
               pvalueCutoff = pval,
               qvalueCutoff  = qval,
               pAdjustMethod = pAdjMethod,
               readable = to_readable)

if (verbose) print("Done.")

universe_found = res@universe

###################
# Compute results #
###################

if(opt$universe=="None"){
  sprintf("%s (default) background genes", length(res@universe))
} else{
  sprintf("%s provided background genes; %s with a corresponding %s id", nrow(U), length(universe), key)
}
n_genes_found = strsplit(res@result$GeneRatio[[1]],'/')[[1]][2]
sprintf("%s provided genes; %s found by `enrichGO`", nrow(G), n_genes_found)
sprintf("Computed GO enrichment (whether significant or not) for %s distinct GO terms", nrow(res@result))
sprintf("Of those %s GO terms, %s have a %s-adjusted p-val < %s", nrow(res@result), nrow(res@result[res@result$p.adjust<=fdr_threshold,]), pAdjMethod, fdr_threshold)

################
# Write output #
################

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
output = sprintf("%s/%s.%s.%s%s", opt$output_dir, opt$output, opt$ontology, opt$fdr, opt$padjMethod)

if (verbose) print("Writing outputs tables...")

# Full table
write.table(res@result, file=sprintf("%s.tsv", output), quote=F, sep="\t", row.names=F)

# Most significant GO terms
write.table(res@result[res@result$p.adjust<=fdr_threshold,], file=sprintf("%s.significant.tsv", output), quote=F, sep="\t", row.names=F)

if (verbose) print("Done. Writing output images...")

png(file=sprintf("%s.upset.png", output),
    width=600, height=350)
upsetplot(res) + ggtitle("UpSet Plot. Number of genes per significant GO term")
graphics.off()

png(file=sprintf("%s.dotplot.png", output),
    width=600, height=350)
dotplot(res) + ggtitle(sprintf("p.adjust based on %s correction", pAdjMethod))
graphics.off()

png(file=sprintf("%s.gene-concept.png", output),
    width=600, height=350)
cnetplot(res, categorySize="pvalue", foldChange=genes, colorEdge = TRUE) + ggtitle(sprintf("Gene-%s Network, circular", ontology))
graphics.off()

png(file=sprintf("%s.gene-concept.circular.png", output),
    width=600, height=350)
cnetplot(res, categorySize="pvalue", foldChange=genes, circular = TRUE, colorEdge = TRUE) + ggtitle(sprintf("Gene-%s Network", ontology))
graphics.off()

png(file=sprintf("%s.heatplot.png", output),
    width=600, height=350)
heatplot(res) + ggtitle("Heatplot: significant GO terms against contributing genes")
graphics.off()

pbar <- function(x) {
  pi=seq(0, 1, length.out=11)
  
  mutate(x, pp = cut(p.adjust, pi)) %>%
    group_by(pp) %>% 
    summarise(cnt = n()) %>% 
    ggplot(aes(pp, cnt)) + geom_col() + 
    theme_minimal() +
    xlab("p value intervals") +
    ylab("Frequency")
}   

if (verbose) print("Writing last output image (this one might take some time)...")

random_genes = sample(universe_found, as.integer(n_genes_found))
enrich_ref = enrichGO(gene = random_genes,
                      OrgDb = org.Hs.eg.db,
                      keyType = key,
                      ont = ontology,
                      universe = universe,
                      pvalueCutoff = pval,
                      qvalueCutoff  = qval,
                      pAdjustMethod = pAdjMethod,
                      readable = to_readable)

p1 <- pbar(res) + ggtitle("p value distribution")
p2 <- pbar(enrich_ref) + ggtitle("p value distribution (one single random sample of genes)")
png(file=sprintf("%s.pvalues.png", output),
    width=600, height=350)
cowplot::plot_grid(p1, p2, ncol=1, labels = LETTERS[1:2])
graphics.off()

if (verbose) print("Done.")

q(save='no')
