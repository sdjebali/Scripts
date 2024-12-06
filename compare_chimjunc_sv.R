## !!! DISCLAIMER !!!
## This is a custom script to parse VCF made by Manta. VCF made by other
## SV detection tools might break the assumptions made here and might require
## some adjustment when parsing the variants

## To install dependencies
##
## If needed, install the Bioconductor package manager
## > install.packages('BiocManager')
## Then, install the relevant packages
## > BiocManager::install(c("VariantAnnotation", "GenomicRanges"))
##
## If needed, to install dplyr:
## > install.packages('dplyr')

args = commandArgs(TRUE)

if(length(args) < 3){
  stop('Usage: Rscript compare_chimjunc_sv.R chimjunc.SAMP.tsv SAMPT_vs_SAMPC.manta.diploid_sv.vcf.gz output.tsv 10000')
}

## arguments for debugging
## args = c('valid.laura.gaelle.hg38.sveachrow.tsv', 'LMS28T_vs_LMS28C.manta.somatic_sv.vcf.gz', 'test.tsv')

## default max distance between breapoints
max.ol.dist = 10000
if(length(args) > 3){
  max.ol.dist = as.integer(args[4])
}
  
## default without using strand
use.strand = FALSE
if(length(args) > 4){
  use.strand = as.logical(args[5])
}
  
## load packages
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))

##
## read chimeras
##
chim.df = read.table(args[1], as.is=TRUE, header=TRUE)

## add an "ID" to keep track of breakpoints and chimeras
chim.df$id = 1:nrow(chim.df)

## parse the first chimera breakpoint
chim.junc1.df = chim.df %>% mutate(seqnames=gsub('(.+)_.+_.:.+_.+_.', '\\1', junc_id),
                                   start=gsub('.+_(.+)_.:.+_.+_.', '\\1', junc_id),
                                   start=as.integer(start),
                                   end=start+1, bkpt=1,
                                   strand=gsub('.+_.+_(.):.+_.+_.', '\\1', junc_id)) %>%
  select(seqnames, start, end, strand, id, bkpt)

## parse the second chimera breakpoint
chim.junc2.df = chim.df %>% mutate(seqnames=gsub('.+_.+_.:(.+)_.+_.', '\\1', junc_id),
                                   start=gsub('.+_.+_.:.+_(.+)_.', '\\1', junc_id),
                                   start=as.integer(start),
                                   end=start+1, bkpt=2,
                                   strand=gsub('.+_.+_.:.+_.+_(.)', '\\1', junc_id)) %>%
  select(seqnames, start, end, strand, id, bkpt)

## merge into the same Genomic Ranges object
chim.junc.gr = rbind(chim.junc1.df, chim.junc2.df) %>% 
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

## chim.junc.gr %>% as.data.frame %>%
##   arrange(id, bkpt) %>% 
##   group_by(id) %>% filter(length(unique(seqnames))==1, all(strand=='+'), diff(start)<0) %>% head

##
## read SVs
##
sv.vcf = readVcf(args[2])

## add an "ID" to keep track of breakpoints and SVs
rowRanges(sv.vcf)$id = 1:length(sv.vcf)

## to save the Genomic Ranges object for SV "junctions"
sv.junc.gr = GRanges()

## prepare breakpoints for BND (translocation-like signal)
sv.bnd = subset(sv.vcf, SVTYPE=='BND')
updateBNDjunc <- function(al, bnd.df, pattern,
                          bp, bt, sp, st){
  parse.o = gsub(pattern, '\\1', al)
  bnd.df %>% mutate(pos_t=ifelse(parse.o!=al, parse.o, pos_t),
                    bkpt_p=ifelse(parse.o!=al, bp, bkpt_p),
                    bkpt_t=ifelse(parse.o!=al, bt, bkpt_t),
                    strand_p=ifelse(parse.o!=al, sp, strand_p),
                    strand_t=ifelse(parse.o!=al, st, strand_t))
}
if(length(sv.bnd) > 0){
  al = unlist(lapply(rowRanges(sv.bnd)$ALT, '[', 1))
  bnd.df = tibble(pos_t=rep(NA, length(al))) %>%
    mutate(bkpt_p=NA, bkpt_t=NA, strand_p=NA, strand_t=NA)
  ## p[t[ like chimera p=1+ -> t=2+
  bnd.df = updateBNDjunc(al, bnd.df,
                         pattern='.+\\[(.*:.*)\\[',
                         bp=1, bt=2, sp='+', st='+')
  ## ]t]p like chimera t=1+ -> p=2+
  bnd.df = updateBNDjunc(al, bnd.df,
                         pattern='\\](.*:.*)\\].+',
                         bp=2, bt=1, sp='+', st='+')
  ## p]t] like chimera p=1+ -> t=2-
  bnd.df = updateBNDjunc(al, bnd.df,
                         pattern='.+\\](.*:.*)\\]',
                         bp=1, bt=2, sp='+', st='-')
  ## [t[p like chimera p=2+ -> t=1-
  bnd.df = updateBNDjunc(al, bnd.df,
                         pattern='\\[(.*:.*)\\[.+',
                         bp=2, bt=1, sp='+', st='-')
  ## check if some didn't get parsed
  if(any(is.na(bnd.df$pos_t))){
    stop('Error: some BND SVs were not parsed.')
  }
  ## make Genomic Ranges object
  bnd.p.gr = rowRanges(sv.bnd)
  names(bnd.p.gr) = mcols(bnd.p.gr) = NULL
  bnd.p.gr$id = rowRanges(sv.bnd)$id
  bnd.p.gr$bkpt = bnd.df$bkpt_p
  strand(bnd.p.gr) = bnd.df$strand_p
  bnd.t.gr = GRanges(bnd.df$pos_t,
                     strand=bnd.df$strand_t,
                     id=rowRanges(sv.bnd)$id,
                     bkpt=bnd.df$bkpt_t)
  sv.junc.gr = c(sv.junc.gr, bnd.p.gr, bnd.t.gr)
}

## prepare breakpoints for DEL(etions)
sv.del = subset(sv.vcf, SVTYPE=='DEL')
if(length(sv.del) > 0){
  del1.gr = rowRanges(sv.del)
  names(del1.gr) = mcols(del1.gr) = NULL
  del1.gr$id = rowRanges(sv.del)$id
  del1.gr$bkpt = 1
  strand(del1.gr) = '+'
  ## second breakpoint in the same "strand" but further
  del2.gr = del1.gr
  del2.gr$bkpt = 2
  end(del2.gr) = info(sv.del)$END
  start(del2.gr) = info(sv.del)$END
  sv.junc.gr = c(sv.junc.gr, del1.gr, del2.gr)
}

## prepare breakpoints for DUP(lications)
sv.dup = subset(sv.vcf, SVTYPE=='DUP')
if(length(sv.dup) > 0){
  dup2.gr = rowRanges(sv.dup)
  ## check that all the duplication calls are predicted to
  ## be tandemly-duplicated sequences. Should be but just in case
  if(any(unlist(dup2.gr$ALT) != '<DUP:TANDEM>')){
    stop('Error: some duplications are not tandem-duplication. Improve script.')
  }
  names(dup2.gr) = mcols(dup2.gr) = NULL
  dup2.gr$id = rowRanges(sv.dup)$id
  dup2.gr$bkpt = 2
  strand(dup2.gr) = '+'
  ## second position is the 1st breakpoint in the same "strand"
  dup1.gr = dup2.gr
  dup1.gr$bkpt = 1
  end(dup1.gr) = info(sv.dup)$END
  start(dup1.gr) = info(sv.dup)$END
  sv.junc.gr = c(sv.junc.gr, dup1.gr, dup2.gr)
}

## we don't care about insertions (right?)

## invert all the SV junctions to consider also the other orientation
sv.junc.gr$orient = 1
sv.junc.inv = sv.junc.gr
sv.junc.inv$orient = 2
strand(sv.junc.inv) = ifelse(strand(sv.junc.inv)=='+', '-', '+')
sv.junc.inv$bkpt = ifelse(sv.junc.inv$bkpt==1, 2, 1)
sv.junc.gr = c(sv.junc.gr, sv.junc.inv)

## TEMPORARY hack: don't consider strands
if (!use.strand){
  strand(sv.junc.gr) = "+"
  strand(chim.junc.gr) = "+"
}

## overlap everything and propagate information to match later
ol = findOverlaps(chim.junc.gr, sv.junc.gr, maxgap=max.ol.dist) %>%
  as.data.frame %>%
  mutate(chim.bkpt=chim.junc.gr$bkpt[queryHits],
         chim.id=chim.junc.gr$id[queryHits],
         sv.bkpt=sv.junc.gr$bkpt[subjectHits],
         sv.orient=sv.junc.gr$orient[subjectHits],
         sv.id=sv.junc.gr$id[subjectHits],
         dist.bkpt=distance(chim.junc.gr[queryHits], sv.junc.gr[subjectHits]))
## head(ol)

## find cases where bkpt configuration match
match.ol = ol %>% group_by(chim.id, sv.id, sv.orient) %>%
  filter(all(chim.bkpt==sv.bkpt), length(unique(chim.bkpt)) > 1)

## if any matches found, add the information about distance
if(nrow(match.ol) == 0){
  chim.df %>% select(-id) %>% mutate(svid=NA, maxdistsv=NA) %>% 
    write.table(args[3], sep='\t', quote=FALSE, row.names=FALSE)
} else {
  
  match.ol = match.ol %>% group_by(chim.id, sv.id, sv.orient) %>%
    summarize(maxdistsv=max(dist.bkpt), .groups='drop')

  ## head(match.ol)

  ## add some info from the chimeras
  match.ol = chim.df %>% mutate(chim.id=id) %>% select(chim.id, junc_id) %>%
    merge(match.ol)

  ## add some info from the SVs
  sv.alts = unlist(lapply(rowRanges(sv.vcf)$ALT, '[', 1))
  sv.df = tibble(sv.id=rowRanges(sv.vcf)$id,
                 sv.label=paste(as.character(rowRanges(sv.vcf)),
                                ifelse(grepl('<', sv.alts), 
                                       abs(info(sv.vcf)$END),
                                       sv.alts),
                                info(sv.vcf)$SVTYPE,
                                sep="_")) %>%
    mutate(sv.label=make.names(sv.label))
  match.ol = merge(match.ol, sv.df)

  ## combine into one (potentially list) of SVs per chimera
  match.ol = match.ol %>% group_by(junc_id) %>%
    summarize(svid=paste(sv.label, collapse=','),
              maxdistsv=paste(maxdistsv, collapse=','))

  chim.df %>% select(-id) %>% merge(match.ol, all.x=TRUE) %>%
    write.table(args[3], sep='\t', quote=FALSE, row.names=FALSE)
}

## NEXT: take CIPOS/CIEND into account
