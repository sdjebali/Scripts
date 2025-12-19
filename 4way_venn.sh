#!/bin/bash
set -Eexo pipefail

# 4way_venn.sh accepts 6 inputs
# - file for 1st set of the venn (string)
# - file for 2nd set of the venn (string)
# - file for 3rd set of the venn (string)
# - file for 4th set of the venn (string)
# - title for the venn (string)
# - file for venn plot (string)
# !!! on Dec 17th 2025 has been changed to accept tsv files with several columns of which the 1st is of interest !!!
# !!! and also to have ids for the 4 sets so a total of 10 inputs !!!
# TODO: extract the legend as well

# has been tested on genologin with R 3.3.3 and fragencode configuration
# module use /work/project/fragencode/.local/privatemodules/
# module load fragencode-2017.10.06

# pgm=/work/project/fragencode/tools/multi/Scripts/4way_venn.sh
# set1=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/connected.components/sniffles.minimap.father.INV.cc.of.run1.txt 
# set2=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/connected.components/sniffles.minimap.father.INV.cc.of.run2.txt 
# set3=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/connected.components/sniffles.minimap.father.INV.cc.of.run3.txt 
# setall=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/connected.components/sniffles.minimap.father.INV.cc.of.runall.txt
# output=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/venn.diagrams/sniffles.minimap.father.INV.cc.3runs.and.merged.venn.png
# time $pgm $set1 $set2 $set3 $setall "trio1-father-minimap+sniffles-INV" $output run1 run2 run3 runall


echo '
     library("futile.logger")
     library("VennDiagram")
     png("'$6'")
     run1=unlist(read.delim("'$1'", sep="\t", h=FALSE)[,1])
     run2=unlist(read.delim("'$2'", sep="\t", h=FALSE)[,1])
     run3=unlist(read.delim("'$3'", sep="\t", h=FALSE)[,1])
     runall=unlist(read.delim("'$4'", sep="\t", h=FALSE)[,1])
     venn.plot <- venn.diagram(
     x=list( 
     A=run1,
     B=run2,
     C=run3,
     D=runall
),
	main="'$5'",   
	filename = NULL,
	fill=c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4"), 
	alpha=0.50,
	category = c("'$7'","'$8'","'$9'","'${10}'"), #legend
	fontfamily="arial",
	main.cex=2,
	main.fontfamily="arial",
	sub.cex=1.3,
	sub.fontfamily="arial",
	cat.fontfamily=c("arial","arial","arial","arial") 
	)
	grid.draw(venn.plot)
	dev.off()
' | R --vanilla
